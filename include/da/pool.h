/**
 * @file pool.h
 * @brief Templated RAII memory pool with O(1) free-list.
 *
 * @details Port of the ad_reserve / ad_assign / ad_alloc / ad_free /
 *   ad_pool_clean / ad_reset mechanics from
 *   ref/tpsa/src/tpsa_extend.cc, turned into a type-safe RAII class.
 *
 * The pool owns one contiguous block of T[full_len * n], divided into
 * n slots of full_len elements each.
 *
 * Free-list structure (ported from ref/tpsa/src/tpsa_extend.cc):
 *   - free_[i]  : "next slot after i" (mirrors adlist[i]).
 *   - head_     : index of the next slot to hand out (mirrors ad_flag).
 *   - tail_     : index of the last slot in the free-list (mirrors ad_end).
 *   - free_[tail_] holds the sentinel value (= size_).
 *   - Exhaustion: head_ == free_[tail_].
 *
 * POD fast-path vs non-POD path (if constexpr):
 *   - If T is trivially copyable, zero_slot uses memset and
 *     copy_slot uses memcpy.
 *   - Otherwise, per-element T{} assignment / operator= is used.
 *     memset/memcpy are NEVER used on the non-POD path to avoid
 *     corrupting reference counts (e.g. SymEngine::Expression).
 *
 * @note assign() throws std::runtime_error when exhausted (the
 *   reference printed "Run out of vectors" and called exit(-1)).
 */

#pragma once

#include <vector>
#include <cstring>
#include <stdexcept>
#include <type_traits>
#include <cmath>     // std::abs

namespace da {

template<class T>
class Pool {
public:
    // ------------------------------------------------------------------ //
    //  Construction / destruction                                          //
    // ------------------------------------------------------------------ //

    /// Default-constructed pool is empty; call reserve() before use.
    Pool() noexcept
        : block_(nullptr), head_(0), tail_(0), full_len_(0), size_(0)
    {}

    ~Pool() noexcept { destroy(); }

    // Non-copyable
    Pool(const Pool&)            = delete;
    Pool& operator=(const Pool&) = delete;

    // Movable
    Pool(Pool&& other) noexcept { move_from(other); }
    Pool& operator=(Pool&& other) noexcept {
        if (this != &other) {
            destroy();
            move_from(other);
        }
        return *this;
    }

    // ------------------------------------------------------------------ //
    //  Allocation / free-list management                                   //
    // ------------------------------------------------------------------ //

    /**
     * @brief Allocate one contiguous block of full_len*n T's, wire slots,
     *        build the free-list, zero all lengths.
     *
     * If the pool was previously reserved, the old block is released first
     * (matching the ad_reserve "if advecpool != nullptr, delete it" logic).
     *
     * @param full_len  Number of T elements per slot.
     * @param n         Number of slots.
     */
    void reserve(unsigned full_len, unsigned n) {
        if (n == 0) return;

        // Release previous allocation if any
        destroy();

        full_len_ = full_len;
        size_     = n;

        // One contiguous block; value-initialize to T{} (zeros doubles,
        // default-constructs non-POD types properly).
        block_ = new T[static_cast<std::size_t>(full_len) * n]();

        slot_.resize(n);
        len_.resize(n, 0u);

        // free_ has n entries (indices 0..n-1), each holding the "next" index.
        // Mirrors adlist in the reference (which has exactly n entries).
        // free_[i] = i+1 for i in [0, n-2]; free_[n-1] = n  (sentinel value = n)
        free_.resize(n);
        for (unsigned i = 0; i < n - 1; ++i)
            free_[i] = i + 1;
        free_[n - 1] = n;  // tail's "next" = sentinel value n

        // Wire slot pointers into the contiguous block
        for (unsigned i = 0; i < n; ++i)
            slot_[i] = block_ + static_cast<std::size_t>(i) * full_len;

        // head_ = 0 (ad_flag), tail_ = n-1 (ad_end)
        head_ = 0;
        tail_ = n - 1;
    }

    /**
     * @brief Pop the next free slot.  Returns its index.  Sets len_[i] = 0.
     *
     * Throws std::runtime_error on exhaustion.
     * (Reference: "Run out of vectors" + exit(-1).)
     */
    unsigned assign() {
        // Exhaustion: head_ == free_[tail_]  (mirrors: ad_flag == adlist[ad_end])
        if (head_ == free_[tail_]) {
            throw std::runtime_error("Pool::assign: Run out of vectors");
        }
        unsigned i = head_;
        len_[i]    = 0;
        head_      = free_[i];  // advance head (mirrors: ad_flag = adlist[ad_flag])
        return i;
    }

    /**
     * @brief Like assign(), but also zeros the slot and sets len = 1.
     * (Port of ad_alloc semantics.)
     */
    unsigned alloc() {
        unsigned i = assign();
        zero_slot(slot_[i]);
        len_[i] = 1;
        return i;
    }

    /**
     * @brief Return slot i to the free-list (O(1)).
     * (Port of ad_free — appends to the tail.)
     *
     * Zeros the slot and resets its length before recycling.
     */
    void free(unsigned i) {
        if (size_ == 0) return;  // pool already destroyed; no-op
        reset(i);
        // Append i to the tail of the free-list.
        // Before: free_[tail_] = sentinel; adlist[ad_end] = n
        // After:  free_[i] = old sentinel, free_[tail_] = i, tail_ = i
        free_[i]     = free_[tail_];  // i's next = old sentinel (= size_)
        free_[tail_]  = i;            // old tail now points to i
        tail_         = i;            // i is the new tail
    }

    // ------------------------------------------------------------------ //
    //  Slot management                                                     //
    // ------------------------------------------------------------------ //

    /**
     * @brief Zero-fill slot i and set its length to 0.
     * (Port of ad_reset.)
     */
    void reset(unsigned i) {
        zero_slot(slot_[i]);
        len_[i] = 0;
    }

    /**
     * @brief Remove coefficients smaller than |eps| in slot i.
     * (Port of ad_clean.)
     */
    void clean(unsigned i, double eps) {
        T* p = slot_[i];
        unsigned N = 0;
        double abseps = std::abs(eps);
        for (unsigned k = 0; k < len_[i]; ++k) {
            if (std::abs(static_cast<double>(p[k])) < abseps)
                p[k] = T{};
            else
                N = k;
        }
        if (len_[i] > N + 1)
            len_[i] = N + 1;
    }

    /**
     * @brief Release all slots from index idx onward (zero + rebuild free-list).
     * (Port of ad_pool_clean(idx).)
     */
    void pool_clean(unsigned idx) {
        unsigned n = size_;
        // Zero data for slots idx..n-1
        for (unsigned k = idx; k < n; ++k) {
            zero_slot(slot_[k]);
            len_[k] = 0;
        }
        // Rebuild free-list from idx..n-1
        for (unsigned k = idx; k < n - 1; ++k)
            free_[k] = k + 1;
        free_[n - 1] = n;  // sentinel
        head_ = idx;
        tail_ = n - 1;
    }

    // ------------------------------------------------------------------ //
    //  Slot primitives — POD vs non-POD dispatch via if constexpr         //
    // ------------------------------------------------------------------ //

    /// Zero all full_len_ elements in the slot pointed to by p.
    void zero_slot(T* p) noexcept {
        if constexpr (std::is_trivially_copyable_v<T>) {
            std::memset(p, 0, static_cast<std::size_t>(full_len_) * sizeof(T));
        } else {
            for (unsigned k = 0; k < full_len_; ++k)
                p[k] = T{};
        }
    }

    /// Copy len elements from src to dst.
    void copy_slot(const T* src, T* dst, unsigned len) noexcept {
        if constexpr (std::is_trivially_copyable_v<T>) {
            std::memcpy(dst, src, static_cast<std::size_t>(len) * sizeof(T));
        } else {
            for (unsigned k = 0; k < len; ++k)
                dst[k] = src[k];
        }
    }

    // ------------------------------------------------------------------ //
    //  Accessors                                                           //
    // ------------------------------------------------------------------ //

    T*       slot(unsigned i)       { return slot_[i]; }
    const T* slot(unsigned i) const { return slot_[i]; }

    unsigned len(unsigned i)        const { return len_[i]; }
    void     set_len(unsigned i, unsigned l) { len_[i] = l; }

    unsigned poolsize() const { return size_; }
    unsigned full_len() const { return full_len_; }

    /**
     * @brief Number of slots currently in use (total - remaining).
     */
    unsigned count() const { return size_ - remain(); }

    /**
     * @brief Number of free slots remaining (walks the free-list).
     */
    unsigned remain() const {
        unsigned cnt  = 0;
        unsigned cur  = head_;
        unsigned sentinel = free_[tail_];
        while (cur != sentinel) {
            ++cnt;
            cur = free_[cur];
        }
        return cnt;
    }

private:
    // ------------------------------------------------------------------ //
    //  Internal helpers                                                    //
    // ------------------------------------------------------------------ //

    void destroy() noexcept {
        // Zero size_ FIRST so any lingering Pool::free() calls see size_==0
        // and return early (no-op) — this is important for the pattern where
        // da_clear() is called before local DAVectors go out of scope.
        size_    = 0;
        full_len_ = 0;
        delete[] block_;
        block_ = nullptr;
        slot_.clear();
        len_.clear();
        free_.clear();
        head_ = tail_ = 0;
    }

    void move_from(Pool& other) noexcept {
        block_    = other.block_;
        slot_     = std::move(other.slot_);
        len_      = std::move(other.len_);
        free_     = std::move(other.free_);
        head_     = other.head_;
        tail_     = other.tail_;
        full_len_ = other.full_len_;
        size_     = other.size_;

        // Leave source in a safe, empty state
        other.block_    = nullptr;
        other.head_     = 0;
        other.tail_     = 0;
        other.full_len_ = 0;
        other.size_     = 0;
    }

    // ------------------------------------------------------------------ //
    //  Members                                                             //
    // ------------------------------------------------------------------ //

    T*                    block_;    ///< Contiguous allocation (owned)
    std::vector<T*>       slot_;     ///< slot_[i] = block_ + i*full_len_
    std::vector<unsigned> len_;      ///< Current used length of each slot
    std::vector<unsigned> free_;     ///< Linked-list "next" pointers (size n)
    unsigned              head_;     ///< Next slot to assign (mirrors ad_flag)
    unsigned              tail_;     ///< Last slot in free-list (mirrors ad_end)
    unsigned              full_len_; ///< Elements per slot
    unsigned              size_;     ///< Total number of slots
};

} // namespace da
