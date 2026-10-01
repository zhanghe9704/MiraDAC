/**
 * @file test_pool.cc
 * @brief Tests for Pool<T> — both the POD path (double) and the non-POD path (std::string).
 */

#include "catch.hpp"
#include "da/pool.h"

#include <string>
#include <stdexcept>
#include <cstring>  // memcmp

// ======================================================================
// Helpers
// ======================================================================

// A custom non-trivial type that tracks copies/moves to help detect
// memset/memcpy misuse on the non-POD path.
struct Tracker {
    std::string value;
    Tracker() : value("") {}
    explicit Tracker(const std::string& s) : value(s) {}
    // Non-trivially copyable (has a non-trivial copy ctor via std::string)
};

static_assert(!std::is_trivially_copyable_v<std::string>,
              "std::string must be non-trivially-copyable for the non-POD path test");
static_assert( std::is_trivially_copyable_v<double>,
              "double must be trivially-copyable for the POD path test");

// ======================================================================
// Pool<double>  — POD path
// ======================================================================

TEST_CASE("Pool<double> basic reserve and assign", "[pool][double]") {
    da::Pool<double> p;
    p.reserve(35, 100);

    REQUIRE(p.poolsize() == 100);
    REQUIRE(p.full_len() == 35);
    REQUIRE(p.remain()   == 100);
    REQUIRE(p.count()    == 0);

    // Assign 100 slots
    std::vector<unsigned> ids;
    for (int i = 0; i < 100; ++i) {
        REQUIRE_NOTHROW([&]{ ids.push_back(p.assign()); }());
    }
    REQUIRE(p.remain() == 0);
    REQUIRE(p.count()  == 100);

    // 101st assign must throw
    REQUIRE_THROWS_AS(p.assign(), std::runtime_error);
}

TEST_CASE("Pool<double> alloc sets len=1 and zeros", "[pool][double]") {
    da::Pool<double> p;
    p.reserve(35, 10);

    unsigned idx = p.alloc();
    REQUIRE(p.len(idx) == 1);

    // Check slot is zeroed
    const double* s = p.slot(idx);
    for (unsigned k = 0; k < 35; ++k)
        REQUIRE(s[k] == 0.0);
}

TEST_CASE("Pool<double> free then assign round-trip", "[pool][double]") {
    da::Pool<double> p;
    p.reserve(35, 5);

    unsigned a = p.assign();
    unsigned b = p.assign();
    // Write something to slot a
    p.slot(a)[0] = 42.0;
    p.set_len(a, 5);

    // Free a — should zero it and return to pool
    p.free(a);
    REQUIRE(p.remain() == 4);

    // Assign again — we get some free slot back (the reference returns the
    // most-recently-freed slot from the tail of the free-list).
    unsigned c = p.assign();
    REQUIRE(p.len(c) == 0);  // assign sets len = 0
    (void)b; (void)c;
}

TEST_CASE("Pool<double> zero_slot and copy_slot", "[pool][double]") {
    da::Pool<double> p;
    p.reserve(10, 4);

    unsigned s0 = p.assign();
    unsigned s1 = p.assign();

    // Fill s0 with known values
    for (unsigned k = 0; k < 10; ++k)
        p.slot(s0)[k] = static_cast<double>(k + 1);
    p.set_len(s0, 10);

    // copy_slot: s0 -> s1 (5 elements)
    p.copy_slot(p.slot(s0), p.slot(s1), 5);
    for (unsigned k = 0; k < 5; ++k)
        REQUIRE(p.slot(s1)[k] == p.slot(s0)[k]);

    // zero_slot: s0
    p.zero_slot(p.slot(s0));
    for (unsigned k = 0; k < 10; ++k)
        REQUIRE(p.slot(s0)[k] == 0.0);
}

TEST_CASE("Pool<double> reset zeros slot and clears len", "[pool][double]") {
    da::Pool<double> p;
    p.reserve(10, 4);

    unsigned idx = p.alloc();
    for (unsigned k = 0; k < 10; ++k)
        p.slot(idx)[k] = 1.0;
    p.set_len(idx, 10);

    p.reset(idx);
    REQUIRE(p.len(idx) == 0);
    for (unsigned k = 0; k < 10; ++k)
        REQUIRE(p.slot(idx)[k] == 0.0);
}

TEST_CASE("Pool<double> move constructor", "[pool][double]") {
    da::Pool<double> p1;
    p1.reserve(35, 50);

    unsigned idx = p1.alloc();
    p1.slot(idx)[0] = 99.0;

    // Move-construct p2 from p1
    da::Pool<double> p2(std::move(p1));

    // p2 is usable
    REQUIRE(p2.poolsize() == 50);
    REQUIRE(p2.full_len() == 35);
    REQUIRE(p2.slot(idx)[0] == 99.0);

    // p1 is in a safe empty state — assigning to it should be fine
    // (but it has no capacity)
    p1.reserve(5, 5);
    REQUIRE(p1.poolsize() == 5);
}

TEST_CASE("Pool<double> move assignment", "[pool][double]") {
    da::Pool<double> p1;
    p1.reserve(35, 50);
    unsigned idx = p1.alloc();
    p1.slot(idx)[1] = 7.0;

    da::Pool<double> p2;
    p2.reserve(10, 10);

    p2 = std::move(p1);
    REQUIRE(p2.poolsize() == 50);
    REQUIRE(p2.slot(idx)[1] == 7.0);
    // p1 is empty now; can safely be destructed
}

TEST_CASE("Pool<double> pool_clean resets from idx", "[pool][double]") {
    da::Pool<double> p;
    p.reserve(35, 10);

    // Assign first 3 slots and write data
    for (int i = 0; i < 3; ++i) {
        unsigned idx = p.assign();
        p.slot(idx)[0] = static_cast<double>(idx + 1);
        p.set_len(idx, 1);
    }
    REQUIRE(p.remain() == 7);

    // Clean from index 1
    p.pool_clean(1);

    // After clean from 1, slots 1..9 are free, slot 0 is untouched
    REQUIRE(p.remain() == 9);
}

// ======================================================================
// Pool<std::string>  — non-POD path
// ======================================================================

TEST_CASE("Pool<string> basic reserve and assign", "[pool][string]") {
    da::Pool<std::string> p;
    p.reserve(8, 100);

    REQUIRE(p.poolsize() == 100);
    REQUIRE(p.remain()   == 100);

    // Assign all 100
    for (int i = 0; i < 100; ++i)
        p.assign();

    REQUIRE(p.remain() == 0);
    REQUIRE_THROWS_AS(p.assign(), std::runtime_error);
}

TEST_CASE("Pool<string> zero_slot assigns T{} not memset", "[pool][string]") {
    da::Pool<std::string> p;
    p.reserve(4, 5);

    unsigned idx = p.assign();
    // Fill slot with non-empty strings
    for (unsigned k = 0; k < 4; ++k)
        p.slot(idx)[k] = "hello_" + std::to_string(k);

    // zero_slot must call ~std::string + construct empty strings,
    // not blindly memset (which would corrupt refcount internals).
    p.zero_slot(p.slot(idx));

    for (unsigned k = 0; k < 4; ++k) {
        // After zero_slot, each element should be a default-constructed string
        REQUIRE(p.slot(idx)[k] == std::string{});
    }
}

TEST_CASE("Pool<string> copy_slot survives round-trip", "[pool][string]") {
    da::Pool<std::string> p;
    p.reserve(4, 5);

    unsigned src = p.assign();
    unsigned dst = p.assign();

    // Fill src with known values
    for (unsigned k = 0; k < 4; ++k)
        p.slot(src)[k] = "item_" + std::to_string(k);

    // copy_slot must use operator= not memcpy
    p.copy_slot(p.slot(src), p.slot(dst), 4);

    for (unsigned k = 0; k < 4; ++k)
        REQUIRE(p.slot(dst)[k] == p.slot(src)[k]);

    // Now zero the source and confirm dst is independent
    p.zero_slot(p.slot(src));
    for (unsigned k = 0; k < 4; ++k)
        REQUIRE(p.slot(dst)[k] == ("item_" + std::to_string(k)));
}

TEST_CASE("Pool<string> free then assign resets string slot", "[pool][string]") {
    da::Pool<std::string> p;
    p.reserve(4, 5);

    unsigned idx = p.assign();
    p.slot(idx)[0] = "persistent";
    p.set_len(idx, 1);

    p.free(idx);
    REQUIRE(p.remain() == 5);  // all 5 back

    // After free, the slot was zeroed (strings are empty again)
    for (unsigned k = 0; k < 4; ++k)
        REQUIRE(p.slot(idx)[k] == std::string{});
}

TEST_CASE("Pool<string> move constructor — no double free", "[pool][string]") {
    da::Pool<std::string> p1;
    p1.reserve(4, 20);

    unsigned idx = p1.alloc();
    p1.slot(idx)[0] = "moved";

    da::Pool<std::string> p2(std::move(p1));
    REQUIRE(p2.slot(idx)[0] == "moved");

    // p1 should be empty and safe to destroy
    REQUIRE(p1.poolsize() == 0);
    // No crash on destruction of either
}

TEST_CASE("Pool<string> assign len is 0; alloc len is 1", "[pool][string]") {
    da::Pool<std::string> p;
    p.reserve(4, 5);

    unsigned a = p.assign();
    REQUIRE(p.len(a) == 0);

    unsigned b = p.alloc();
    REQUIRE(p.len(b) == 1);
}

// ===========================================================================
// A slot freed while the pool is exhausted must go back on the free list.
// ===========================================================================
TEST_CASE("Pool<double> free after exhaustion recycles the slot", "[pool][double]") {
    da::Pool<double> p;
    p.reserve(35, 3);
    unsigned a = p.assign(), b = p.assign(), c = p.assign();
    REQUIRE(p.remain() == 0);
    REQUIRE_THROWS_AS(p.assign(), std::runtime_error);

    p.free(b);
    REQUIRE(p.remain() == 1);
    REQUIRE(p.assign() == b);
    REQUIRE(p.remain() == 0);

    p.free(a);
    p.free(b);
    p.free(c);
    REQUIRE(p.remain() == 3);
    REQUIRE(p.count() == 0);
    unsigned x = p.assign(), y = p.assign(), z = p.assign();
    REQUIRE(x != y);
    REQUIRE(y != z);
    REQUIRE(x != z);
}

TEST_CASE("Pool exhaustion throws da::PoolExhausted", "[pool][double]") {
    da::Pool<double> p;
    p.reserve(35, 2);
    p.assign();
    p.assign();
    REQUIRE_THROWS_AS(p.assign(), da::PoolExhausted);
    REQUIRE_THROWS_AS(p.alloc(), std::runtime_error);
    REQUIRE_THROWS_WITH(p.assign(), "Pool::assign: Run out of vectors");
}

// ===========================================================================
// A slot from assign()/alloc() is always zero, whichever data it held before
// it was freed. Doubles are zeroed when the slot is handed out (it is about to
// be written, so the cache is warm); non-trivial types are cleared when the
// slot is freed, so their resources (e.g. SymEngine references) go at once.
// ===========================================================================
TEST_CASE("Pool<double> reused slot comes back zero", "[pool][double]") {
    da::Pool<double> p;
    p.reserve(8, 2);
    unsigned a = p.alloc();
    for (unsigned k = 0; k < 8; ++k) p.slot(a)[k] = 1.0 + k;
    p.set_len(a, 8);
    p.free(a);
    unsigned b = p.assign();
    unsigned c = p.assign();
    for (unsigned s : {b, c})
        for (unsigned k = 0; k < 8; ++k) REQUIRE(p.slot(s)[k] == 0.0);
    REQUIRE(p.len(b) == 0);
    p.free(b);
    unsigned d = p.alloc();
    for (unsigned k = 0; k < 8; ++k) REQUIRE(p.slot(d)[k] == 0.0);
    REQUIRE(p.len(d) == 1);
}

TEST_CASE("Pool<string> free releases the elements at once", "[pool][string]") {
    da::Pool<std::string> p;
    p.reserve(4, 2);
    unsigned a = p.alloc();
    p.slot(a)[2] = std::string(1000, 'x');
    p.free(a);
    for (unsigned k = 0; k < 4; ++k) REQUIRE(p.slot(a)[k].empty());
    unsigned b = p.alloc();
    for (unsigned k = 0; k < 4; ++k) REQUIRE(p.slot(b)[k].empty());
}

// ===========================================================================
// The most recently freed slot is handed out first (a stack), so a kernel's
// temporaries keep reusing a few warm slots whatever order the caller freed
// its vectors in (a garbage-collected caller frees them in scattered order).
// ===========================================================================
TEST_CASE("Pool hands out the most recently freed slot first", "[pool][double]") {
    da::Pool<double> p;
    p.reserve(8, 6);
    unsigned a = p.assign(), b = p.assign(), c = p.assign();
    p.free(a);
    p.free(c);
    REQUIRE(p.assign() == c);
    REQUIRE(p.assign() == a);
    p.free(b);
    REQUIRE(p.assign() == b);
    REQUIRE(p.remain() == 3);
    REQUIRE(p.count() == 3);
}
