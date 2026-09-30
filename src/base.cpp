/**
 * @file base.cpp
 * @brief Base struct, da::base global, cd_composition, inv_map, file I/O.
 *
 * @details Ported from ref/tpsa/src/da.cc.  Algorithms unchanged.
 */

#include "da/da.h"

#include <algorithm>
#include <cassert>
#include <cmath>
#include <fstream>
#include <limits>
#include <set>
#include <sstream>
#include <vector>

namespace da {

// ===========================================================================
// Global bases object
// ===========================================================================

Base base; ///< da::base[i] — the i-th DA base variable.

// ===========================================================================
// da_init / da_clear  (defined here so they can touch da::base)
// ===========================================================================

int da_init(unsigned order, unsigned num_vars, unsigned pool_size, bool table) {
    // Clear base vectors BEFORE destroying the old env so their destructors
    // can safely return slots to the still-live old pool.
    base.base_vecs.clear();
    da_init_env_only(order, num_vars, pool_size, table);
    base.set_base(num_vars);
    if (table) da_current_env().layout().generate_order_table();
    return 0;
}

void da_clear() {
    base.base_vecs.clear();
    da_clear_env_only();
}

// ===========================================================================
// Base implementation
// ===========================================================================

void Base::set_base(unsigned int n) {
    // base_vecs should already be empty (cleared in da_init before env swap),
    // but clear defensively to avoid any stale slots.
    base_vecs.clear();
    base_vecs.resize(n);   // default-constructs n NDA objects from current (new) env
    for (unsigned int i = 0; i < n; ++i) {
        DAEnv& env = da_current_env();
        detail::ad_var(env.layout(), env.pool<double>(), base_vecs[i].slot_, 0.0, i);
    }
}

NDA da_base(unsigned int i) {
    NDA v;
    detail::ad_var(v.env_->layout(), v.env_->pool<double>(), v.slot_, 0.0, i);
    return v;
}

void Base::set_base() {
    set_base(static_cast<unsigned int>(NDA::dim()));
}

// ===========================================================================
// LU decomposition helpers for inv_map (ported from da.cc verbatim)
// ===========================================================================

static void _ludcmp(std::vector<std::vector<double>>& a, const int n,
                    std::vector<int>& idx, int& d)
{
    int imax = 0;
    double big, dum, sum, temp;
    std::vector<double> vv(n);
    auto tiny = std::numeric_limits<double>::min();
    d = 1;
    for (int i = 0; i < n; ++i) {
        big = 0.0;
        for (int j = 0; j < n; ++j)
            if ((temp = std::fabs(a[i][j])) > big) big = temp;
        if (big < tiny)
            throw std::invalid_argument("inv_map: the linear part of the map is singular");
        vv[i] = 1.0 / big;
    }
    for (int j = 0; j < n; ++j) {
        for (int i = 0; i < j; ++i) {
            sum = a[i][j];
            for (int k = 0; k < i; ++k) sum -= a[i][k] * a[k][j];
            a[i][j] = sum;
        }
        big = 0.0;
        for (int i = j; i < n; ++i) {
            sum = a[i][j];
            for (int k = 0; k < j; ++k) sum -= a[i][k] * a[k][j];
            a[i][j] = sum;
            if ((dum = vv[i] * std::fabs(sum)) >= big) {
                big = dum;
                imax = i;
            }
        }
        if (j != imax) {
            for (int k = 0; k < n; ++k) { dum = a[imax][k]; a[imax][k] = a[j][k]; a[j][k] = dum; }
            d = -d;
            vv[imax] = vv[j];
        }
        idx[j] = imax;
        // The reference replaced a zero pivot with tiny, which turned a
        // singular map into huge or inf coefficients without any error.
        if (std::fabs(a[j][j]) < tiny)
            throw std::invalid_argument("inv_map: the linear part of the map is singular");
        if (j != (n - 1)) {
            dum = 1.0 / (a[j][j]);
            for (int i = j + 1; i < n; ++i) a[i][j] *= dum;
        }
    }
}

static void _lubksb(std::vector<std::vector<double>>& a, const int n,
                    std::vector<int>& idx, std::vector<double>& b)
{
    int ii = -1, ip;
    double sum;
    for (int i = 0; i < n; ++i) {
        ip = idx[i]; sum = b[ip]; b[ip] = b[i];
        if (ii + 1)
            for (int j = ii; j <= i - 1; ++j) sum -= a[i][j] * b[j];
        else if (sum) ii = i;
        b[i] = sum;
    }
    for (int i = n - 1; i >= 0; --i) {
        sum = b[i];
        for (int j = i + 1; j < n; ++j) sum -= a[i][j] * b[j];
        b[i] = sum / a[i][i];
    }
}

static void _inv_matrix(std::vector<std::vector<double>>& a, const int n,
                        std::vector<std::vector<double>>& y)
{
    int d = 0;
    std::vector<double> col(n);
    std::vector<int> idx(n);
    _ludcmp(a, n, idx, d);
    for (int j = 0; j < n; ++j) {
        for (int i = 0; i < n; ++i) col[i] = 0;
        col[j] = 1;
        _lubksb(a, n, idx, col);
        for (int i = 0; i < n; ++i) y[i][j] = col[i];
    }
}

// ===========================================================================
// inv_map
// ===========================================================================

void inv_map(std::vector<NDA>& ivecs, int dim, std::vector<NDA>& ovecs)
{
    assert(dim <= NDA::dim() && "Wrong dimension of map in inv_map!");
    // Note: the reference writes "v.con() < min()", which is satisfied by any
    // negative constant part; the magnitude is what must be zero.
    for (auto& v : ivecs)
        assert(std::fabs(v.con()) < std::numeric_limits<double>::min() &&
               "Constant part of the input map is NOT zero in inv_map!");

    std::vector<std::vector<double>> lin_matrix(dim, std::vector<double>(dim));
    std::vector<int> c(dim, 0);
    for (int i = 0; i < dim; ++i) {
        for (int j = 0; j < dim; ++j) {
            c[j] = 1;
            lin_matrix[i][j] = ivecs[i].element(c);
            c[j] = 0;
        }
    }

    // Restore the caller's order, not the original one: restore_order() would
    // discard an order the caller lowered with da_change_order().
    Layout& layout = da_current_env().layout();
    const unsigned int caller_order = layout.max_order();
    std::vector<NDA> nlin_map;
    for (auto& v : ivecs) {
        layout.change_order(1);
        NDA t = v;
        layout.change_order(caller_order);
        t = v - t;
        nlin_map.push_back(t);
    }

    // Base vectors of the current env (the global da::base belongs to the
    // default env only).
    std::vector<NDA> bases;
    for (int j = 0; j < dim; ++j) bases.push_back(da_base(static_cast<unsigned int>(j)));

    std::vector<std::vector<double>> inv_lin_matrix(dim, std::vector<double>(dim));
    _inv_matrix(lin_matrix, dim, inv_lin_matrix);

    std::vector<NDA> inv_lin_map;
    for (int i = 0; i < dim; ++i) {
        NDA t(0.0);
        for (int j = 0; j < dim; ++j) {
            if (std::fabs(inv_lin_matrix[i][j]) > 1e-16)
                t = t + inv_lin_matrix[i][j] * bases[j];
        }
        inv_lin_map.push_back(t);
    }

    std::vector<NDA> tmp_ovecs;
    for (auto t : inv_lin_map) tmp_ovecs.push_back(t);
    std::vector<NDA> tmp(dim);

    for (int i = 0; i < NDA::order() + 1; ++i) {
        da_composition(nlin_map, tmp_ovecs, tmp);
        for (int j = 0; j < dim; ++j)
            tmp_ovecs[j] = bases[j] - tmp[j];
        da_composition(inv_lin_map, tmp_ovecs, tmp);
        for (int j = 0; j < dim; ++j)
            tmp_ovecs[j] = tmp[j];
    }
    for (int i = 0; i < dim; ++i)
        ovecs[i] = tmp_ovecs[i];
}

// ===========================================================================
// Helper: "is the cached power slot already computed?"
// For NDA: use abs() > 0 (the original check).
// For SDA: use a std::vector<bool> initialized table (passed by reference).
// ===========================================================================

// cdav_norm: returns double norm for T=double (the only path where inited==nullptr).
// This is only used in the T=double instantiation since SDA always passes inited!=nullptr.
static double cdav_norm(const std::complex<NDA>& v) {
    return abs(v);  // da::abs(complex<NDA>) -> double
}

// ===========================================================================
// complex_da_pow_int_pos — template over DAVector<T>
// The memoization table "power_v" is shared; "inited" tracks filled slots for SDA.
// For T=double the original abs-check is preserved (inited not used).
// ===========================================================================

template<class T>
std::complex<DAVector<T>>& complex_da_pow_int_pos(
    const std::complex<DAVector<T>>& iv,
    std::vector<std::complex<DAVector<T>>>& power_v,
    int order, int order_rec,
    std::vector<bool>* inited)
{
    std::size_t key = static_cast<std::size_t>(order * order_rec);
    // Check already-computed:
    // For T=double (inited==nullptr), use abs-check. For symbolic (inited!=nullptr), use table.
    bool already;
    if constexpr (std::is_same<T, double>::value) {
        already = (inited != nullptr)
            ? (*inited)[key]
            : (cdav_norm(power_v.at(key)) > 0);
    } else {
        already = (inited != nullptr) && (*inited)[key];
    }
    if (already) return power_v.at(key);

    if (order == 0) {
        get_real(power_v.at(0)) = DAVector<T>(1.0);
        get_imag(power_v.at(0)) = DAVector<T>(0.0);
        if (inited) (*inited)[0] = true;
        return power_v.at(0);
    } else if (order == 1) {
        std::size_t idx1 = static_cast<std::size_t>(1 * order_rec);
        bool a1;
        if constexpr (std::is_same<T, double>::value) {
            a1 = (inited != nullptr) ? (*inited)[idx1] : (cdav_norm(power_v.at(idx1)) > 0);
        } else {
            a1 = (inited != nullptr) && (*inited)[idx1];
        }
        if (!a1) {
            get_real(power_v.at(idx1)) = get_real(iv);
            get_imag(power_v.at(idx1)) = get_imag(iv);
            if (inited) (*inited)[idx1] = true;
        }
        return power_v.at(idx1);
    } else if (order == 2) {
        std::size_t idx2 = static_cast<std::size_t>(2 * order_rec);
        bool a2;
        if constexpr (std::is_same<T, double>::value) {
            a2 = (inited != nullptr) ? (*inited)[idx2] : (cdav_norm(power_v.at(idx2)) > 0);
        } else {
            a2 = (inited != nullptr) && (*inited)[idx2];
        }
        if (!a2) {
            get_real(power_v.at(idx2)) =
                get_real(iv) * get_real(iv) - get_imag(iv) * get_imag(iv);
            get_imag(power_v.at(idx2)) =
                2.0 * get_real(iv) * get_imag(iv);
            if (inited) (*inited)[idx2] = true;
        }
        return power_v.at(idx2);
    } else if (order & 1) {
        int order_idx = order_rec * 2;
        std::size_t idxE = static_cast<std::size_t>(order_idx);
        bool aE;
        if constexpr (std::is_same<T, double>::value) {
            aE = (inited != nullptr) ? (*inited)[idxE] : (cdav_norm(power_v.at(idxE)) > 0);
        } else {
            aE = (inited != nullptr) && (*inited)[idxE];
        }
        if (!aE) {
            get_real(power_v.at(idxE)) =
                get_real(iv) * get_real(iv) - get_imag(iv) * get_imag(iv);
            get_imag(power_v.at(idxE)) =
                2.0 * get_real(iv) * get_imag(iv);
            if (inited) (*inited)[idxE] = true;
        }
        std::complex<DAVector<T>>& vres = power_v.at(idxE);
        vres = complex_da_pow_int_pos(vres, power_v, order / 2, order_idx, inited);
        std::size_t next_idx = static_cast<std::size_t>(order_rec + (order / 2) * order_idx);
        get_real(power_v.at(next_idx)) =
            get_real(vres) * get_real(iv) - get_imag(vres) * get_imag(iv);
        get_imag(power_v.at(next_idx)) =
            get_real(vres) * get_imag(iv) + get_imag(vres) * get_real(iv);
        if (inited) (*inited)[next_idx] = true;
        return power_v.at(next_idx);
    } else {
        int order_idx = order_rec * 2;
        std::size_t idxE = static_cast<std::size_t>(order_idx);
        bool aE;
        if constexpr (std::is_same<T, double>::value) {
            aE = (inited != nullptr) ? (*inited)[idxE] : (cdav_norm(power_v.at(idxE)) > 0);
        } else {
            aE = (inited != nullptr) && (*inited)[idxE];
        }
        if (!aE) {
            get_real(power_v.at(idxE)) =
                get_real(iv) * get_real(iv) - get_imag(iv) * get_imag(iv);
            get_imag(power_v.at(idxE)) =
                2.0 * get_real(iv) * get_imag(iv);
            if (inited) (*inited)[idxE] = true;
        }
        std::complex<DAVector<T>>& vres = power_v.at(idxE);
        return complex_da_pow_int_pos(vres, power_v, order / 2, order_idx, inited);
    }
}

// ===========================================================================
// get_element_T: access the i-th raw coefficient of a DAVector<T> as type T.
// Unlike DAVector<T>::element(int i) which returns double (calls eval_double),
// this returns the native T — needed for symbolic coefficients.
// ===========================================================================
template<class T>
static const T& get_element_T(const DAVector<T>& v, unsigned i) {
    return v.env_->template pool<T>().slot(v.slot_)[i];
}

// ===========================================================================
// coeff_is_near_zero: compile-time dispatched "skip this coefficient" check.
// For T=double: use std::abs < min.  For T=Expression: use da::is_zero.
// ===========================================================================
static bool coeff_is_near_zero_impl(double x, std::true_type) {
    return std::abs(x) < std::numeric_limits<double>::min();
}

#ifdef DA_WITH_SYMBOLIC
static bool coeff_is_near_zero_impl(const SymEngine::Expression& x, std::false_type) {
    return da::is_zero(x);
}
#endif

template<class T>
static bool coeff_is_near_zero(const T& x) {
    return coeff_is_near_zero_impl(x, std::is_same<T, double>{});
}

// ===========================================================================
// cd_composition — template over DAVector<T>
// Primary overload: vector<DAVector<T>> ivecs composed into complex<DAVector<T>>.
// ===========================================================================

template<class T>
void cd_composition(std::vector<DAVector<T>>& ivecs,
                    std::vector<std::complex<DAVector<T>>>& v,
                    std::vector<std::complex<DAVector<T>>>& ovecs)
{
    int gnv = DAVector<T>::dim();
    int gnd = DAVector<T>::order();
    assert(static_cast<int>(v.size()) == gnv &&
           "Error in cd_composition: No. of DA vectors NOT EQUAL to No. of bases!");
    assert(ivecs.size() == ovecs.size() &&
           "Error in cd_composition: No. of input vectors NOT EQUAL to No. of output vectors!");

    // Initialize power_vv with DAVector<T>(0) for each slot.
    // We cannot use std::complex<double>(0,0) as initializer for T=Expression.
    std::vector<std::vector<std::complex<DAVector<T>>>> power_vv(gnv);
    for (int i = 0; i < gnv; ++i) {
        power_vv[i].resize(gnd + 1); // default-constructs each complex<DAVector<T>>
    }

    // For symbolic, maintain an initialized table (numeric uses abs-check directly)
    const bool need_inited = !std::is_same<T, double>::value;
    std::vector<std::vector<bool>> inited_vv;
    if (need_inited)
        inited_vv.assign(gnv, std::vector<bool>(gnd + 1, false));

    for (int i = 0; i < gnv; ++i) {
        cd_copy(1.0, power_vv[i][0]);
        cd_copy(v[i], power_vv[i][1]);
        if (need_inited) {
            inited_vv[i][0] = true;
            inited_vv[i][1] = true;
        }
    }

    int veclen_max = 0;
    for (auto& vv : ivecs)
        if (static_cast<int>(vv.length()) > veclen_max)
            veclen_max = static_cast<int>(vv.length());

    int vec_size = static_cast<int>(ivecs.size());
    for (int iv = 0; iv < vec_size; ++iv)
        cd_copy(ivecs[iv].con(), ovecs[iv]);

    DAVector<T> product_real, product_imag, tmp;

    if (!da_current_env().layout().has_order_table())
        da_current_env().layout().generate_order_table();

    for (int i = 1; i < veclen_max; ++i) {
        auto& orders = da_element_orders(i);
        bool product_flag = true;
        product_real = DAVector<T>(1.0);
        product_imag = DAVector<T>(0.0);
        bool c_flag = true;
        for (int iv = 0; iv < vec_size; ++iv) {
            if (i >= static_cast<int>(ivecs[iv].length())) continue;
            // Skip near-zero coefficients
            const T& coef_t = get_element_T(ivecs[iv], static_cast<unsigned>(i));
            if (coeff_is_near_zero(coef_t)) continue;
            if (c_flag) {
                for (int id = 0; id < gnv; ++id) {
                    std::vector<bool>* ip = need_inited ? &inited_vv[id] : nullptr;
                    complex_da_pow_int_pos(v[id], power_vv[id], orders[id], 1, ip);
                }
                c_flag = false;
            }
            const T& coef = coef_t;
            if (product_flag) {
                for (int id = 0; id < gnv; ++id) {
                    int ord = orders[id];
                    if (ord > 0) {
                        tmp = product_real * get_real(power_vv[id][ord])
                              - product_imag * get_imag(power_vv[id][ord]);
                        product_imag = product_real * get_imag(power_vv[id][ord])
                                       + product_imag * get_real(power_vv[id][ord]);
                        product_real = tmp;
                    }
                }
                product_flag = false;
            }
            tmp = coef * product_real;
            get_real(ovecs[iv]) += tmp;
            tmp = coef * product_imag;
            get_imag(ovecs[iv]) += tmp;
        }
    }
    for (auto& vv : ovecs) {
        get_real(vv).clean();
        get_imag(vv).clean();
    }
}

template<class T>
void cd_composition(std::vector<std::complex<DAVector<T>>>& ivecs,
                    std::vector<std::complex<DAVector<T>>>& v,
                    std::vector<std::complex<DAVector<T>>>& ovecs)
{
    assert(DAVector<T>::dim() == static_cast<int>(v.size()) &&
           "Error in cd_composition: No. of DA vectors NOT EQUAL to No. of bases!");
    assert(ivecs.size() == ovecs.size() &&
           "Error in cd_composition: No. of input vectors NOT EQUAL to No. of output vectors!");
    int n = static_cast<int>(ivecs.size());
    std::vector<DAVector<T>> iv(2 * n);
    for (int i = 0; i < n; ++i) {
        iv[2 * i]     = get_real(ivecs[i]);
        iv[2 * i + 1] = get_imag(ivecs[i]);
    }
    std::vector<std::complex<DAVector<T>>> ov(2 * n);
    cd_composition(iv, v, ov);
    for (int i = 0; i < n; ++i) {
        get_real(ovecs[i]) = get_real(ov[2 * i])     - get_imag(ov[2 * i + 1]);
        get_imag(ovecs[i]) = get_imag(ov[2 * i])     + get_real(ov[2 * i + 1]);
    }
}

template<class T>
void cd_composition(std::vector<std::complex<DAVector<T>>>& ivecs,
                    std::vector<DAVector<T>>& v,
                    std::vector<std::complex<DAVector<T>>>& ovecs)
{
    assert(DAVector<T>::dim() == static_cast<int>(v.size()) &&
           "Error in cd_composition: No. of DA vectors NOT EQUAL to No. of bases!");
    assert(ivecs.size() == ovecs.size() &&
           "Error in cd_composition: No. of input vectors NOT EQUAL to No. of output vectors!");
    int n = static_cast<int>(ivecs.size());
    std::vector<DAVector<T>> iv(2 * n);
    for (int i = 0; i < n; ++i) {
        iv[2 * i]     = get_real(ivecs[i]);
        iv[2 * i + 1] = get_imag(ivecs[i]);
    }
    std::vector<DAVector<T>> ov(2 * n);
    da_composition(iv, v, ov);
    for (int i = 0; i < n; ++i) {
        get_real(ovecs[i]) = ov[2 * i];
        get_imag(ovecs[i]) = ov[2 * i + 1];
    }
}

// ===========================================================================
// Explicit instantiations for T = double (always)
// ===========================================================================
template std::complex<NDA>& complex_da_pow_int_pos<double>(
    const std::complex<NDA>&, std::vector<std::complex<NDA>>&, int, int, std::vector<bool>*);
template void cd_composition<double>(std::vector<NDA>&,
    std::vector<std::complex<NDA>>&, std::vector<std::complex<NDA>>&);
template void cd_composition<double>(std::vector<std::complex<NDA>>&,
    std::vector<std::complex<NDA>>&, std::vector<std::complex<NDA>>&);
template void cd_composition<double>(std::vector<std::complex<NDA>>&,
    std::vector<NDA>&, std::vector<std::complex<NDA>>&);

// ===========================================================================
// String / file I/O helpers (ported from da.cc)
// ===========================================================================

static std::string ltrim_ws(std::string s) {
    const std::string ws = " \n\r\t\v\f";
    auto st = s.find_first_not_of(ws);
    return (st == std::string::npos) ? "" : s.substr(st);
}
static std::string rtrim_ws(std::string s) {
    const std::string ws = " \n\r\t\v\f";
    auto fi = s.find_last_not_of(ws);
    return (fi == std::string::npos) ? "" : s.substr(0, fi + 1);
}
std::string trim_whitespace(std::string s) {
    return rtrim_ws(ltrim_ws(s));
}

bool read_da_from_file(std::string filename, NDA& d)
{
    d.reset();
    std::fstream input;
    input.open(filename, std::ios::in | std::ios::out | std::ios::app);
    std::string line;
    bool reading = false, read_success = false;
    int n_terms = 0;
    while (std::getline(input, line)) {
        if (!line.empty() && line.back() == '\r') line.pop_back();
        if (!line.empty()) {
            line = trim_whitespace(line);
            if (reading) {
                int idx_term = std::stoi(line.substr(line.length() - 6));
                line.erase(line.end() - 6, line.end());
                double elem = 0;
                std::vector<int> idx;
                std::istringstream iss(line);
                std::string word;
                iss >> word; iss >> word;
                elem = std::stod(word);
                while (!iss.eof()) {
                    iss >> word;
                    if (!word.empty()) idx.push_back(std::stoi(word));
                }
                d.set_element(idx, elem);
                if (n_terms - idx_term == 1) read_success = true;
            } else {
                auto pi = line.find_last_of('/');
                auto pf = line.find_last_of('[');
                if (pi != std::string::npos && pf != std::string::npos) {
                    n_terms = std::stoi(line.substr(pi + 1, pf));
                    continue;
                }
                line.erase(std::remove(line.begin(), line.end(), '-'), line.end());
                if (line.empty()) reading = true;
            }
        }
    }
    return read_success;
}

bool read_cd_from_file(std::string filename, std::complex<NDA>& cd)
{
    NDA& rl = get_real(cd);
    NDA& img = get_imag(cd);
    rl.reset();
    img.reset();
    std::fstream input;
    input.open(filename, std::ios::in | std::ios::out | std::ios::app);
    std::string line;
    bool reading = false, read_success = false;
    int n_terms = 0;
    while (std::getline(input, line)) {
        if (!line.empty() && line.back() == '\r') line.pop_back();
        if (!line.empty()) {
            line = trim_whitespace(line);
            if (reading) {
                int idx_term = std::stoi(line.substr(line.length() - 6));
                line.erase(line.end() - 6, line.end());
                double elem_rl = 0, elem_img = 0;
                std::vector<int> idx;
                std::istringstream iss(line);
                std::string word;
                iss >> word; iss >> word;
                elem_rl = std::stod(word);
                iss >> word;
                elem_img = std::stod(word);
                while (!iss.eof()) {
                    iss >> word;
                    if (!word.empty()) idx.push_back(std::stoi(word));
                }
                if (std::abs(elem_rl) >= std::numeric_limits<double>::min())
                    rl.set_element(idx, elem_rl);
                if (std::abs(elem_img) >= std::numeric_limits<double>::min())
                    img.set_element(idx, elem_img);
                if (n_terms - idx_term == 1) read_success = true;
            } else {
                auto pi = line.find_last_of('/');
                auto pf = line.find_last_of('[');
                if (pi != std::string::npos && pf != std::string::npos) {
                    n_terms = std::stoi(line.substr(pi + 1, pf));
                    continue;
                }
                line.erase(std::remove(line.begin(), line.end(), '-'), line.end());
                if (line.empty()) reading = true;
            }
        }
    }
    return read_success;
}

NDA devide_by_element(NDA& t, NDA& b)
{
    NDA r;
    int l = static_cast<int>(t.length());
    if (l < static_cast<int>(b.length())) l = static_cast<int>(b.length());
    for (int i = 0; i < l; ++i) {
        double te = t.element(i);
        double be = b.element(i);
        if (std::abs(be) > std::numeric_limits<double>::min())
            r.set_element(const_cast<std::vector<int>&>(r.element_orders(i)), te / be);
        else
            r.set_element(const_cast<std::vector<int>&>(r.element_orders(i)), te);
    }
    return r;
}

bool compare_da_vectors(NDA& a, NDA& b, double eps)
{
    NDA r = a - b;
    r = devide_by_element(r, b);
    return r.iszero(eps);
}

bool compare_da_with_file(std::string filename, NDA& d, double eps)
{
    NDA r;
    bool ok = read_da_from_file(filename, r);
    if (!ok) {
        std::cout << "read_da_from_file failed in compare_da_with_file." << std::endl;
        return false;
    }
    return compare_da_vectors(r, d, eps);
}

bool compare_cd_vectors(std::complex<NDA>& a, std::complex<NDA>& b, double eps)
{
    NDA& ar = get_real(a); NDA& ai = get_imag(a);
    NDA& br = get_real(b); NDA& bi = get_imag(b);
    return compare_da_vectors(ar, br, eps) && compare_da_vectors(ai, bi, eps);
}

bool compare_cd_with_file(std::string filename, std::complex<NDA>& d, double eps)
{
    std::complex<NDA> r;
    bool ok = read_cd_from_file(filename, r);
    if (!ok) {
        std::cout << "read_cd_from_file failed in compare_cd_with_file." << std::endl;
        return false;
    }
    return compare_cd_vectors(r, d, eps);
}

#ifdef DA_WITH_SYMBOLIC
// ===========================================================================
// Explicit instantiations for T = SymEngine::Expression (symbolic complex DA)
// ===========================================================================
template std::complex<SDA>& complex_da_pow_int_pos<SymEngine::Expression>(
    const std::complex<SDA>&, std::vector<std::complex<SDA>>&, int, int, std::vector<bool>*);
template void cd_composition<SymEngine::Expression>(std::vector<SDA>&,
    std::vector<std::complex<SDA>>&, std::vector<std::complex<SDA>>&);
template void cd_composition<SymEngine::Expression>(std::vector<std::complex<SDA>>&,
    std::vector<std::complex<SDA>>&, std::vector<std::complex<SDA>>&);
template void cd_composition<SymEngine::Expression>(std::vector<std::complex<SDA>>&,
    std::vector<SDA>&, std::vector<std::complex<SDA>>&);
#endif // DA_WITH_SYMBOLIC

} // namespace da
