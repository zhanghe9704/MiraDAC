/**
 * @file davector.cpp
 * @brief Implementation of DAVector<T> methods that cannot be inline:
 *   - print() / operator<<  (need print_vec logic)
 *   - Explicit instantiation for double.
 */

#include "da/da.h"

#include <iomanip>
#include <iostream>
#include <cmath>
#include <limits>
#include <cstring>
#ifdef DA_WITH_SYMBOLIC
#  include <symengine/expression.h>
#  include <symengine/printers.h>
#  include <sstream>
#endif

namespace da {

// ===========================================================================
// Internal print helper (mirrors print_vec in tpsa_extend.cc)
// ===========================================================================

static void print_vec_impl(DAEnv* env, unsigned slot, std::ostream& os)
{
    Layout& layout = env->layout();
    Pool<double>& pool = env->pool<double>();

    const unsigned gnv = layout.num_vars();
    const unsigned gnd = layout.max_order();
    const unsigned full_len = layout.full_len();
    const unsigned* bptr = layout.base();

    const double* v = pool.slot(slot);
    const unsigned l = pool.len(slot);

    int width_base = 2;
    if (gnd > 9) ++width_base;

    int cnt_width = 1;
    if (l > 9) cnt_width = static_cast<int>(std::ceil(std::log10(static_cast<double>(l))));
    ++cnt_width;

    std::string start(static_cast<std::size_t>(cnt_width), ' ');
    std::string sep(static_cast<std::size_t>(cnt_width), '-');
    // 'I' marker at the end of start
    start[start.size() - 1] = 'I';

    std::ios::fmtflags prevflags = os.flags();

    os << start;
    os << "          V [" << slot << "]              Base  [ "
       << l << " / " << full_len << " ]" << std::endl
       << sep << "----------------------------------------------" << std::endl;

    int cnt = 0;
    const unsigned* p = bptr;
    for (unsigned i = 0; i < l; ++i) {
        if (std::abs(v[i]) < std::numeric_limits<double>::min()) {
            p += gnv;
            continue;
        }
        ++cnt;
        os << std::setw(cnt_width) << cnt;
        os << ' ' << std::setprecision(15) << std::scientific
           << std::setw(15 + 8) << v[i] << "    ";
        for (unsigned j = 0; j < gnv - 1; ++j) {
            os << std::setw(width_base) << static_cast<unsigned>(*p - *(p + 1));
            ++p;
        }
        os << std::setw(width_base) << static_cast<unsigned>(*p++) << std::setw(6) << i
           << std::endl;
    }
    os << std::endl;
    os.flags(prevflags);
}

// Two-component (complex DA) print helper
static void print_vec2_impl(DAEnv* env, unsigned slot_r, unsigned slot_i, std::ostream& os)
{
    Layout& layout = env->layout();
    Pool<double>& pool = env->pool<double>();

    const unsigned gnv = layout.num_vars();
    const unsigned gnd = layout.max_order();
    const unsigned full_len = layout.full_len();
    const unsigned* bptr = layout.base();

    const double* vr = pool.slot(slot_r);
    const double* vi = pool.slot(slot_i);

    unsigned lr = pool.len(slot_r);
    unsigned li = pool.len(slot_i);
    unsigned l = std::max(lr, li);

    int width_base = 2;
    if (gnd > 9) ++width_base;

    int cnt_width = 1;
    if (l > 9) cnt_width = static_cast<int>(std::ceil(std::log10(static_cast<double>(l))));
    ++cnt_width;

    std::string start(static_cast<std::size_t>(cnt_width), ' ');
    std::string sep(static_cast<std::size_t>(cnt_width), '-');
    start[start.size() - 1] = 'I';

    std::ios::fmtflags prevflags = os.flags();

    os << start;
    os << "          V [" << slot_r << "/" << slot_i << "]              Base  [ "
       << l << " / " << full_len << " ]" << std::endl
       << sep << "----------------------------------------------" << std::endl;

    int cnt = 0;
    const unsigned* p = bptr;
    for (unsigned i = 0; i < l; ++i) {
        double absv = std::max(
            (i < lr ? std::abs(vr[i]) : 0.0),
            (i < li ? std::abs(vi[i]) : 0.0));
        if (absv < std::numeric_limits<double>::min()) {
            p += gnv;
            continue;
        }
        ++cnt;
        os << std::setw(cnt_width) << cnt;
        double rv = (i < lr ? vr[i] : 0.0);
        double iv = (i < li ? vi[i] : 0.0);
        os << ' ' << std::setprecision(15) << std::scientific
           << std::setw(15 + 8) << rv << ' '
           << std::setw(15 + 8) << iv << "    ";
        for (unsigned j = 0; j < gnv - 1; ++j) {
            os << std::setw(width_base) << static_cast<unsigned>(*p - *(p + 1));
            ++p;
        }
        os << std::setw(width_base) << static_cast<unsigned>(*p++) << std::setw(6) << i
           << std::endl;
    }
    os << std::endl;
    os.flags(prevflags);
}

// ===========================================================================
// DAVector<double>::print
// ===========================================================================

template<>
void DAVector<double>::print() const
{
    print_vec_impl(env_, slot_, std::cout);
}

// ===========================================================================
// operator<<
// ===========================================================================

template<>
std::ostream& operator<<(std::ostream& os, const DAVector<double>& v)
{
    print_vec_impl(v.env_, v.slot_, os);
    return os;
}

template<>
std::ostream& operator<<(std::ostream& os, const std::complex<DAVector<double>>& cd)
{
    const DAVector<double>& r = get_real(cd);
    const DAVector<double>& im = get_imag(cd);
    print_vec2_impl(r.env_, r.slot_, im.slot_, os);
    return os;
}

// ===========================================================================
// Explicit instantiation
// ===========================================================================

#ifdef DA_WITH_SYMBOLIC
// operator<< for SDA (DAVector<SymEngine::Expression>)
template<>
std::ostream& operator<<(std::ostream& os, const DAVector<SymEngine::Expression>& v)
{
    using E = SymEngine::Expression;
    DAEnv* env = v.env_;
    Layout& layout = env->layout();
    Pool<E>& pool = env->pool<E>();

    const unsigned gnv = layout.num_vars();
    const unsigned gnd = layout.max_order();
    const unsigned full_len = layout.full_len();
    const unsigned* bptr = layout.base();

    const E* vals = pool.slot(v.slot_);
    const unsigned l = pool.len(v.slot_);

    int width_base = 2;
    if (gnd > 9) ++width_base;

    int cnt_width = 1;
    if (l > 9) cnt_width = static_cast<int>(std::ceil(std::log10(static_cast<double>(l))));
    ++cnt_width;

    std::string start(static_cast<std::size_t>(cnt_width), ' ');
    std::string sep(static_cast<std::size_t>(cnt_width), '-');
    start[start.size() - 1] = 'I';

    os << start;
    os << "          V [" << v.slot_ << "]              Base  [ "
       << l << " / " << full_len << " ]" << std::endl
       << sep << "----------------------------------------------" << std::endl;

    int cnt = 0;
    const unsigned* p = bptr;
    for (unsigned i = 0; i < l; ++i) {
        // Skip zero terms
        bool is_z = da::is_zero(vals[i]);
        if (is_z) {
            p += gnv;
            continue;
        }
        ++cnt;
        std::ostringstream ss;
        ss << vals[i];
        os << std::setw(cnt_width) << cnt;
        os << ' ' << std::setw(32) << ss.str() << "    ";
        for (unsigned j = 0; j < gnv - 1; ++j) {
            os << std::setw(width_base) << static_cast<unsigned>(*p - *(p + 1));
            ++p;
        }
        os << std::setw(width_base) << static_cast<unsigned>(*p++) << std::setw(6) << i
           << std::endl;
    }
    os << std::endl;
    return os;
}

// operator<< for CSDA: the two parts one after the other (symbolic
// coefficients are too wide for the side-by-side layout of CNDA).
template<>
std::ostream& operator<<(std::ostream& os, const std::complex<DAVector<SymEngine::Expression>>& cd)
{
    os << "Real part:" << std::endl << get_real(cd)
       << "Imaginary part:" << std::endl << get_imag(cd);
    return os;
}

// DAVector<SymEngine::Expression>::print — reuse the SDA operator<< above.
// Declared before the explicit instantiation of DAVector<Expression> below so
// the specialization is picked up (otherwise print() would be undefined for SDA).
template<>
void DAVector<SymEngine::Expression>::print() const
{
    std::cout << *this;
}
#endif // DA_WITH_SYMBOLIC

template struct DAVector<double>;

template DAVector<double> operator+(const DAVector<double>&, const DAVector<double>&);
template DAVector<double> operator-(const DAVector<double>&, const DAVector<double>&);
template DAVector<double> operator*(const DAVector<double>&, const DAVector<double>&);
template DAVector<double> operator/(const DAVector<double>&, const DAVector<double>&);
template DAVector<double> operator+(const DAVector<double>&, double);
template DAVector<double> operator+(double, const DAVector<double>&);
template DAVector<double> operator-(const DAVector<double>&, double);
template DAVector<double> operator-(double, const DAVector<double>&);
template DAVector<double> operator*(const DAVector<double>&, double);
template DAVector<double> operator*(double, const DAVector<double>&);
template DAVector<double> operator/(const DAVector<double>&, double);
template DAVector<double> operator/(double, const DAVector<double>&);
template DAVector<double> operator+(const DAVector<double>&);
template DAVector<double> operator-(const DAVector<double>&);
template bool operator==(const DAVector<double>&, const DAVector<double>&);
template bool same_env(const DAVector<double>&, const DAVector<double>&) noexcept;

#ifdef DA_WITH_SYMBOLIC
template struct DAVector<SymEngine::Expression>;

template DAVector<SymEngine::Expression> operator+(const DAVector<SymEngine::Expression>&, const DAVector<SymEngine::Expression>&);
template DAVector<SymEngine::Expression> operator-(const DAVector<SymEngine::Expression>&, const DAVector<SymEngine::Expression>&);
template DAVector<SymEngine::Expression> operator*(const DAVector<SymEngine::Expression>&, const DAVector<SymEngine::Expression>&);
template DAVector<SymEngine::Expression> operator/(const DAVector<SymEngine::Expression>&, const DAVector<SymEngine::Expression>&);
template DAVector<SymEngine::Expression> operator+(const DAVector<SymEngine::Expression>&, double);
template DAVector<SymEngine::Expression> operator+(double, const DAVector<SymEngine::Expression>&);
template DAVector<SymEngine::Expression> operator-(const DAVector<SymEngine::Expression>&, double);
template DAVector<SymEngine::Expression> operator-(double, const DAVector<SymEngine::Expression>&);
template DAVector<SymEngine::Expression> operator*(const DAVector<SymEngine::Expression>&, double);
template DAVector<SymEngine::Expression> operator*(double, const DAVector<SymEngine::Expression>&);
template DAVector<SymEngine::Expression> operator/(const DAVector<SymEngine::Expression>&, double);
template DAVector<SymEngine::Expression> operator/(double, const DAVector<SymEngine::Expression>&);
template DAVector<SymEngine::Expression> operator+(const DAVector<SymEngine::Expression>&);
template DAVector<SymEngine::Expression> operator-(const DAVector<SymEngine::Expression>&);
template bool operator==(const DAVector<SymEngine::Expression>&, const DAVector<SymEngine::Expression>&);
template bool same_env(const DAVector<SymEngine::Expression>&, const DAVector<SymEngine::Expression>&) noexcept;
#endif // DA_WITH_SYMBOLIC

} // namespace da
