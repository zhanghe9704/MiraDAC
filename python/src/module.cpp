#include "common.h"

NB_MODULE(_core, m) {
    m.attr("__version__") = DA_VERSION;
#ifdef DA_WITH_SYMBOLIC
    m.attr("HAS_SYMBOLIC") = true;
#else
    m.attr("HAS_SYMBOLIC") = false;
#endif
    // Overloads taking NDAList, CNDAList, ... try to convert a plain list to
    // each in turn; a failed try is normal overload resolution, not an error.
    nb::set_implicit_cast_warnings(false);
    register_exceptions(m);
    bind_env(m);
    auto nda = bind_nda(m);
    bind_cnda(m, nda);
#ifdef DA_WITH_SYMBOLIC
    bind_expr(m);
    bind_sda(m, nda);
    bind_csda(m);
#endif
}
