/**
 * @file main.cc
 * @brief Consumer smoke test: what a downstream find_package(da) user sees.
 *
 * Regression guard for compile definitions that used to be applied with
 * add_compile_definitions() (build-tree only) instead of being PUBLIC usage
 * requirements on the library targets. A consumer then compiled the headers
 * with DA_CHECK_ENV undefined — the preprocessor reads an undefined
 * identifier as 0, so the cross-environment guard vanished without a word —
 * and, against a symbolic install, without DA_WITH_SYMBOLIC, giving a DAEnv
 * that lacks pool_e_ and disagrees with the library about sizeof(DAEnv).
 */
#include "da/da.h"

#include <cstdio>
#include <exception>

static int failures = 0;

static void check(bool ok, const char* what) {
    printf("%-46s %s\n", what, ok ? "ok" : "FAILED");
    if (!ok) ++failures;
}

int main() {
    printf("sizeof(da::DAEnv) = %zu\n", sizeof(da::DAEnv));

#ifdef DA_CHECK_ENV
    check(DA_CHECK_ENV == 1, "DA_CHECK_ENV reaches the consumer as 1");
#else
    check(false, "DA_CHECK_ENV is defined for the consumer");
#endif

    da::da_init(4, 2, 500);
    da::DAEnv* main_env = &da::da_current_env();

    // The guard must fire on an operation mixing two environments.
    {
        da::DAEnv& env2 = da::da_make_env(3, 2, 200);
        bool threw = false;
        {
            da::NDA v2;                         // belongs to env2
            try {
                da::NDA bad = v2 + da::base[0]; // base[] belongs to main env
                (void)bad;
            } catch (const std::exception&) {
                threw = true;
            }
        }
        check(threw, "cross-environment operation throws");
        da::da_select_env(*main_env);
        delete &env2;
    }

    // A symbolic install must expose the symbolic half to consumers too.
#ifdef DA_WITH_SYMBOLIC
    printf("DA_WITH_SYMBOLIC: defined (symbolic install)\n");
    {
        da::SDA s;                              // must compile and run
        check(s.con() == SymEngine::Expression(0), "SDA usable in consumer code");
    }
#else
    printf("DA_WITH_SYMBOLIC: not defined (numeric-only install)\n");
#endif

    da::da_clear();

    printf(failures ? "\nconsumer smoke FAILED\n" : "\nconsumer smoke passed\n");
    return failures ? 1 : 0;
}
