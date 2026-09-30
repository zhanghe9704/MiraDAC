/**
 * @file test_env.cc
 * @brief Tests for DAEnv + environment management functions.
 */

#include "catch.hpp"
#include "da/env.h"
#include "da/layout.h"

#include <stdexcept>

// ======================================================================
// Basic construction and current-env
// ======================================================================

TEST_CASE("da_init creates env with correct layout", "[env]") {
    da::da_init(4, 3, 100);
    da::DAEnv& env = da::da_current_env();

    REQUIRE(env.layout().num_vars()  == 3);
    REQUIRE(env.layout().max_order() == 4);
    REQUIRE(env.layout().full_len()  == 35);  // C(3+4,4) = 35

    da::da_clear();
}

TEST_CASE("da_current_env throws before da_init", "[env]") {
    // Ensure no env is selected by clearing first (may already be clear)
    da::da_clear();
    REQUIRE_THROWS_AS(da::da_current_env(), std::runtime_error);
}

TEST_CASE("da_init sets current env; da_current_env returns it", "[env]") {
    da::da_init(4, 3, 100);
    da::DAEnv* p = &da::da_current_env();
    REQUIRE(p != nullptr);

    // Call again: same pointer (same default env)
    REQUIRE(&da::da_current_env() == p);
    da::da_clear();
}

// ======================================================================
// Pool access via env
// ======================================================================

TEST_CASE("da_init reserves double pool with correct size", "[env]") {
    da::da_init(4, 3, 100);
    da::DAEnv& env = da::da_current_env();

    auto& p = env.pool<double>();
    REQUIRE(p.poolsize()  == 100);
    REQUIRE(p.full_len()  == 35);
    // da_init creates num_vars base vectors in the pool, so remain = pool_size - num_vars
    REQUIRE(p.remain()    == 97);  // 100 - 3 base vectors for 3 variables

    // Assign and free a slot
    unsigned idx = p.assign();
    REQUIRE(p.remain() == 96);
    p.free(idx);
    REQUIRE(p.remain() == 97);

    da::da_clear();
}

// ======================================================================
// da_make_env + da_select_env
// ======================================================================

TEST_CASE("da_make_env creates a separate env and selects it", "[env]") {
    // Start with one env
    da::da_init(4, 3, 100);
    da::DAEnv* default_p = &da::da_current_env();

    // Make a second env
    da::DAEnv& env2 = da::da_make_env(6, 2, 50);
    REQUIRE(env2.layout().num_vars()  == 2);
    REQUIRE(env2.layout().max_order() == 6);
    REQUIRE(env2.layout().full_len()  == 28);  // C(2+6,6) = 28

    // da_make_env selects env2 as current
    REQUIRE(&da::da_current_env() == &env2);
    REQUIRE(&da::da_current_env() != default_p);

    // Switch back to the default env
    da::da_select_env(*default_p);
    REQUIRE(&da::da_current_env() == default_p);

    // Clean up: delete env2 (da_make_env returned a raw-new pointer)
    delete &env2;

    da::da_clear();
}

TEST_CASE("da_select_env switches current env", "[env]") {
    da::da_init(4, 3, 100);
    da::DAEnv& env_a = da::da_current_env();

    // Create a second env on the stack
    da::DAEnv env_b(2, 2, 50, false);  // C(2+2,2) = 6

    da::da_select_env(env_b);
    REQUIRE(&da::da_current_env() == &env_b);
    REQUIRE(da::da_current_env().layout().full_len() == 6);

    // Switch back
    da::da_select_env(env_a);
    REQUIRE(&da::da_current_env() == &env_a);
    REQUIRE(da::da_current_env().layout().full_len() == 35);

    da::da_clear();
}

TEST_CASE("da_exchange_env swaps the current env and returns the previous one", "[env]") {
    da::da_init(4, 3, 100);
    da::DAEnv* env_a = &da::da_current_env();
    da::DAEnv env_b(2, 2, 50, false);

    REQUIRE(da::da_exchange_env(&env_b) == env_a);
    REQUIRE(&da::da_current_env() == &env_b);
    REQUIRE(da::da_exchange_env(env_a) == &env_b);
    REQUIRE(&da::da_current_env() == env_a);

    REQUIRE(da::da_exchange_env(nullptr) == env_a);
    REQUIRE_THROWS_AS(da::da_current_env(), std::runtime_error);
    REQUIRE(da::da_exchange_env(env_a) == nullptr);

    da::da_clear();
}

// ======================================================================
// Tear down and re-init
// ======================================================================

TEST_CASE("da_clear then da_init in same process", "[env]") {
    da::da_init(4, 3, 100);
    REQUIRE(da::da_current_env().layout().full_len() == 35);

    da::da_clear();
    REQUIRE_THROWS_AS(da::da_current_env(), std::runtime_error);

    // Re-init with different parameters
    da::da_init(2, 2, 50);
    REQUIRE(da::da_current_env().layout().full_len() == 6);  // C(2+2,2)

    da::da_clear();
}

TEST_CASE("multiple da_init calls replace the default env", "[env]") {
    da::da_init(4, 3, 100);
    da::DAEnv* p1 = &da::da_current_env();

    da::da_init(3, 2, 50);   // re-initialize
    da::DAEnv* p2 = &da::da_current_env();

    // Both calls succeed; the second created a fresh env
    REQUIRE(p2 != nullptr);
    REQUIRE(da::da_current_env().layout().max_order() == 3);
    REQUIRE(da::da_current_env().layout().num_vars()  == 2);
    // p1 was deleted by the second da_init (unique_ptr reset)
    (void)p1;

    da::da_clear();
}

// ======================================================================
// import<double>: same-env slot copy
// ======================================================================

TEST_CASE("DAEnv::import copies a slot within compatible envs", "[env]") {
    da::da_init(4, 3, 100);
    da::DAEnv& env = da::da_current_env();

    auto& p = env.pool<double>();

    // Write known data into slot src
    unsigned src = p.assign();
    p.slot(src)[0] = 1.0;
    p.slot(src)[1] = 2.5;
    p.set_len(src, 2);

    // Import from same env (legal: layouts are identical)
    unsigned dst = env.import<double>(env, src);
    REQUIRE(p.len(dst)     == 2);
    REQUIRE(p.slot(dst)[0] == 1.0);
    REQUIRE(p.slot(dst)[1] == 2.5);

    p.free(src);
    p.free(dst);
    da::da_clear();
}

TEST_CASE("DAEnv::import throws on mismatched layouts", "[env]") {
    da::DAEnv env_a(4, 3, 50, false);  // full_len = 35
    da::DAEnv env_b(2, 2, 50, false);  // full_len = 6

    auto& pa = env_a.pool<double>();
    unsigned src = pa.assign();
    pa.slot(src)[0] = 3.14;
    pa.set_len(src, 1);

    REQUIRE_THROWS_AS(env_b.import<double>(env_a, src),
                      std::invalid_argument);
}

// ======================================================================
// same_env and check_env helpers
// ======================================================================

TEST_CASE("same_env returns true for identical pointers", "[env]") {
    da::da_init(4, 3, 100);
    da::DAEnv* p = &da::da_current_env();
    REQUIRE(da::same_env(p, p));
    da::da_clear();
}

TEST_CASE("same_env returns false for different envs", "[env]") {
    da::DAEnv env_a(4, 3, 50, false);
    da::DAEnv env_b(4, 3, 50, false);
    REQUIRE_FALSE(da::same_env(&env_a, &env_b));
}

#if DA_CHECK_ENV
TEST_CASE("check_env throws for different envs when DA_CHECK_ENV=1", "[env]") {
    da::DAEnv env_a(4, 3, 50, false);
    da::DAEnv env_b(4, 3, 50, false);
    REQUIRE_THROWS_AS(da::check_env(&env_a, &env_b), std::logic_error);
}

TEST_CASE("check_env does not throw for same env when DA_CHECK_ENV=1", "[env]") {
    da::DAEnv env_a(4, 3, 50, false);
    REQUIRE_NOTHROW(da::check_env(&env_a, &env_a));
}
#endif
