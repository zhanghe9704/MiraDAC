/**
 * @file test_multienv.cc
 * @brief Tests for the multi-environment policy (Stage 8).
 *
 * @details Proves the multi-env design from plan A.4 / A.5:
 *   1. Two DAEnvs with different (order,nv) coexist; NDAs carry the right env.
 *   2. import/convert moves a vector between two SAME-layout envs (copy+retag)
 *      and the result computes correctly in the target env.
 *   3. Binary ops between vectors from DIFFERENT envs THROW when DA_CHECK_ENV=1.
 *   4. A different-layout import throws.
 *
 * NUMERICAL ONLY — no SymEngine — runs in both OFF and ON builds.
 * Guard-specific assertions (item 3) are gated with #if DA_CHECK_ENV so they
 * compile cleanly out when DA_CHECK_ENV=0.
 */

#include "catch.hpp"
#include "da/da.h"

#include <cmath>
#include <stdexcept>
#include <vector>

// ============================================================================
// 1. Two DAEnvs with different (order,nv) coexist
// ============================================================================

TEST_CASE("Two DAEnvs with different layouts coexist", "[multienv]") {
    // env_a: order=4, nv=3  -> full_len = C(7,4) = 35
    // env_b: order=6, nv=2  -> full_len = C(8,6) = 28
    da::DAEnv env_a(4, 3, 100, false);
    da::DAEnv env_b(6, 2, 50,  false);

    REQUIRE(env_a.layout().num_vars()  == 3);
    REQUIRE(env_a.layout().max_order() == 4);
    REQUIRE(env_a.layout().full_len()  == 35);

    REQUIRE(env_b.layout().num_vars()  == 2);
    REQUIRE(env_b.layout().max_order() == 6);
    REQUIRE(env_b.layout().full_len()  == 28);

    REQUIRE(&env_a != &env_b);
    REQUIRE(env_a.pool<double>().poolsize() == 100);
    REQUIRE(env_b.pool<double>().poolsize() == 50);
}

TEST_CASE("NDA constructed in env_a carries env_a pointer", "[multienv]") {
    da::DAEnv env_a(4, 3, 200, false);
    da::DAEnv env_b(2, 2, 100, false);

    da::da_select_env(env_a);
    da::NDA v_a;     // alloc'd from current env = env_a
    REQUIRE(v_a.env_ == &env_a);

    da::da_select_env(env_b);
    da::NDA v_b;     // alloc'd from current env = env_b
    REQUIRE(v_b.env_ == &env_b);

    REQUIRE(v_a.env_ != v_b.env_);
}

// ============================================================================
// 2. import: SAME-layout cross-env copy + retag; correct computation
// ============================================================================

TEST_CASE("import between same-layout envs copies data correctly", "[multienv]") {
    da::DAEnv env_a(4, 3, 200, false);
    da::DAEnv env_b(4, 3, 200, false);

    REQUIRE(env_a.layout().full_len() == env_b.layout().full_len());

    da::Pool<double>& pa = env_a.pool<double>();
    unsigned src_slot = pa.assign();
    pa.slot(src_slot)[0] = 3.14;
    pa.slot(src_slot)[1] = 2.72;
    pa.set_len(src_slot, 2);

    // Same layout -> import must succeed
    unsigned dst_slot = env_b.import<double>(env_a, src_slot);

    da::Pool<double>& pb = env_b.pool<double>();
    REQUIRE(pb.len(dst_slot)     == 2);
    REQUIRE(pb.slot(dst_slot)[0] == Approx(3.14));
    REQUIRE(pb.slot(dst_slot)[1] == Approx(2.72));

    // The imported slot lives in env_b's pool; env_b can operate on it
    da::da_select_env(env_b);
    da::NDA wrapper;
    da::detail::ad_copy(env_b.layout(), pb, dst_slot, wrapper.slot_);
    REQUIRE(wrapper.env_  == &env_b);
    REQUIRE(wrapper.con() == Approx(3.14));

    pa.free(src_slot);
    pb.free(dst_slot);
}

TEST_CASE("NDA imported into env_b computes correctly in env_b", "[multienv]") {
    // order=2, nv=2 for both: layout [1, x, y, x^2, xy, y^2]
    da::DAEnv env_a(2, 2, 100, false);
    da::DAEnv env_b(2, 2, 100, false);

    da::Pool<double>& pa = env_a.pool<double>();
    da::Pool<double>& pb = env_b.pool<double>();

    // Build v = 1 + x in env_a (index 0 = const, index 1 = x)
    unsigned src_slot = pa.assign();
    pa.slot(src_slot)[0] = 1.0;
    pa.slot(src_slot)[1] = 1.0;
    pa.set_len(src_slot, 2);

    // Import into env_b
    unsigned dst_slot = env_b.import<double>(env_a, src_slot);

    // Wrap in an NDA living in env_b
    da::da_select_env(env_b);
    da::NDA v_b;
    da::detail::ad_copy(env_b.layout(), pb, dst_slot, v_b.slot_);
    REQUIRE(v_b.env_ == &env_b);
    REQUIRE(v_b.con() == Approx(1.0));

    // Build constant 2.0 in env_b and add -> should yield 3 + x
    da::NDA two;
    da::detail::ad_const(env_b.layout(), pb, two.slot_, 2.0);
    da::NDA result = v_b + two;
    REQUIRE(result.env_ == &env_b);
    REQUIRE(result.con() == Approx(3.0));
    // Linear term (index 1 = x) preserved
    REQUIRE(pb.slot(result.slot_)[1] == Approx(1.0));

    pa.free(src_slot);
    pb.free(dst_slot);
}

// ============================================================================
// 3. Cross-env binary operations throw when DA_CHECK_ENV=1
// ============================================================================

#if DA_CHECK_ENV

TEST_CASE("NDA + NDA from different envs throws (DA_CHECK_ENV=1)", "[multienv][guard]") {
    da::DAEnv env_a(4, 3, 100, false);
    da::DAEnv env_b(4, 3, 100, false);

    da::da_select_env(env_a);
    da::NDA v_a;
    da::da_select_env(env_b);
    da::NDA v_b;

    REQUIRE_THROWS_AS(v_a + v_b, std::logic_error);
    REQUIRE_THROWS_AS(v_a - v_b, std::logic_error);
}

TEST_CASE("NDA * NDA from different envs throws (DA_CHECK_ENV=1)", "[multienv][guard]") {
    da::DAEnv env_a(3, 2, 100, false);
    da::DAEnv env_b(3, 2, 100, false);

    da::da_select_env(env_a);
    da::NDA v_a;
    da::da_select_env(env_b);
    da::NDA v_b;

    REQUIRE_THROWS_AS(v_a * v_b, std::logic_error);
    REQUIRE_THROWS_AS(v_a / v_b, std::logic_error);
}

TEST_CASE("compound op between different-env NDAs throws (DA_CHECK_ENV=1)", "[multienv][guard]") {
    da::DAEnv env_a(2, 2, 100, false);
    da::DAEnv env_b(2, 2, 100, false);

    da::da_select_env(env_a);
    da::NDA v_a;
    da::da_select_env(env_b);
    da::NDA v_b;

    // operator+= calls check_env internally
    REQUIRE_THROWS_AS(v_a += v_b, std::logic_error);
}

TEST_CASE("check_env free function throws for different envs (DA_CHECK_ENV=1)", "[multienv][guard]") {
    da::DAEnv env_a(2, 2, 50, false);
    da::DAEnv env_b(2, 2, 50, false);

    REQUIRE_THROWS_AS(da::check_env(&env_a, &env_b), std::logic_error);
    REQUIRE_NOTHROW(da::check_env(&env_a, &env_a));
}

TEST_CASE("same_env returns false for vectors in different envs", "[multienv][guard]") {
    da::DAEnv env_a(2, 2, 50, false);
    da::DAEnv env_b(2, 2, 50, false);

    da::da_select_env(env_a);
    da::NDA v_a;
    da::da_select_env(env_b);
    da::NDA v_b;

    REQUIRE_FALSE(da::same_env(v_a, v_b));
}

#endif // DA_CHECK_ENV

// ============================================================================
// 4. Different-layout import throws
// ============================================================================

TEST_CASE("import between envs with different full_len throws", "[multienv]") {
    da::DAEnv env_a(4, 3, 50, false);   // full_len = 35
    da::DAEnv env_b(2, 2, 50, false);   // full_len = 6

    da::Pool<double>& pa = env_a.pool<double>();
    unsigned src = pa.assign();
    pa.slot(src)[0] = 1.0;
    pa.set_len(src, 1);

    REQUIRE_THROWS_AS(env_b.import<double>(env_a, src), std::invalid_argument);

    pa.free(src);
}

TEST_CASE("import from env with different order throws", "[multienv]") {
    da::DAEnv env_a(6, 3, 50, false);   // full_len = C(9,6) = 84
    da::DAEnv env_b(4, 3, 50, false);   // full_len = C(7,4) = 35

    da::Pool<double>& pa = env_a.pool<double>();
    unsigned src = pa.assign();
    pa.slot(src)[0] = 2.0;
    pa.set_len(src, 1);

    REQUIRE_THROWS_AS(env_b.import<double>(env_a, src), std::invalid_argument);

    pa.free(src);
}

TEST_CASE("import from env with different num_vars throws", "[multienv]") {
    da::DAEnv env_a(4, 6, 50, false);   // 6 variables
    da::DAEnv env_b(4, 3, 50, false);   // 3 variables

    da::Pool<double>& pa = env_a.pool<double>();
    unsigned src = pa.assign();
    pa.slot(src)[0] = 1.5;
    pa.set_len(src, 1);

    REQUIRE_THROWS_AS(env_b.import<double>(env_a, src), std::invalid_argument);

    pa.free(src);
}

// ============================================================================
// 5. Independent computation in each env is correct
// ============================================================================

TEST_CASE("Independent computation in two envs gives correct results", "[multienv]") {
    da::DAEnv env_a(2, 2, 200, false);
    da::DAEnv env_b(3, 3, 200, false);

    // ---- env_a: compute x^2 in order-2, nv=2 env ----
    da::da_select_env(env_a);
    da::Pool<double>& pa = env_a.pool<double>();

    da::NDA x_a;
    // x: constant=0, coeff of x (index 1) = 1
    da::detail::ad_const(env_a.layout(), pa, x_a.slot_, 0.0);
    {
        std::vector<int> idx{1, 0};
        da::detail::ad_pok(env_a.layout(), pa, x_a.slot_, idx.data(), 2, 1.0);
    }
    da::NDA x2_a = x_a * x_a;
    REQUIRE(x2_a.env_ == &env_a);
    {
        std::vector<int> idx{2, 0};
        REQUIRE(x2_a.element(idx) == Approx(1.0));
    }

    // ---- env_b: compute y^2 in order-3, nv=3 env ----
    da::da_select_env(env_b);
    da::Pool<double>& pb = env_b.pool<double>();

    da::NDA y_b;
    da::detail::ad_const(env_b.layout(), pb, y_b.slot_, 0.0);
    {
        std::vector<int> idx{0, 1, 0};
        da::detail::ad_pok(env_b.layout(), pb, y_b.slot_, idx.data(), 3, 1.0);
    }
    da::NDA y2_b = y_b * y_b;
    REQUIRE(y2_b.env_ == &env_b);
    {
        std::vector<int> idx{0, 2, 0};
        REQUIRE(y2_b.element(idx) == Approx(1.0));
    }

    // Envs are completely independent
    REQUIRE(x2_a.env_ != y2_b.env_);
}

// ============================================================================
// 6. same_env positive case
// ============================================================================

TEST_CASE("same_env returns true for two vectors in the same env", "[multienv]") {
    da::DAEnv env_a(2, 2, 100, false);
    da::da_select_env(env_a);

    da::NDA v1, v2;

    REQUIRE(da::same_env(v1, v2));
    REQUIRE(v1.env_ == v2.env_);
    REQUIRE(v1.env_ == &env_a);
}

// ============================================================================
// 7. da_make_env lifecycle
// ============================================================================

TEST_CASE("da_make_env lifecycle: two heap envs, switch, delete", "[multienv]") {
    da::da_init(4, 3, 200);
    da::DAEnv* default_env = &da::da_current_env();

    // Make a second heap env (da_make_env owns nothing; caller must delete)
    da::DAEnv& env2 = da::da_make_env(2, 2, 100);
    REQUIRE(&da::da_current_env() == &env2);
    REQUIRE(env2.layout().num_vars()  == 2);
    REQUIRE(env2.layout().max_order() == 2);

    // A DAVector holds a raw pointer to its env, so the env MUST outlive the
    // vector. Scope every vector so it is destroyed before its env goes away.
    {
        da::NDA v2;                       // lives in env2
        REQUIRE(v2.env_ == &env2);
    }                                     // v2 destroyed here (env2 still alive)

    // Switch back to default env and use it in its own scope.
    da::da_select_env(*default_env);
    {
        da::NDA v_def;
        REQUIRE(v_def.env_ == default_env);
        REQUIRE(v_def.env_ != &env2);
    }                                     // v_def destroyed here (default env alive)

    // No vectors reference env2 now -> safe to delete.
    delete &env2;

    // No vectors reference the default env now -> safe to clear.
    da::da_clear();
}

// ===========================================================================
// 8. Environment retirement
//
// da_clear() used to delete the default DAEnv outright, while ~DAVector
// reads poolsize() through env_. Any vector still in scope at that point
// was a use-after-free. da_clear() now retires an env that still has live
// vectors: its memory is released, the shell stays valid, and the late
// destructors see poolsize() == 0. Run these under ASan.
// ===========================================================================

TEST_CASE("da_clear with live vectors is safe", "[multienv][retire]") {
    da::da_init(4, 2, 500);

    da::NDA x = da::base[0];
    da::NDA y = 1.0 + da::base[1];
    da::NDA z = x * y;
    REQUIRE(z.n_element() > 0);

    da::da_clear();                 // x, y, z are still in scope

    // Their destructors run at the end of this TEST_CASE, against the
    // retired env. Nothing here may touch the released pool.
}

TEST_CASE("re-init with live vectors is safe", "[multienv][retire]") {
    da::da_init(4, 2, 500);
    da::NDA old_vec = 1.0 + da::base[0];
    da::DAEnv* first = &da::da_current_env();

    da::da_init(3, 2, 200);         // replaces the default env
    REQUIRE(&da::da_current_env() != first);

    {
        da::NDA fresh = 2.0 + da::base[0];
        REQUIRE(fresh.env_ == &da::da_current_env());
        REQUIRE(old_vec.env_ == first);   // still points at the retired env
    }

    da::da_clear();                 // old_vec is still alive here too
}

TEST_CASE("live_slots reports outstanding vectors", "[multienv][retire]") {
    da::da_init(4, 2, 500);
    da::DAEnv& env = da::da_current_env();

    // da_init created the base vectors, so start from that baseline.
    unsigned base_slots = env.live_slots();
    {
        da::NDA a, b, c;
        REQUIRE(env.live_slots() == base_slots + 3);
    }
    REQUIRE(env.live_slots() == base_slots);
    REQUIRE_FALSE(env.retired());

    da::da_clear();
}

TEST_CASE("da_destroy_env handles both cases", "[multienv][retire]") {
    da::da_init(4, 2, 500);
    da::DAEnv* main_env = &da::da_current_env();

    SECTION("no live vectors: env is deleted") {
        da::DAEnv& env2 = da::da_make_env(3, 2, 100);
        { da::NDA v; REQUIRE(v.env_ == &env2); }
        da::da_select_env(*main_env);
        da::da_destroy_env(env2);          // nothing references it
    }

    SECTION("live vectors: env is retired, not deleted") {
        da::DAEnv& env2 = da::da_make_env(3, 2, 100);
        da::NDA v;                          // deliberately outlives the call
        REQUIRE(v.env_ == &env2);
        da::da_select_env(*main_env);
        da::da_destroy_env(env2);
        REQUIRE(env2.retired());            // shell still readable
        REQUIRE(env2.live_slots() == 0);    // pool released
    }                                       // v destroyed here, safely

    da::da_select_env(*main_env);
    da::da_clear();
}

// ============================================================================
// Library functions must not depend on the default env's global da::base
// ============================================================================

static void check_inv_map_identity() {
    da::NDA x = da::da_base(0);
    da::NDA y = da::da_base(1);
    std::vector<da::NDA> map = { 2.0 * x + 0.3 * y + 0.1 * x * x + 0.05 * x * y,
                                -0.4 * x + 1.5 * y + 0.2 * y * y };
    std::vector<da::NDA> inv(2), composed(2);
    REQUIRE_NOTHROW(da::inv_map(map, 2, inv));
    da::da_composition(map, inv, composed);
    REQUIRE((composed[0] - x).norm() < 1e-12);
    REQUIRE((composed[1] - y).norm() < 1e-12);
}

TEST_CASE("inv_map works in a non-default env", "[multienv][inv_map]") {
    da::da_init(4, 3, 500);                   // default env: different layout
    da::DAEnv* main_env = &da::da_current_env();
    da::DAEnv& env2 = da::da_make_env(5, 2, 1000);   // selected as current
    check_inv_map_identity();
    da::da_select_env(*main_env);
    da::da_destroy_env(env2);
    da::da_clear();
}

TEST_CASE("inv_map works when no default env exists", "[multienv][inv_map]") {
    da::da_clear();                           // global da::base is now empty
    da::DAEnv& env2 = da::da_make_env(5, 2, 1000);
    check_inv_map_identity();
    da::da_destroy_env(env2);
}
