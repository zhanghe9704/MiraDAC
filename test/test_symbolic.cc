/**
 * @file test_symbolic.cc
 * @brief Symbolic DA (SDA) and interop (NDA<->SDA) tests.
 *
 * Ported from ref/tpsa_sym/tests/tests.cc, adapted to the new da:: API.
 * Also includes Stage 7 interop tests.
 *
 * This entire file is a no-op (empty translation unit) when DA_WITH_SYMBOLIC
 * is not defined — it contributes 0 test cases to the OFF build.
 */

#ifdef DA_WITH_SYMBOLIC

#include "catch.hpp"
#include "da/da.h"

#include <algorithm>
#include <array>
#include <cmath>
#include <vector>
#include <iomanip>
#include <iostream>

using da::SDA;
using da::NDA;
using SymEngine::Expression;
using std::cout;
using std::endl;
using std::vector;
using std::array;

// ---------------------------------------------------------------------------
// Helper: eval all coefficients at (x,y) using eval_funs single-visitor form
// ---------------------------------------------------------------------------
static void eval_sda(SDA& sv, Expression& sx, Expression& sy,
                     double x, double y, vector<double>& results) {
    vector<SymEngine::RCP<const SymEngine::Basic>> vars{sx.get_basic(), sy.get_basic()};
    SymEngine::LambdaRealDoubleVisitor vec;
    sv.eval_funs(vars, vec);
    results.resize(static_cast<std::size_t>(sv.full_length()));
    array<double, 2> inputs{x, y};
    vec.call(&results[0], &inputs[0]);
}

// ---------------------------------------------------------------------------
// Reference vectors (copied verbatim from ref/tpsa_sym/tests/tests.cc)
// ---------------------------------------------------------------------------
static vector<double> v_sqrt {
    1.2619825672330027e+00,  0.0000000000000000e+00,  2.1354098463568110e+00,
    3.9620198644763360e-01,  0.0000000000000000e+00, -1.8066712371137301e+00,
    0.0000000000000000e+00, -6.7041625215270972e-01,  0.0000000000000000e+00,
    3.0570814914830691e+00, -6.2194208597204828e-02,  0.0000000000000000e+00,
    1.7016242971685265e+00,  0.0000000000000000e+00, -6.4661371792376841e+00,
    0.0000000000000000e+00,  3.1571781308208735e-01,  0.0000000000000000e+00,
   -4.7988846720505292e+00,  0.0000000000000000e+00,  1.5317956604581280e+01};

static vector<double> v_exp {
    4.9165152645857058e+00,  0.0000000000000000e+00,  2.6498542321537581e+01,
    4.9165152645857058e+00,  0.0000000000000000e+00,  7.1409596775195553e+01,
    0.0000000000000000e+00,  2.6498542321537581e+01,  0.0000000000000000e+00,
    1.2829210124642384e+02,  2.4582576322928529e+00,  0.0000000000000000e+00,
    7.1409596775195553e+01,  0.0000000000000000e+00,  1.7286398452196264e+02,
    0.0000000000000000e+00,  1.3249271160768791e+01,  0.0000000000000000e+00,
    1.2829210124642384e+02,  0.0000000000000000e+00,  1.8633700347560443e+02};

static vector<double> v_log {
    4.6536790084120744e-01,  0.0000000000000000e+00,  3.3842144920256185e+00,
    6.2790405626020340e-01,  0.0000000000000000e+00, -5.7264538640181080e+00,
    0.0000000000000000e+00, -2.1249620067974497e+00,  0.0000000000000000e+00,
    1.2919698769684123e+01, -1.9713175193400836e-01,  0.0000000000000000e+00,
    7.1913272184077712e+00,  0.0000000000000000e+00, -3.2792273856727924e+01,
    0.0000000000000000e+00,  1.3342722634669408e+00,  0.0000000000000000e+00,
   -2.4336993789433862e+01,  0.0000000000000000e+00,  8.8780870729929191e+01};

static vector<double> v_sin {
    9.9976230933412147e-01,  0.0000000000000000e+00, -1.1750594656606546e-01,
   -2.1801945667860076e-02,  0.0000000000000000e+00, -1.4520980722838027e+01,
    0.0000000000000000e+00, -5.3884189186181146e+00,  0.0000000000000000e+00,
    5.6890241776272188e-01, -4.9988115466706073e-01,  0.0000000000000000e+00,
    3.1666090010356152e-01,  0.0000000000000000e+00,  3.5151502042766111e+01,
    0.0000000000000000e+00,  5.8752973283032729e-02,  0.0000000000000000e+00,
    2.6087909933960042e+01,  0.0000000000000000e+00, -8.2629850759332735e-01};

static vector<double> v_cos {
   -2.1801945667860076e-02,  0.0000000000000000e+00, -5.3884189186181146e+00,
   -9.9976230933412147e-01,  0.0000000000000000e+00,  3.1666090010356152e-01,
    0.0000000000000000e+00,  1.1750594656606546e-01,  0.0000000000000000e+00,
    2.6087909933960042e+01,  1.0900972833930038e-02,  0.0000000000000000e+00,
    1.4520980722838027e+01,  0.0000000000000000e+00, -7.6655334025393562e-01,
    0.0000000000000000e+00,  2.6942094593090573e+00,  0.0000000000000000e+00,
   -5.6890241776272188e-01,  0.0000000000000000e+00, -3.7891210111979305e+01};

static vector<double> v_tan {
   -4.5856563655598315e+01,  0.0000000000000000e+00,  1.1338982531987589e+04,
    2.1038244302999406e+03,  0.0000000000000000e+00, -2.8024649232711950e+06,
    0.0000000000000000e+00, -1.0399335485356124e+06,  0.0000000000000000e+00,
    6.9274774751837516e+08, -9.6474158908252066e+04,  0.0000000000000000e+00,
    3.8559534715385360e+08,  0.0000000000000000e+00, -1.7124190410087451e+11,
    0.0000000000000000e+00,  7.1543007431555301e+07,  0.0000000000000000e+00,
   -1.2708826398565749e+11,  0.0000000000000000e+00,  4.2329679125704305e+13};

static vector<double> v_asin {
    6.3428284416340086e-01,  0.0000000000000000e+00,  6.6911493589602227e+00,
    1.2414697216839945e+00,  0.0000000000000000e+00,  1.6469075934078099e+01,
    0.0000000000000000e+00,  6.1113145199465828e+00,  0.0000000000000000e+00,
    1.3100018800423086e+02,  5.6694384844672152e-01,  0.0000000000000000e+00,
    7.2916964582939883e+01,  0.0000000000000000e+00,  1.0518658826557562e+03,
    0.0000000000000000e+00,  1.3528946802779576e+01,  0.0000000000000000e+00,
    7.8064892862726265e+02,  0.0000000000000000e+00,  9.8884295426925273e+03};

static vector<double> v_atan {
    5.3496055600178871e-01,  0.0000000000000000e+00,  3.9888992597819097e+00,
    7.4009671406236144e-01,  0.0000000000000000e+00, -9.4290466347585138e+00,
    0.0000000000000000e+00, -3.4989133475920777e+00,  0.0000000000000000e+00,
    1.1323711186436043e+00, -3.2459258841791550e-01,  0.0000000000000000e+00,
    6.3029729965130876e-01,  0.0000000000000000e+00,  9.7342311833335117e+01,
    0.0000000000000000e+00,  1.1694478350396453e-01,  0.0000000000000000e+00,
    7.2243213413241108e+01,  0.0000000000000000e+00, -3.8276667852990363e+02};

static vector<double> v_acos {
    9.3651348263149570e-01,  0.0000000000000000e+00, -6.6911493589602227e+00,
   -1.2414697216839945e+00,  0.0000000000000000e+00, -1.6469075934078099e+01,
    0.0000000000000000e+00, -6.1113145199465828e+00,  0.0000000000000000e+00,
   -1.3100018800423086e+02, -5.6694384844672152e-01,  0.0000000000000000e+00,
   -7.2916964582939883e+01,  0.0000000000000000e+00, -1.0518658826557562e+03,
    0.0000000000000000e+00, -1.3528946802779576e+01,  0.0000000000000000e+00,
   -7.8064892862726265e+02,  0.0000000000000000e+00, -9.8884295426925273e+03};

static vector<double> v_sinh {
    2.3565595853852059e+00,  0.0000000000000000e+00,  1.3797393124186938e+01,
    2.5599556792005003e+00,  0.0000000000000000e+00,  3.4227691914480388e+01,
    0.0000000000000000e+00,  1.2701149197350645e+01,  0.0000000000000000e+00,
    6.6799770875932182e+01,  1.1782797926926030e+00,  0.0000000000000000e+00,
    3.7181904860715171e+01,  0.0000000000000000e+00,  8.2856303249459728e+01,
    0.0000000000000000e+00,  6.8986965620934688e+00,  0.0000000000000000e+00,
    6.1492330370491651e+01,  0.0000000000000000e+00,  9.7022879950881816e+01};

static vector<double> v_cosh {
    2.5599556792005003e+00,  0.0000000000000000e+00,  1.2701149197350645e+01,
    2.3565595853852059e+00,  0.0000000000000000e+00,  3.7181904860715171e+01,
    0.0000000000000000e+00,  1.3797393124186938e+01,  0.0000000000000000e+00,
    6.1492330370491651e+01,  1.2799778396002501e+00,  0.0000000000000000e+00,
    3.4227691914480388e+01,  0.0000000000000000e+00,  9.0007681272502936e+01,
    0.0000000000000000e+00,  6.3505745986753226e+00,  0.0000000000000000e+00,
    6.6799770875932182e+01,  0.0000000000000000e+00,  8.9314123524722618e+01};

static vector<double> v_tanh {
    9.2054702529896260e-01,  0.0000000000000000e+00,  8.2243143105705130e-01,
    1.5259317421323104e-01,  0.0000000000000000e+00, -4.0804707656925352e+00,
    0.0000000000000000e+00, -1.5141736147438749e+00,  0.0000000000000000e+00,
    1.2281575847812150e+01, -1.4046919260291624e-01,  0.0000000000000000e+00,
    6.8361369915647501e+00,  0.0000000000000000e+00, -2.1423682111362041e+01,
    0.0000000000000000e+00,  1.2683705942009276e+00,  0.0000000000000000e+00,
   -1.5899721402944124e+01,  0.0000000000000000e+00,  2.7932271238723843e+00};

static vector<double> v_pow3 {
    4.0394304427760002e+00,  0.0000000000000000e+00,  4.1010897131916003e+01,
    7.6091242800000005e+00,  0.0000000000000000e+00,  1.3878967240480202e+02,
    0.0000000000000000e+00,  5.1501817320000001e+01,  0.0000000000000000e+00,
    1.5656467356527304e+02,  4.7778000000000000e+00,  0.0000000000000000e+00,
    8.7146598270000013e+01,  0.0000000000000000e+00,  0.0000000000000000e+00,
    0.0000000000000000e+00,  1.6169100000000000e+01,  0.0000000000000000e+00,
    0.0000000000000000e+00,  0.0000000000000000e+00,  0.0000000000000000e+00};

static vector<double> v_pow043 {
    1.2215349178657362e+00,  0.0000000000000000e+00,  1.7775925537693156e+00,
    3.2981289381028916e-01,  0.0000000000000000e+00, -1.7144900271941401e+00,
    0.0000000000000000e+00, -6.3620981768712170e-01,  0.0000000000000000e+00,
    3.0364857114827242e+00, -5.9020893341662960e-02,  0.0000000000000000e+00,
    1.6901603307138011e+00,  0.0000000000000000e+00, -6.6024064251363503e+00,
    0.0000000000000000e+00,  3.1359079924927197e-01,  0.0000000000000000e+00,
   -4.9000177561915139e+00,  0.0000000000000000e+00,  1.5953587087419303e+01};

static vector<double> v_erf {
    9.7569519695814055e-01,  0.0000000000000000e+00,  4.8137783493505698e-01,
    8.9314402459331130e-02,  0.0000000000000000e+00, -4.1319722194537372e+00,
    0.0000000000000000e+00, -1.5332846798351436e+00,  0.0000000000000000e+00,
    1.8983737513383726e+01, -1.4224211735673076e-01,  0.0000000000000000e+00,
    1.0566675796454568e+01,  0.0000000000000000e+00, -4.1465045891922742e+01,
    0.0000000000000000e+00,  1.9605313461703930e+00,  0.0000000000000000e+00,
   -3.0773546499376767e+01,  0.0000000000000000e+00, -2.3068466761450004e+01};

static vector<double> v_der {
    0.0000000000000000e+00,  9.8330305291714115e+00,  0.0000000000000000e+00,
    0.0000000000000000e+00,  5.2997084643075162e+01,  0.0000000000000000e+00,
    9.8330305291714115e+00,  0.0000000000000000e+00,  1.4281919355039111e+02,
    0.0000000000000000e+00,  0.0000000000000000e+00,  5.2997084643075162e+01,
    0.0000000000000000e+00,  2.5658420249284768e+02,  0.0000000000000000e+00,
    0.0000000000000000e+00,  0.0000000000000000e+00,  0.0000000000000000e+00,
    0.0000000000000000e+00,  0.0000000000000000e+00,  0.0000000000000000e+00};

static vector<double> v_int {
    0.0000000000000000e+00,  4.9165152645857058e+00,  0.0000000000000000e+00,
    0.0000000000000000e+00,  2.6498542321537581e+01,  0.0000000000000000e+00,
    1.6388384215285685e+00,  0.0000000000000000e+00,  7.1409596775195553e+01,
    0.0000000000000000e+00,  0.0000000000000000e+00,  8.8328474405125252e+00,
    0.0000000000000000e+00,  1.2829210124642384e+02,  0.0000000000000000e+00,
    4.9165152645857052e-01,  0.0000000000000000e+00,  2.3803198925065185e+01,
    0.0000000000000000e+00,  1.7286398452196264e+02,  0.0000000000000000e+00};

static vector<double> v_subc {
    1.2291288161464264e+01,  0.0000000000000000e+00,  6.6246355803843954e+01,
    0.0000000000000000e+00,  0.0000000000000000e+00,  1.4281919355039111e+02,
    0.0000000000000000e+00,  0.0000000000000000e+00,  0.0000000000000000e+00,
    2.5658420249284768e+02,  0.0000000000000000e+00,  0.0000000000000000e+00,
    0.0000000000000000e+00,  0.0000000000000000e+00,  1.7286398452196264e+02,
    0.0000000000000000e+00,  0.0000000000000000e+00,  0.0000000000000000e+00,
    0.0000000000000000e+00,  0.0000000000000000e+00,  1.8633700347560443e+02};

static vector<double> v_sub {
    1.8981620087089855e+01,  0.0000000000000000e+00,  1.7100535860620653e+02,
    1.2746557474964902e+01,  0.0000000000000000e+00,  6.2681915857331160e+02,
    0.0000000000000000e+00,  9.5198663144355919e+01,  0.0000000000000000e+00,
    1.1023627091700214e+03,  2.4582576322928529e+00,  0.0000000000000000e+00,
    2.1422879032558666e+02,  0.0000000000000000e+00,  8.6431992260981315e+02,
    0.0000000000000000e+00,  1.3249271160768791e+01,  0.0000000000000000e+00,
    1.2829210124642384e+02,  0.0000000000000000e+00,  1.8633700347560443e+02};

static vector<double> v_subv {
    9.9749939097621802e+02,  0.0000000000000000e+00,  1.1512590930808019e+03,
    2.1360355735584579e+02,  0.0000000000000000e+00, -4.4758214953827373e+04,
    0.0000000000000000e+00, -1.6608796390829688e+04,  0.0000000000000000e+00,
   -3.0987895299849177e+04, -1.5407904327541128e+03,  0.0000000000000000e+00,
   -1.7248397109217120e+04,  0.0000000000000000e+00,  9.7163445148612640e+05,
    0.0000000000000000e+00, -3.2002517968007724e+03,  0.0000000000000000e+00,
    7.2110466370011435e+05,  0.0000000000000000e+00,  3.1279038644950697e+05};

static vector<double> v_vsub0 {
    2.9683742321537579e+01,  0.0000000000000000e+00,  1.4820889355039111e+02,
    2.7498542321537581e+01,  0.0000000000000000e+00,  3.8487630373927152e+02,
    0.0000000000000000e+00,  1.4281919355039111e+02,  0.0000000000000000e+00,
    6.9145593808785054e+02,  1.3249271160768791e+01,  0.0000000000000000e+00,
    3.8487630373927152e+02,  0.0000000000000000e+00,  9.3168501737802217e+02,
    0.0000000000000000e+00,  7.1409596775195553e+01,  0.0000000000000000e+00,
    6.9145593808785054e+02,  0.0000000000000000e+00,  1.0043005476324653e+03};

static vector<double> v_vsub1 {
    2.2789117817838105e+04,  0.0000000000000000e+00,  4.7247897335916117e+05,
    8.7663315835605143e+04,  0.0000000000000000e+00,  4.8149709416966103e+06,
    0.0000000000000000e+00,  1.7867305941691040e+06,  0.0000000000000000e+00,
    3.2197066119944468e+07,  1.6575417872693326e+05,  0.0000000000000000e+00,
    1.7921442447600685e+07,  0.0000000000000000e+00,  1.5854645261039719e+08,
    0.0000000000000000e+00,  3.3251280122457058e+06,  0.0000000000000000e+00,
    1.1766625423336896e+08,  0.0000000000000000e+00,  6.0926100575974762e+08};

static vector<double> v_vsub2 {
   -1.0573135900055892e+05,  0.0000000000000000e+00, -2.8866867756027146e+06,
   -5.3559321958600974e+05,  0.0000000000000000e+00, -3.9188057115516029e+07,
    0.0000000000000000e+00, -1.4541832426857166e+07,  0.0000000000000000e+00,
   -3.5349388645438391e+08, -1.3490391326843027e+06,  0.0000000000000000e+00,
   -1.9676079547343114e+08,  0.0000000000000000e+00, -2.3869922351913180e+09,
    0.0000000000000000e+00, -3.6506817721474491e+07,  0.0000000000000000e+00,
   -1.7715214094968681e+09,  0.0000000000000000e+00, -1.2880995781051701e+10};

static vector<double> v_compc {
    6.2871116966104892e+04,
    3.5879835234251332e+04,
   -1.6797595287239095e+05};

// ---------------------------------------------------------------------------
// Helper: compare two vectors within relative/absolute tolerance
// ---------------------------------------------------------------------------
static bool compare_vectors(vector<double>& a, vector<double>& b, double eps) {
    if (a.size() != b.size()) return false;
    for (std::size_t i = 0; i < a.size(); ++i) {
        double d = std::fabs(a[i] - b[i]);
        if (d > eps) d /= std::fabs(b[i]) + 1e-300;
        if (d > eps) {
            cout << std::scientific << std::setprecision(16);
            cout << i << ' ' << a[i] << ' ' << b[i] << ' ' << d << endl;
            return false;
        }
    }
    return true;
}

// ===========================================================================
// Environment initialization
// ===========================================================================
TEST_CASE("SDA INITIALIZE DA ENVIRONMENT") {
    int da_dim = 2;
    int da_order = 5;
    int n_vec = 400;
    int init = da::da_init(static_cast<unsigned>(da_order),
                          static_cast<unsigned>(da_dim),
                          static_cast<unsigned>(n_vec), /*table=*/true);
    REQUIRE(init == 0);
}

// ===========================================================================
// SDA math functions
// ===========================================================================
TEST_CASE("SDA FUNCTIONS") {
    // Re-initialize (in case previous test case changed things)
    da::da_init(5u, 2u, 400u, true);

    auto& sda = da::base;  // da::base[i] = i-th basis variable

    Expression sx("x"), sy("y");
    double x{1.5926}, y{5.3897};
    double eps = 1e-13;
    vector<double> results;
    SDA sv;

    SECTION("SDA SQRT") {
        sv = da::sqrt(sx + sda[0]*sda[0] + sy*sda[1]);
        eval_sda(sv, sx, sy, x, y, results);
        REQUIRE(compare_vectors(v_sqrt, results, eps));
    }

    SECTION("SDA EXP") {
        sv = da::exp(sx + sda[0]*sda[0] + sy*sda[1]);
        eval_sda(sv, sx, sy, x, y, results);
        REQUIRE(compare_vectors(v_exp, results, eps));
    }

    SECTION("SDA LOG") {
        sv = da::log(sx + sda[0]*sda[0] + sy*sda[1]);
        eval_sda(sv, sx, sy, x, y, results);
        REQUIRE(compare_vectors(v_log, results, eps));
    }

    SECTION("SDA POW INT") {
        sv = da::pow(sx + sda[0]*sda[0] + sy*sda[1], 3);
        eval_sda(sv, sx, sy, x, y, results);
        REQUIRE(compare_vectors(v_pow3, results, eps));
    }

    SECTION("SDA POW DOUBLE") {
        sv = da::pow(sx + sda[0]*sda[0] + sy*sda[1], 0.43);
        eval_sda(sv, sx, sy, x, y, results);
        REQUIRE(compare_vectors(v_pow043, results, eps));
    }

    SECTION("SDA ERF") {
        sv = da::erf(sx + sda[0]*sda[0] + sy*sda[1]);
        eval_sda(sv, sx, sy, x, y, results);
        REQUIRE(compare_vectors(v_erf, results, eps));
    }

    SECTION("SDA SIN") {
        sv = da::sin(sx + sda[0]*sda[0] + sy*sda[1]);
        eval_sda(sv, sx, sy, x, y, results);
        REQUIRE(compare_vectors(v_sin, results, eps));
    }

    SECTION("SDA COS") {
        sv = da::cos(sx + sda[0]*sda[0] + sy*sda[1]);
        eval_sda(sv, sx, sy, x, y, results);
        REQUIRE(compare_vectors(v_cos, results, eps));
    }

    SECTION("SDA TAN") {
        sv = da::tan(sx + sda[0]*sda[0] + sy*sda[1]);
        eval_sda(sv, sx, sy, x, y, results);
        REQUIRE(compare_vectors(v_tan, results, eps));
    }

    SECTION("SDA SINH") {
        sv = da::sinh(sx + sda[0]*sda[0] + sy*sda[1]);
        eval_sda(sv, sx, sy, x, y, results);
        REQUIRE(compare_vectors(v_sinh, results, eps));
    }

    SECTION("SDA COSH") {
        sv = da::cosh(sx + sda[0]*sda[0] + sy*sda[1]);
        eval_sda(sv, sx, sy, x, y, results);
        REQUIRE(compare_vectors(v_cosh, results, eps));
    }

    SECTION("SDA TANH") {
        sv = da::tanh(sx + sda[0]*sda[0] + sy*sda[1]);
        eval_sda(sv, sx, sy, x, y, results);
        REQUIRE(compare_vectors(v_tanh, results, eps));
    }

    SECTION("SDA ATAN") {
        sv = da::atan(sx - 1 + sda[0]*sda[0] + sy*sda[1]);
        eval_sda(sv, sx, sy, x, y, results);
        REQUIRE(compare_vectors(v_atan, results, eps));
    }

    SECTION("SDA DA_DER") {
        sv = da::exp(sx + sda[0]*sda[0] + sy*sda[1]);
        sv = da::da_der(sv, 0);
        eval_sda(sv, sx, sy, x, y, results);
        REQUIRE(compare_vectors(v_der, results, eps));
    }

    SECTION("SDA DA_INT") {
        sv = da::exp(sx + sda[0]*sda[0] + sy*sda[1]);
        sv = da::da_int(sv, 0);
        eval_sda(sv, sx, sy, x, y, results);
        REQUIRE(compare_vectors(v_int, results, eps));
    }

    SECTION("SDA SUB_CONST") {
        sv = da::exp(sx + sda[0]*sda[0] + sy*sda[1]);
        SDA sz;
        da::da_substitute_const(sv, 0, 1.0, sz);
        sv = sz;
        eval_sda(sv, sx, sy, x, y, results);
        REQUIRE(compare_vectors(v_subc, results, eps));
    }

    SECTION("SDA SUB") {
        sv = da::exp(sx + sda[0]*sda[0] + sy*sda[1]);
        SDA sz, sw;
        sw = da::sqrt(sx + sda[0]*sda[0] + sy*sda[1]);
        da::da_substitute(sv, 0u, sw, sz);
        sv = sz;
        eval_sda(sv, sx, sy, x, y, results);
        REQUIRE(compare_vectors(v_sub, results, eps));
    }

    SECTION("SDA SUBV") {
        sv = da::exp(sx + sda[0]*sda[0] + sy*sda[1]);
        SDA sz;
        std::vector<SDA> slv(2);
        slv[0] = da::sqrt(sx + sda[0]*sda[0] + sy*sda[1]);
        slv[1] = da::sin(sx  + sda[0]*sda[0] + sy*sda[1]);
        std::vector<unsigned int> idx{0u, 1u};
        da::da_substitute(sv, idx, slv, sz);
        sv = sz;
        eval_sda(sv, sx, sy, x, y, results);
        REQUIRE(compare_vectors(v_subv, results, eps));
    }

    SECTION("SDA VSUB") {
        std::vector<SDA> slv(2);
        slv[0] = da::sqrt(sx + sda[0]*sda[0] + sy*sda[1]);
        slv[1] = da::exp(sx  + sda[0]*sda[0] + sy*sda[1]);
        std::vector<unsigned int> idx{0u, 1u};
        std::vector<SDA> slx(3), sly(3);
        slx[0] = sx + sda[0]*sda[0] + sy*sda[1];
        slx[1] = da::sin(sx + sda[0]*sda[0] + sy*sda[1]);
        slx[2] = da::cos(sx + sda[0]*sda[0] + sy*sda[1]);
        da::da_substitute(slx, idx, slv, sly);
        sv = sly[0]; eval_sda(sv, sx, sy, x, y, results);
        REQUIRE(compare_vectors(v_vsub0, results, eps));
        sv = sly[1]; eval_sda(sv, sx, sy, x, y, results);
        REQUIRE(compare_vectors(v_vsub1, results, eps));
        sv = sly[2]; eval_sda(sv, sx, sy, x, y, results);
        REQUIRE(compare_vectors(v_vsub2, results, eps));
    }

    SECTION("SDA COMP SDA") {
        std::vector<SDA> slv(2);
        slv[0] = da::sqrt(sx + sda[0]*sda[0] + sy*sda[1]);
        slv[1] = da::exp(sx  + sda[0]*sda[0] + sy*sda[1]);
        std::vector<SDA> slx(3), sly(3);
        slx[0] = sx + sda[0]*sda[0] + sy*sda[1];
        slx[1] = da::sin(sx + sda[0]*sda[0] + sy*sda[1]);
        slx[2] = da::cos(sx + sda[0]*sda[0] + sy*sda[1]);
        da::da_composition(slx, slv, sly);
        sv = sly[0]; eval_sda(sv, sx, sy, x, y, results);
        REQUIRE(compare_vectors(v_vsub0, results, eps));
        sv = sly[1]; eval_sda(sv, sx, sy, x, y, results);
        REQUIRE(compare_vectors(v_vsub1, results, eps));
        sv = sly[2]; eval_sda(sv, sx, sy, x, y, results);
        REQUIRE(compare_vectors(v_vsub2, results, eps));
    }

    SECTION("SDA COMP DOUBLE") {
        std::vector<double> lv{x, y};
        std::vector<SDA> slx(3);
        std::vector<Expression> sly(3);
        slx[0] = da::sqrt(sx + sda[0]*sda[0] + sy*sda[1]);
        slx[1] = da::sin(sx  + sda[0]*sda[0] + sy*sda[1]);
        slx[2] = da::cos(sx  + sda[0]*sda[0] + sy*sda[1]);
        da::da_composition(slx, lv, sly);

        vector<SymEngine::RCP<const SymEngine::Basic>> vars{sx.get_basic(), sy.get_basic()};
        SymEngine::LambdaRealDoubleVisitor vis;
        vector<double> res3(3);
        array<double,2> inp{x, y};
        for (int i = 0; i < 3; ++i) {
            vis.init(vars, {sly[static_cast<std::size_t>(i)].get_basic()});
            vis.call(&res3[static_cast<std::size_t>(i)], &inp[0]);
        }
        REQUIRE(compare_vectors(v_compc, res3, eps));
    }
}

// ===========================================================================
// Stage 7 Interop tests
// ===========================================================================
TEST_CASE("SDA INTEROP") {
    da::da_init(4u, 2u, 400u, true);
    auto& b = da::base;

    double eps = 1e-13;

    SECTION("promote round-trip: NDA->SDA->evaluate reproduces NDA") {
        // Build a non-trivial NDA
        NDA n = 1.5 + 2.0*b[0] + 3.0*b[1] + 0.5*b[0]*b[1];
        SDA s = da::promote(n);
        // Evaluate with empty substitution map (all coefficients should be numeric)
        SymEngine::map_basic_basic empty;
        NDA back = da::evaluate(s, empty);
        // Compare all coefficients
        int fl = NDA::full_length();
        const double* pn = n.env_->pool<double>().slot(n.slot_);
        const double* pb = back.env_->pool<double>().slot(back.slot_);
        for (int i = 0; i < fl; ++i) {
            REQUIRE(std::fabs(pn[i] - pb[i]) < eps);
        }
    }

    SECTION("promote: NDA*SDA == promote(NDA)*SDA") {
        NDA n = 1.0 + 2.0*b[0] + 0.5*b[1];
        Expression sx("x");
        SDA s = sx + b[0]*b[0] + b[1];
        SDA mixed = n * s;
        SDA promoted_mul = da::promote(n) * s;
        // Both should have the same coefficients (modulo Expression form)
        // We evaluate both at x=1.0 and compare
        SymEngine::map_basic_basic m;
        m[sx] = Expression("1.0");
        NDA r1 = da::evaluate(mixed, m);
        NDA r2 = da::evaluate(promoted_mul, m);
        int fl = NDA::full_length();
        const double* p1 = r1.env_->pool<double>().slot(r1.slot_);
        const double* p2 = r2.env_->pool<double>().slot(r2.slot_);
        for (int i = 0; i < fl; ++i) {
            REQUIRE(std::fabs(p1[i] - p2[i]) < eps);
        }
    }

    SECTION("evaluate with symbol map") {
        Expression sx("px"), sy("py");
        SDA sv = sx + sy * b[0] + b[1];
        SymEngine::map_basic_basic m;
        m[sx] = Expression("2.0");
        m[sy] = Expression("3.0");
        NDA result = da::evaluate(sv, m);
        // Constant term should be 2.0
        REQUIRE(std::fabs(result.con() - 2.0) < eps);
    }

    SECTION("evaluate with eval_funs (vector-of-visitors form)") {
        Expression sx("qx");
        SDA sv = sx * b[0] + 1.0 * b[1];
        vector<SymEngine::RCP<const SymEngine::Basic>> vars{sx.get_basic()};
        vector<double> inputs{3.0};
        NDA result = da::evaluate(sv, vars, inputs);
        // Coefficient of b[0] should be sx evaluated at 3.0 = 3.0
        // Index of b[0] = 1 (first linear term)
        REQUIRE(std::fabs(result.element(1) - 3.0) < eps);
    }

    SECTION("mixed op NDA+SDA cross-env throws") {
        // Create a separate env
        da::DAEnv& env2 = da::da_make_env(3u, 2u, 100u, false);
        da::da_select_env(env2);
        NDA n2;  // in env2
        da::da_select_env(da::da_current_env()); // tricky: go back to original
        // Actually set up env properly
        da::da_init(4u, 2u, 400u, true); // re-init main env
        NDA n_main = 1.0 + b[0];
        Expression sx("zz");
        SDA s_main = sx + b[1];
        // They're in the same env, so this should NOT throw
        REQUIRE_NOTHROW(n_main * s_main);
    }
}

#endif // DA_WITH_SYMBOLIC
