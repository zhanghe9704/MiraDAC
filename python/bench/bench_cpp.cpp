// C++ baseline for the Python binding benchmark (plan T0.5).
//
// Times one op per case on dense random DA vectors. Each case first doubles
// its repetition count until one sample takes >= 10 ms (the warm-up). Then
// kRounds rounds each give every case 0.2 s of samples (at least one), and a
// case reports its minimum sample. Interleaving spreads a case's samples over
// the whole run, so a burst of load from other processes hits some samples of
// every case rather than all of one case; short samples make it likely that
// some are clean, and the minimum keeps those, since load only ever adds
// time. Times are thread CPU time. On Linux, transparent huge pages are off
// for the process; with them on, small cases moved by several percent
// between runs.
// Prints {"case": ns_per_op, ...} as JSON on stdout. With --symbolic, runs
// only the symbolic cases (for comparing SymEngine builds); with --numeric,
// only the numeric ones. With --driven, runs one round per line read from
// stdin and answers each with "done", printing the JSON at end of input;
// bench_ops.py uses this to interleave the C++ rounds with its own.

#include "da/da.h"

#include <algorithm>
#include <cstdio>
#include <ctime>
#include <iostream>
#include <limits>
#include <map>
#include <random>
#include <string>
#include <utility>
#include <vector>

#ifdef __linux__
#include <sys/prctl.h>
#endif

namespace {

constexpr double kSample = 0.01;  // seconds per sample, at least
constexpr double kPerRound = 0.2;  // seconds of samples per case and round
constexpr int kRounds = 15;

std::vector<std::string> names;  // in first-run order
std::map<std::string, long> reps_of;    // set by the warm-up
std::map<std::string, double> best_ns;  // minimum sample

double cpu_now() {
    timespec ts;
    clock_gettime(CLOCK_THREAD_CPUTIME_ID, &ts);
    return ts.tv_sec + 1e-9 * ts.tv_nsec;
}

template <class F>
double seconds(F& op, long reps) {
    double t0 = cpu_now();
    for (long i = 0; i < reps; ++i) op();
    return cpu_now() - t0;
}

template <class F>
void bench(const std::string& name, F op) {
    auto it = reps_of.find(name);
    if (it == reps_of.end()) {
        long reps = 1;
        while (seconds(op, reps) < kSample) reps *= 2;
        it = reps_of.emplace(name, reps).first;
        names.push_back(name);
        best_ns[name] = std::numeric_limits<double>::infinity();
    }
    long reps = it->second;
    double& best = best_ns[name];
    for (double spent = 0; spent < kPerRound;) {
        double t = seconds(op, reps);
        best = std::min(best, t * 1e9 / reps);
        spent += t;
    }
}

std::uniform_real_distribution<double> uni(-1.0, 1.0);

// Same seed every round, so every round times the same data.
da::NDA random_nda(std::mt19937& rng) {
    da::NDA v;
    for (int i = 0; i < da::da_full_length(); ++i)
        v.set_element(da::da_element_orders(i), uni(rng));
    return v;
}

void numeric(unsigned nvars, unsigned order) {
    da::da_init(order, nvars, 400, true);
    {
        std::mt19937 rng(12345);
        std::string p = "n" + std::to_string(nvars) + "o" + std::to_string(order) + "/";
        da::NDA a = random_nda(rng), b = random_nda(rng), c;
        bench(p + "add", [&] { c = a + b; });
        bench(p + "mul", [&] { c = a * b; });
        bench(p + "iadd", [&] { c += b; });
        bench(p + "mul_const", [&] { c = a * 2.0; });
        bench(p + "exp", [&] { c = da::exp(a); });

        da::CNDA ca(a, b), cb(b, a), cc;
        bench(p + "cmul", [&] { cc = ca * cb; });
        bench(p + "cexp", [&] { cc = da::exp(ca); });

        std::vector<da::NDA> m, n, out(nvars);
        for (unsigned i = 0; i < nvars; ++i) {
            m.push_back(random_nda(rng));
            n.push_back(random_nda(rng));
        }
        bench(p + "composition", [&] { da::da_composition(m, n, out); });
    }
    da::da_clear();
}

void symbolic(unsigned nvars, unsigned order) {
    using SymEngine::Expression;
    da::da_init(order, nvars, 400, true);
    {
        // Every coefficient is its own symbol.
        auto random_sda = [](const char* prefix) {
            da::SDA v;
            for (int i = 0; i < da::da_full_length(); ++i)
                v.set_element(i, Expression(prefix + std::to_string(i)));
            return v;
        };
        std::string p = "sym_n" + std::to_string(nvars) + "o" + std::to_string(order) + "/";
        da::SDA a = random_sda("a"), b = random_sda("b"), c;
        bench(p + "mul", [&] { c = a * b; });
        bench(p + "exp", [&] { c = da::exp(a); });
    }
    da::da_clear();
}

}  // namespace

int main(int argc, char** argv) {
#ifdef __linux__
    prctl(PR_SET_THP_DISABLE, 1, 0, 0, 0);
#endif
    std::vector<std::string> args(argv + 1, argv + argc);
    auto has = [&](const char* flag) {
        return std::find(args.begin(), args.end(), flag) != args.end();
    };
    bool numeric_cases = !has("--symbolic");
    bool symbolic_cases = !has("--numeric");
    auto round = [&] {
        if (numeric_cases)
            for (auto [nvars, order] : {std::pair{3u, 4u}, {6u, 6u}, {6u, 10u}})
                numeric(nvars, order);
        if (symbolic_cases) symbolic(3, 3);
    };
    if (has("--driven")) {
        for (std::string line; std::getline(std::cin, line);) {
            round();
            std::printf("done\n");
            std::fflush(stdout);
        }
    } else {
        for (int r = 0; r < kRounds; ++r) {
            round();
            std::fprintf(stderr, "round %d/%d done\n", r + 1, kRounds);
        }
    }

    std::printf("{\n");
    for (std::size_t i = 0; i < names.size(); ++i)
        std::printf("  \"%s\": %.1f%s\n", names[i].c_str(), best_ns[names[i]],
                    i + 1 < names.size() ? "," : "");
    std::printf("}\n");
}
