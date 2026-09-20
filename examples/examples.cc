/**
 * @file examples.cc
 * @brief Port of ref/tpsa/examples/examples.cc to the new da:: API.
 *
 * Changes:
 *   - #include "da.h" -> #include "da/da.h"
 *   - DAVector -> da::NDA  (or using da::NDA; with using namespace da;)
 *   - da[i]    -> da::base[i]
 *   - da_init, da_clear, etc. -> da:: namespace
 *   - DAVector::eps -> da::NDA::eps
 */
#include "da/da.h"
#include <cmath>
#include <fstream>
#include <iostream>
#include <vector>
#include <complex>

using namespace std::complex_literals;
using std::complex;
using da::NDA;

int main() {
    unsigned int da_dim   = 3;
    unsigned int da_order = 4;
    unsigned int n_vec    = 100;

    da::da_init(da_order, da_dim, n_vec);

    {  // scope: all DA vectors destroyed here, before da_clear()
    std::cout << "Print out the base vector." << std::endl << std::endl;
    da::base[0].print();
    da::base[1].print();
    da::base[2].print();

    std::cout << "Fundamental calculations of DA vectors." << std::endl << std::endl;
    NDA x = NDA(1.0) + da::base[0] + NDA(2.0)*da::base[1] + NDA(5.0)*da::base[2];
    x.print();

    NDA y = da::exp(x);
    y.print();

    // Substitute a number for a base.
    std::cout << "Substitute a number for a base." << std::endl << std::endl;
    NDA z;
    da::da_substitute(y, 0, 1.0, z);
    z.print();

    // Substitute a DA vector for a base.
    std::cout << "Substitute a DA vector for a base." << std::endl << std::endl;
    da::da_substitute(y, 0, x, z);
    z.print();

    // Substitute multiple DA vectors for bases at once.
    std::cout << "Substitute multiple DA vectors for bases at once." << std::endl << std::endl;
    std::vector<NDA> lv(2);
    lv[0] = da::sin(x);
    lv[1] = da::cos(x);
    std::vector<unsigned int> idx{0, 1};
    da::da_substitute(y, idx, lv, z);
    z.print();

    std::cout << "The norm of z is " << z.norm() << std::endl;
    std::cout << "The weighted norm of z is " << z.weighted_norm(0.1) << std::endl;

    // Bunch processing for substitutions.
    std::cout << "Bunch processing for substitutions." << std::endl << std::endl;
    std::vector<NDA> lx(3);
    std::vector<NDA> ly(3);
    lx[0] = x;
    lx[1] = y;
    lx[2] = da::sinh(x);

    da::da_substitute(lx, idx, lv, ly);
    ly[0].print();
    ly[1].print();
    ly[2].print();

    // Composition of DA vectors with numbers.
    std::cout << "Composition of DA vectors with numbers." << std::endl << std::endl;
    std::vector<double> lm{0.1, 2, 1};
    std::vector<double> ln(3);
    da::da_composition(lx, lm, ln);
    for (auto& val : ln) std::cout << val << ' ';
    std::cout << std::endl << std::endl;

    // Composition of DA vectors with DA vectors.
    std::cout << "Composition of DA vectors with DA vectors." << std::endl << std::endl;
    std::vector<NDA> lu(3);
    lu[0] = da::sin(x);
    lu[1] = da::cos(x);
    lu[2] = da::tan(x);
    da::da_composition(lx, lu, ly);
    ly[0].print();
    ly[1].print();
    ly[2].print();

    // Output a DA vector to file.
    std::ofstream fout;
    fout.open("da_output.txt");
    fout << z;
    fout.close();

    z /= 1e5;
    std::cout << "z:" << std::endl;
    z.print();
    std::cout << "norm of z: " << z.norm() << std::endl;

    std::cout << "DA eps: " << NDA::eps << std::endl;
    da::da_set_eps(1e-20);
    std::cout << "Reset DA eps: " << NDA::eps << std::endl;

    }  // end of DA scope: all NDA variables destroyed before da_clear()
    da::da_clear();
    return 0;
}
