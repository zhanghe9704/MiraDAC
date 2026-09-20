/**
 * @file example_complex_da.cc
 * @brief Port of ref/tpsa/examples/example_complex_da.cc to the new da:: API.
 */
#include "da/da.h"
#include <cmath>
#include <iostream>
#include <vector>
#include <complex>

using namespace std::complex_literals;
using std::complex;
using da::NDA;

int main() {
    unsigned int da_dim   = 3;
    unsigned int da_order = 4;
    unsigned int n_vec    = 400;

    da::da_init(da_order, da_dim, n_vec);

    {  // scope: all DA vectors destroyed before da_clear()
    NDA x1, x2, x3, x4;
    x1 = da::base[0] + NDA(2.0)*da::base[1] + NDA(3.0)*da::base[2];
    x2 = da::sin(x1);
    x1 = da::cos(x1);

    x3 = NDA(0.5)*da::base[0] + NDA(4.0)*da::base[1] + NDA(2.7)*da::base[2];
    x4 = da::sin(x3);
    x3 = da::cos(x3);

    auto y1 = x1 + x2 * 1i;
    auto y2 = x3 + x4 * 1i;

    std::cout << "y1: " << std::endl << y1 << std::endl;
    std::cout << "y2: " << std::endl << y2 << std::endl;

    std::cout << "y1+y2: " << std::endl << y1 + y2 << std::endl;
    std::cout << "y1-y2: " << std::endl << y1 - y2 << std::endl;
    std::cout << "y1*y2: " << std::endl << y1 * y2 << std::endl;
    std::cout << "y1/y2: " << std::endl << y1 / y2 << std::endl;

    std::vector<NDA> mmap;
    mmap.push_back(x1);
    mmap.push_back(x2);

    std::vector<std::complex<double>> nmap;
    nmap.push_back(4.2 + 0.3*1i);
    nmap.push_back(1/3.0 + std::sqrt(2.0)*1i);
    nmap.push_back(std::sin(0.7) + std::cos(0.4)*1i);

    std::vector<std::complex<double>> omap(2);
    da::da_composition(mmap, nmap, omap);
    std::cout << "Composition of DA vectors with complex numbers." << std::endl << std::endl;
    for (auto& o : omap) std::cout << o << std::endl;

    std::vector<complex<NDA>> cnmap;
    cnmap.push_back(y1);
    cnmap.push_back(y2);
    cnmap.push_back(y1 * y2);

    std::vector<complex<NDA>> comap(2);
    da::cd_composition(mmap, cnmap, comap);
    std::cout << "Composition of DA vectors with complex DA vectors." << std::endl << std::endl;
    for (auto& o : comap) std::cout << o << std::endl;

    std::vector<complex<NDA>> cmmap;
    cmmap.push_back(x1 + 1i * da::exp(x1));
    cmmap.push_back(x2 + 1i * da::exp(x2));

    da::cd_composition(cmmap, cnmap, comap);
    std::cout << "Composition of complex DA vectors with complex DA vectors." << std::endl << std::endl;
    for (auto& o : comap) std::cout << o << std::endl;

    mmap.push_back(x1 + NDA(0.33)*x2);
    da::cd_composition(cmmap, mmap, comap);
    std::cout << "Composition of complex DA vectors with DA vectors." << std::endl << std::endl;
    for (auto& o : comap) std::cout << o << std::endl;

    }  // end of DA scope
    da::da_clear();
    return 0;
}
