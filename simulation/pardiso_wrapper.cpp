#include<iostream>
#include<sstream>

#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include "GteSymmetricEigensolver3x3.h"

using namespace gte;






// solver
SymmetricEigensolver3x3<double> solver;
void Eigenvalues(double xx, double yy, double zz,
    double xy, double yz, double zx,
    double* eigenvalues) {
    std::array<double, 3> eval;
    std::array<std::array<double, 3>, 3> evec;
    solver(xx, xy, zx, yy, yz, zz, false, -1, eval, evec);
    eigenvalues[0] = eval[0];
    eigenvalues[1] = eval[1];
    eigenvalues[2] = eval[2];
}
