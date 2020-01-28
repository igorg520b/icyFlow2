#ifndef LINEARSYSTEM_H
#define LINEARSYSTEM_H

#include "csrdictionary.h"

namespace icy {
class LinearSystem;
}

class icy::LinearSystem
{
public:
    CSRDictionary csrd;
    double *vals, *rhs, *dx;   // non-zero values, right-hand side and solution

    LinearSystem();
    void CreateStructure();
    void Solve();
    double NormOfDx();
    void AddToRHS(int atWhichIndex, double d0, double d1, double d2);
    void AddToLHS_Symmetric(int row, int column,
            double a00, double a01, double a02,
            double a10, double a11, double a12,
            double a20, double a21, double a22);
    void Assert();
    int dvalsSize() { return csrd.nnz*9; }
    int dxSize() { return csrd.N*3; }
private:
    int SolveDouble3(int *ja, int *ia, double *a,
        const int n, double *b, double *x,
                     int _mtype = -2,
                     int iparam4 = 0,
                     int dim = 3,
                     int msglvl_ = 0,
                     int check = 0);
    int dx_length = 0; // number of currently allocated elements for dx and rhs
    int vals_length = 0; // currently allocated vals
};

#endif // LINEARSYSTEM_H
