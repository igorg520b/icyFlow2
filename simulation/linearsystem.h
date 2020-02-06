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
    double NormOfRHS();
    double NormOfLHS();
    void AddToRHS(int atWhichIndex, double d0, double d1, double d2);
    void AddToLHS_Symmetric(int row, int column,
            double a00, double a01, double a02,
            double a10, double a11, double a12,
            double a20, double a21, double a22);
    void Assert();
    int dvalsSize() { return csrd.nnz*9; }
    int dxSize() { return csrd.N*3; }
    void printout(); //for testing
    void testSolve(); // use test data (below)
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

    // testing
    int test_rows[5] ={0, 2, 4, 5, 6};
    int test_cols[6] = {0, 3, 1, 2, 2, 3};
    double test_rhs[12] = {14876.1, 95200.2, 111749, 94820.2, -171439, -217089, 194609, 418886, 729953, -130703, -162140, 255506};
    double test_lhs[54] = {7.99789e+07, -1.50543e+07, 2.51794e+07, -1.50543e+07, 3.18622e+07, -7.22864e+06,
                           2.51794e+07, -7.22864e+06, 3.96307e+07, -8.60798e+06, 3.31117e+07, -1.73241e+07,
                           2.18364e+07, -2.70753e+07, 1.58016e+07, -1.11512e+07, 1.94435e+07, -2.5342e+07,
                           1.42278e+08, -6.09242e+06, -6.72327e+07, -6.09242e+06, 5.61172e+07, 4.73046e+06,
                           -6.72327e+07, 4.73046e+06, 1.07891e+08, -7.22231e+07, 3.55869e+07, 6.99779e+07,
                           2.4375e+07, -4.85718e+07, -2.23346e+07, 5.38297e+07, -2.99038e+07, -1.23901e+08,
                           7.63008e+07, -1.85796e+07, -3.27563e+07, -1.85796e+07, 1.06283e+08, 6.85186e+07,
                           -3.27563e+07, 6.85186e+07, 1.88219e+08, 3.05137e+07, 2.64046e+06, -1.40745e+06,
                           2.64046e+06, 8.94798e+07, -3.14937e+07, -1.40745e+06, -3.14937e+07, 4.71829e+07};
};

#endif // LINEARSYSTEM_H
