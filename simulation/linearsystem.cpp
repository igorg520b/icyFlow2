#include "linearsystem.h"
#include "mkl_pardiso.h"
#include "mkl_types.h"
#include "mkl.h"
#include <stdexcept>
#include <cmath>
#include <cstring>

icy::LinearSystem::LinearSystem()
{
    vals = rhs = dx = nullptr;
}

void icy::LinearSystem::CreateStructure()
{
    csrd.CreateStructure();
    // allocate value arrays
    if (vals_length < dvalsSize()) {
        delete vals;
        vals_length = dvalsSize()*2;
        vals = new double[vals_length];
    }

    if(dx_length < dxSize()) {
        delete rhs;
        delete dx;
        dx_length = dxSize();
        rhs = new double[dx_length];
        dx = new double[dx_length];
    }
    memset(rhs, 0, sizeof(double)*dxSize());
    memset(dx, 0, sizeof(double)*dxSize());
    memset(vals, 0, sizeof(double)*dvalsSize());
    Assert();
}

void icy::LinearSystem::Solve()
{
    const int mklMatrixType = -2; // -2 for symmetric indefinite
    const int check = 0;
    const int verbosity = 0;
    const int dim = 3;
    const int param4 = 0;

    int mklResult = SolveDouble3(csrd.csr_cols, csrd.csr_rows,
                                 vals, csrd.N, rhs, dx,
                                 mklMatrixType, param4, dim, verbosity, check);
    if(mklResult != 0) throw std::runtime_error("MKL solver error");
}

double icy::LinearSystem::NormOfDx()
{
    double result = 0;
    for (int i = 0; i < dxSize(); i++) result += dx[i] * dx[i];
    return result;
}

void icy::LinearSystem::AddToRHS(int atWhichIndex, double d0, double d1, double d2)
{
    if (atWhichIndex < 0) return;
//    if(std::isnan(d0)) throw std::runtime_error("d0 is nan");
//    if(std::isnan(d1)) throw std::runtime_error("d1 is nan");
//    if(std::isnan(d2)) throw std::runtime_error("d2 is nan");

    int i3 = atWhichIndex*3;

    if(i3 + 2 >= csrd.N * 3)
        throw std::runtime_error("i3 + 2 >= csrd.N * 3; index out of range");
    rhs[i3] += d0;
    rhs[i3 + 1] += d1;
    rhs[i3 + 2] += d2;
}

void icy::LinearSystem::AddToLHS_Symmetric(int row, int column,
        double a00, double a01, double a02,
        double a10, double a11, double a12,
        double a20, double a21, double a22)
{
    if (row > column || row < 0 || column < 0) return;
/*
    if(std::isnan(a00) ||
            std::isnan(a01) ||
            std::isnan(a02) ||
            std::isnan(a10) ||
            std::isnan(a11) ||
            std::isnan(a12) ||
            std::isnan(a20) ||
            std::isnan(a21) ||
            std::isnan(a22))
        throw std::runtime_error("lhs element is nan");
*/

    int offset = csrd.offset(row, column);
    if(offset >= csrd.nnz) throw std::runtime_error("offset >=csrd.nnz");
    offset *= 9;
    vals[offset + 0] += a00;
    vals[offset + 1] += a01;
    vals[offset + 2] += a02;
    vals[offset + 3] += a10;
    vals[offset + 4] += a11;
    vals[offset + 5] += a12;
    vals[offset + 6] += a20;
    vals[offset + 7] += a21;
    vals[offset + 8] += a22;
}

void icy::LinearSystem::Assert()
{
    for(int i=0;i<dxSize();i++)
        if(std::isnan(rhs[i])) throw std::runtime_error("rhs constains NaN");
    for(int i=0;i<dvalsSize();i++)
        if(std::isnan(vals[i])) throw std::runtime_error("matrix contains NaN");
    csrd.Assert();
}



// ja = rows; ia = cols; a = vals
int icy::LinearSystem::SolveDouble3(int *ja, int *ia, double *a,
    const int n, double *b, double *x,
                 int _mtype,
                 int iparam4,
                 int dim,
                 int msglvl_,
                 int check) {

    MKL_INT mtype = _mtype;       // Real symmetric matrix: -2;  real unsymmetric: 11
    MKL_INT nrhs = 1;     // Number of right hand sides.
    void *pt[64];
    MKL_INT iparm[64] = {};
    MKL_INT maxfct, mnum, phase, error, msglvl;
    MKL_INT i;
    MKL_INT idum;
    iparm[0] = 1;         // No solver default
    iparm[1] = 3;         // Fill-in reordering from METIS (was 2)
    iparm[3] = iparam4;         // No iterative-direct algorithm
    iparm[4] = 0;         // No user fill-in reducing permutation
    iparm[5] = 0;         // Write solution into x
    iparm[6] = 0;         // Not in use
    iparm[7] = 0;         // Max numbers of iterative refinement steps
    iparm[8] = 0;
    iparm[9] = 8;        // Perturb the pivot elements with 1E-13
    iparm[10] = 0;        // Use nonsymmetric permutation and scaling MPS
    iparm[11] = 0;
    iparm[12] = 0;        // Maximum weighted matching algorithm is switched-off (default for symmetric). Try iparm[12] = 1 in case of inappropriate accuracy
    iparm[13] = 0;        // Output: Number of perturbed pivots
    iparm[14] = 0;
    iparm[15] = 0;
    iparm[16] = 0;
    iparm[17] = 1;       // 1 - disable report; Output: Number of nonzeros in the factor LU
    iparm[18] = 1;			// 1- disable report; output number of operations
    iparm[19] = 0;
    iparm[26] = check; // check matrix structure for errors
    iparm[27] = 0; // 0 double; 1 single
    iparm[34] = 1;// zero-base index
    iparm[36] = dim;// BSR with block size 3
    maxfct = 1;
    mnum = 1;
    msglvl = msglvl_; // use 1 for verbose output
    error = 0;
    for (i = 0; i < 64; i++) pt[i] = 0;
    phase = 13;
    PARDISO(pt, &maxfct, &mnum, &mtype, &phase, &n, a, ia, ja, &idum, &nrhs, iparm, &msglvl, b, x, &error);

    // clean up (?)
    /*
    phase = -1;
    double ddum;
    PARDISO(pt, &maxfct, &mnum, &mtype, &phase, &n, &ddum, ia, ja, &idum, &nrhs, iparm, &msglvl, &ddum, &ddum, &error);
*/
    return error;
}
