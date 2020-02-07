#include "mainwindow.h"

#include <QApplication>
#include <QSurfaceFormat>

#include "mkl_pardiso.h"
#include "mkl_types.h"
#include<chrono>
#include<iostream>
using namespace std;

// ja = rows; ia = cols; a = vals
int solveMKL_BCSR(int *ja, int *ia, double *a,
        const int n, double *b, double *x, int _mtype,
        int iparam4, int dim, int msglvl_, int check, int zeroIdx=1)
{
    auto t1 = std::chrono::high_resolution_clock::now();
    //	mkl_verbose(0);
    //	std::stringstream redirectStream;
    //	std::cout.rdbuf(redirectStream.rdbuf());

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
    iparm[34] = zeroIdx;// zero-base index
    iparm[36] = dim;// BSR with block size 3
    maxfct = 1;
    mnum = 1;
    msglvl = msglvl_; // use 1 for verbose output
    error = 0;
    for (i = 0; i < 64; i++) pt[i] = 0;
    phase = 13;
    PARDISO(pt, &maxfct, &mnum, &mtype, &phase, &n, a, ia, ja, &idum, &nrhs, iparm, &msglvl, b, x, &error);

    phase = -1;
    double ddum;
    PARDISO(pt, &maxfct, &mnum, &mtype, &phase, &n, &ddum, ia, ja, &idum, &nrhs, iparm, &msglvl, &ddum, &ddum, &error);

    auto t2 = std::chrono::high_resolution_clock::now();

    auto duration = std::chrono::duration_cast<std::chrono::microseconds>( t2 - t1 ).count();

    duration/=1000;
    return error;
}


int main(int argc, char *argv[])
{
    QApplication a(argc, argv);

    QSurfaceFormat fmt = QVTKOpenGLNativeWidget::defaultFormat();
    fmt.setAlphaBufferSize(0);
    QSurfaceFormat::setDefaultFormat(fmt);

    MainWindow w;

    cout << "start " << endl;
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

    double x[12]={};
    int result = solveMKL_BCSR(test_cols,
            test_rows, test_lhs, 4, test_rhs, x, -2,0,3,1,0,1);
    cout << "result " << result << endl;

    w.showMaximized();
//    w.show();
    return a.exec();
}
