#ifndef NUMBERCRUNCHING_H
#define NUMBERCRUNCHING_H

#include <omp.h>
#include <cstring>
#include "geometry/meshcollection.h"
#include "geometry/element.h"
#include "simulation/linearsystem.h"
#include "modelprms.h"
#include <vector>
#include <unordered_set>
#include "GteSymmetricEigensolver3x3.h"
using namespace gte;

namespace icy {
class NumberCrunching;
class CPResult;
class ModelPrms;
}

class icy::NumberCrunching
{
public:

    // narrow phase
    static const double EPS;
    static std::vector<CPResult> cprList;

    static void InitializeConstants();

    static void NarrowPhase(std::vector<Element*> &broadList, MeshCollection &mc);

    // collision response (collision data stored in cprList)
    static void CollisionResponse(LinearSystem &ls, double DistanceEpsilon, double k);

    // linear tetrahedron
    static void AssembleElems(LinearSystem &ls, std::vector<Element*> &elasticElements, ModelPrms &prms, double h);

    // cohesive zones
    static void AssembleCZs(
            LinearSystem &ls, std::vector<CZ*> &czs, ModelPrms &prms,
            int &totalFailed, int &totalDamaged);

private:

    // narrow phase

    static std::vector<int> resultingList; // results of NarrowPhaseTwoElems
    static std::unordered_set<long long> NL2set;       // Tuple ( node# inside element, which element)
    static std::vector<long long> NL2vector;
    static SymmetricEigensolver3x3<double> solver;

    inline static void Bvalues(double x0, double y0, double z0,
        double x1, double y1, double z1,
        double x2, double y2, double z2,
        double x3, double y3, double z3,
        double &b11, double &b12, double &b13,
        double &b21, double &b22, double &b23,
        double &b31, double &b32, double &b33);

    inline static bool ctest(double b11, double b12, double b13,
        double b21, double b22, double b23,
        double b31, double b32, double b33,
        double x, double y, double z);

    static int NarrowPhaseTwoElems(Element *tetra1, Element *tetra2);
    static Face* FindClosestFace(Node *nd, Element *elem);

    // distance triangle-to-node
    static double dtn(
                double f1x, double f1y, double f1z,
                double f2x, double f2y, double f2z,
                double f3x, double f3y, double f3z,
                double ndx, double ndy, double ndz);

    static double DOT(double (&v1)[3], double (&v2)[3]) { return (v1[0] * v2[0] + v1[1] * v2[1] + v1[2] * v2[2]); }

    static void SUB(double (&dest)[3], double (&v1)[3], double (&v2)[3])
    {
        dest[0] = v1[0] - v2[0];
        dest[1] = v1[1] - v2[1];
        dest[2] = v1[2] - v2[2];
    }

    static double clamp(double n) { return n <= 0 ? 0 : n >= 1 ? 1 : n; }
    // collision response
    static void OneCollision(double distanceEpsilonSqared, double k, CPResult &res);

    //linear tetrahedron

    // the results of this function subsequently go into the equaiton of motion
    // f[12] = elastic forces acting on nodes
    // Df[12][12] = df/dx
    // V = tetrahedron rest volume
    static void F_and_Df_Corotational(
            const double(&x0)[12], const double(&xc)[12],
    double(&f)[12], double(&Df)[12][12], double &V,
    double(&sigma)[6], double (&principal_str)[3],
    const double (&E)[6][6]);

    static inline void fastRotationMatrix(
            double p0x, double p0y, double p0z,
            double p1x, double p1y, double p1z,
            double p2x, double p2y, double p2z,
            double &r11, double &r12, double &r13,
            double &r21, double &r22, double &r23,
            double &r31, double &r32, double &r33);

    static inline void multABd(
            double a11, double a12, double a13,
            double a21, double a22, double a23,
            double a31, double a32, double a33,
            double b11, double b12, double b13,
            double b21, double b22, double b23,
            double b31, double b32, double b33,
            double &m11, double &m12, double &m13,
            double &m21, double &m22, double &m23,
            double &m31, double &m32, double &m33);

    static inline void multAX(
            double a11, double a12, double a13,
            double a21, double a22, double a23,
            double a31, double a32, double a33,
            double x1, double x2, double x3,
            double &y1, double &y2, double &y3);

    static void Eigenvalues(double xx, double yy, double zz,
        double xy, double yz, double zx,
        double* eigenvalues);

    static void ElementElasticity(
            Element *elem,
            const double (&E)[6][6], const double rho,
    const double dampingMass, const double dampingStiffness, const double h,
    const double NewmarkBeta, const double NewmarkGamma, const double (&M)[12][12]);


    // cohesive zones
    // these are copied from ModelPrms
    static double deln, delt, p_m, p_n, alpha, beta, gam_n, gam_t, tau_max, sigma_max, pMtn, pMnt, lambda_n, lambda_t;
    static double B[3][3][18];
    static double sf[3][3];

    inline static double Tn_(const double Dn, const double Dt);
    inline static double Tt_(const double Dn, const double Dt);
    inline static double Dnn_(const double opn, const double opt);
    inline static double Dtt_(const double opn, const double opt);
    inline static double Dnt_(const double opn, const double opt);

    inline static void cohesive_law(
            const double opn, const double opt,
            bool &cz_contact, bool &cz_failed,
            double &pmax, double &tmax,
            double &Tn, double &Tt, double &Dnn,
            double &Dtt, double &Dnt, double &Dtn);

    inline static void CZRotationMatrix(
        double x0, double y0, double z0,
        double x1, double y1, double z1,
        double x2, double y2, double z2,
        double &r00, double &r01, double &r02,
        double &r10, double &r11, double &r12,
        double &r20, double &r21, double &r22,
        double &a_Jacob);

    static void CZForce(CZ *cz);

};


class icy::CPResult
{
public:
    Node *nd = nullptr;
    Face *fc = nullptr;
    double fi[12] = {};
    double dfi[12][12] = {};
    int idxs[4] = {};
    Node *nds[4] = {};

    CPResult() {}

    void Clear()
    {
        memset(fi, 0, sizeof(double)*12);
        memset(dfi, 0, sizeof(double)*12*12);
        //memset(idxs, 0, sizeof(int)*4);
    }
};

#endif // NUMBERCRUNCHING_H
