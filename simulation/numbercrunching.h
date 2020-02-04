#ifndef NUMBERCRUNCHING_H
#define NUMBERCRUNCHING_H

#include <cstring>
#include "geometry/meshcollection.h"
#include "geometry/element.h"
#include "simulation/linearsystem.h"
#include <vector>
#include <unordered_set>

namespace icy {
class NumberCrunching;
class CPResult;
}

class icy::NumberCrunching
{
public:
    NumberCrunching();

    // narrow phase
    static const double EPS;

    static void NarrowPhase(std::vector<Element*> &broadList, MeshCollection &mc);

    // collision response (collision data stored in cprList)
    static void CollisionResponse(LinearSystem &ls,
                                  double DistanceEpsilon, double k);

private:

    // narrow phase

    static std::vector<int> resultingList; // results of NarrowPhaseTwoElems
    static std::unordered_set<long long> NL2set;       // Tuple ( node# inside element, which element)
    static std::vector<long long> NL2vector;
    static std::vector<CPResult> cprList;

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
