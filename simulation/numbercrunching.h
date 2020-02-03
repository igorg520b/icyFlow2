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

    static void NarrowPhase(std::vector<Element*> &broadList, MeshCollection &mc,
        std::vector<CPResult> &cprList);

    // collision response
    static void CollisionResponse(LinearSystem &ls, std::vector<CPResult> &cprList,
                                  double DistanceEpsilon, double k);

private:

    // narrow phase

    static std::vector<int> resultingList; // results of NarrowPhaseTwoElems
    static std::unordered_set<long long> NL2set;       // Tuple ( node# inside element, which element)
    static std::vector<long long> NL2vector;

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
    static double dtn(
                double f1x, double f1y, double f1z,
                double f2x, double f2y, double f2z,
                double f3x, double f3y, double f3z,
                double ndx, double ndy, double ndz);

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
