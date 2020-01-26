#ifndef BVHT_H
#define BVHT_H

#include "geometry/element.h"
#include <vector>

namespace icy {
class BVHT;
class BVHN;
class kDOP24;
}

class icy::BVHT
{
public:
    BVHT();
};

class icy::kDOP24
{
public:
    double d0, d1, d2, d3, d4, d5, d6, d7, d8, d9, d10, d11; // lower boundaries
    double d12, d13, d14, d15, d16, d17, d18, d19, d20, d21, d22, d23; // higher boundaries
    Element *elem = nullptr; // element that is enveloped by this kDOP, if exists
    kDOP24();
    kDOP24(kDOP24 &k);
    void Reset();
    void UpdateTentative(Element &e);
    bool Overlaps(kDOP24 &b);
    void Expand(double x, double y, double z);
    void Expand(kDOP24 &b);
    void Dimensions(double &dx, double &dy, double &dz);

private:
    void MinMax(double p, double &mi, double &ma);
    void MinMax(double a, double b, double &mi, double &ma);
    double min(double a, double b) { return a < b ? a : b; }
    double max(double a, double b) { return a > b ? a : b; }
};



class icy::BVHN
{
public:
    static std::vector<Element*> broad_list; // list of pairs of elements, reuslt of broad phase

    kDOP24 box;
    BVHN *child1, *child2, *parent;
    int level;
    bool isLeaf;

    BVHN(BVHN &parent, std::vector<kDOP24> *bvs, int level);
    void Initialize(BVHN &parent, std::vector<kDOP24> *bvs, int level);
    void FinalizeConstruction();
    void Update();
    void SelfCollide();
    void Collide(BVHN &b);

private:
    std::vector<kDOP24> *bvs = nullptr;
};





#endif // BVHT_H
