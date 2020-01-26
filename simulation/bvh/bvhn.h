#ifndef BVHN_H
#define BVHN_H

#include <vector>
#include "geometry/element.h"
#include "kdop24.h"
#include "AwesomeObjectFactory.h"

namespace icy {
class BVHN;
class kDOP24;
}


class icy::BVHN
{
public:
    static std::vector<Element*> broad_list; // list of pairs of elements, reuslt of broad phase
    static AwesomeObjectFactory<std::vector<kDOP24>> kDopFactory(50);

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

#endif // BVHN_H
