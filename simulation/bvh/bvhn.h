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
    static AwesomeObjectFactory<std::vector<kDOP24*>> kDopVectorFactory;
    static AwesomeObjectFactory<BVHN> BVHNFactory;

    kDOP24 box;
    bool isLeaf;
    Element *elem = nullptr; // element that is enveloped by this kDOP, if leaf

    void Initialize(std::vector<kDOP24*> *bvs);
    void Update();
    void UpdateLeaf();  // use tentative coordinates from elem
    void SelfCollide();
    void Collide(BVHN *b);

private:
    BVHN *child1, *child2;
};

#endif // BVHN_H
