#ifndef BVHN_H
#define BVHN_H

#include <vector>
#include "geometry/element.h"
#include "kdop24.h"
#include "SimpleObjectPool.h"
#include "bvht.h"

namespace icy {
class BVHN;
class BVHT;
class kDOP24;
}


class icy::BVHN
{
public:
    static SimpleObjectPool<std::vector<BVHN*>> VectorFactory;
    static SimpleObjectPool<BVHN> BVHNFactory;

    kDOP24 box;
    bool isLeaf;
    int level;
    Element *elem = nullptr; // element that is enveloped by this kDOP, if leaf

    BVHN();
    void Initialize(std::vector<BVHN*> *bvs, int level_);
    void Update();
    void UpdateLeaf();  // use tentative coordinates from elem
    void SelfCollide();
    void Collide(BVHN *b);

private:
    BVHN *child1, *child2;
};

#endif // BVHN_H
