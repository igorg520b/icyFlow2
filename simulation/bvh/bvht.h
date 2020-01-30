#ifndef BVHT_H
#define BVHT_H

#include "bvhn.h"
#include "geometry/element.h"
#include "SimpleObjectPool.h"
#include <vector>

namespace icy {
class BVHT;
class BVHN;
class kDOP24;
}

class icy::BVHT
{
public:
    static SimpleObjectPool<BVHN> BVHN_Leaf_Factory; // separate list of leaf nodes
    static std::vector<Element*> broad_list; // list of pairs of elements, reuslt of broad phase

    icy::BVHN *root;

    BVHT();
    void Construct(std::vector<Element*> &elems); // construct from a list of surface elems
    void Traverse();
};










#endif // BVHT_H
