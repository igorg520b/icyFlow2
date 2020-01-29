//#include <limits>
#include "bvht.h"

std::vector<icy::Element*> icy::BVHT::broad_list; // populated during traversal with pairs
icy::SimpleObjectPool<icy::BVHN> icy::BVHT::BVHN_Leaf_Factory(10000);

icy::BVHT::BVHT()
{
    broad_list.reserve(2000000);
    root = new BVHN();
}

void icy::BVHT::Construct(std::vector<Element*> elems)
{
    BVHN_Leaf_Factory.releaseAll();
    BVHN::BVHNFactory.releaseAll();
    BVHN::VectorFactory.releaseAll();

    // construct boundary volume hierarchy from a list of elements
    std::vector<BVHN*> *root_vec = BVHN::VectorFactory.take();

    for(auto const &elem : elems)
    {
        // initialize from element
        BVHN *leafNode = BVHN_Leaf_Factory.take();
        leafNode->isLeaf = true;

        kDOP24 &box = leafNode->box;
        box.Reset();
        Node *nd = elem->vrts[0];
        box.Expand(nd->tx,nd->ty,nd->tz);
        nd = elem->vrts[1];
        box.Expand(nd->tx,nd->ty,nd->tz);
        nd = elem->vrts[2];
        box.Expand(nd->tx,nd->ty,nd->tz);
        nd = elem->vrts[3];
        box.Expand(nd->tx,nd->ty,nd->tz);
    }
    root->Initialize(root_vec);

    BVHN::VectorFactory.release(root_vec);
}

void icy::BVHT::Traverse()
{
    broad_list.clear();
    root->SelfCollide();
}




