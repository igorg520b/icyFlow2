#include <stdexcept>
#include <cfloat>
#include <algorithm>
#include "bvhn.h"

std::vector<icy::Element*> icy::BVHN::broad_list;
icy::AwesomeObjectFactory<std::vector<icy::kDOP24*>> icy::BVHN::kDopVectorFactory(50);
icy::AwesomeObjectFactory<icy::BVHN> icy::BVHN::BVHNFactory(10000);

void icy::BVHN::Initialize(std::vector<kDOP24*> *bvs)
{
    auto count = bvs->size();
    if(count == 0) throw new std::runtime_error("bvs->size is zero in BVHN::Initialize");

    box.Reset();
    if(count == 1) {
        // current node becomes a leaf
        box.Expand(*bvs->front());
        isLeaf = true;
        return;
    }

    isLeaf = false;
    for(auto const &bv : *bvs) box.Expand(*bv); // expand box to the size of bvs collection

    // get length, width and height of the resulting box
    double dX, dY, dZ;
    box.Dimensions(dX, dY, dZ);

    std::vector<icy::kDOP24*> *left = kDopVectorFactory.take();
    std::vector<icy::kDOP24*> *right = kDopVectorFactory.take();
    left->clear();
    right->clear();
    left->reserve(bvs->size());
    right->reserve(bvs->size());

    if (dX >= dY && dX >= dZ)
    {
        double center = box.centerX();
        for(auto const &bv : *bvs) {
            if(bv->centerX() < center) left->push_back(bv);
            else right->push_back(bv);
        }

        // make sure that there is at least one element on each side
        if(left->size() == 0)
        {
            kDOP24 *selected;
            double min1 = DBL_MAX;
            for(auto const &bv : *right)
                if(min1 >= bv->centerX()) {
                    min1 = bv->centerX();
                    selected = bv;
                }
            // move "selected" from left to right
            left->push_back(selected);
            std::remove(right->begin(), right->end(),selected);
        } else if(right->size() == 0)
        {
            kDOP24 *selected;
            double min1 = DBL_MAX;
            for(auto const &bv : *left)
                if(min1 >= bv->centerX()) {
                    min1 = bv->centerX();
                    selected = bv;
                }
            // move "selected" from right to left
            right->push_back(selected);
            std::remove(left->begin(), left->end(),selected);
        }
    }
    else if(dY >= dX && dY >= dZ)
    {
        double center = box.centerY();
        for(auto const &bv : *bvs) {
            if(bv->centerY() < center) left->push_back(bv);
            else right->push_back(bv);
        }

        // make sure that there is at least one element on each side
        if(left->size() == 0)
        {
            kDOP24 *selected;
            double min1 = DBL_MAX;
            for(auto const &bv : *right)
                if(min1 >= bv->centerY()) {
                    min1 = bv->centerY();
                    selected = bv;
                }
            // move "selected" from left to right
            left->push_back(selected);
            std::remove(right->begin(), right->end(),selected);
        }
        else if(right->size() == 0)
        {
            kDOP24 *selected;
            double min1 = DBL_MAX;
            for(auto const &bv : *left)
                if(min1 >= bv->centerY()) {
                    min1 = bv->centerY();
                    selected = bv;
                }
            // move "selected" from right to left
            right->push_back(selected);
            std::remove(left->begin(), left->end(),selected);
        }
    }
    else
    {
        double center = box.centerZ();
        for(auto const &bv : *bvs) {
            if(bv->centerZ() < center) left->push_back(bv);
            else right->push_back(bv);
        }

        // make sure that there is at least one element on each side
        if(left->size() == 0)
        {
            kDOP24 *selected;
            double min1 = DBL_MAX;
            for(auto const &bv : *right)
                if(min1 >= bv->centerZ()) {
                    min1 = bv->centerZ();
                    selected = bv;
                }
            // move "selected" from left to right
            left->push_back(selected);
            std::remove(right->begin(), right->end(),selected);
        }
        else if(right->size() == 0)
        {
            kDOP24 *selected;
            double min1 = DBL_MAX;
            for(auto const &bv : *left)
                if(min1 >= bv->centerZ()) {
                    min1 = bv->centerZ();
                    selected = bv;
                }
            // move "selected" from right to left
            right->push_back(selected);
            std::remove(left->begin(), left->end(),selected);
        }
    }

    child1 = BVHNFactory.take();
    child1->Initialize(left);
    kDopVectorFactory.release(left);

    child2 = BVHNFactory.take();
    child2->Initialize(right);
    kDopVectorFactory.release(right);
}

void icy::BVHN::UpdateLeaf()
{
    if(!isLeaf || elem==nullptr) throw std::runtime_error("BVHN::UpdateLeaf error");
    box.Reset();
    Node *nd = elem->vrts[0];
    box.Expand(nd->tx, nd->ty, nd->tz);
    nd = elem->vrts[1];
    box.Expand(nd->tx, nd->ty, nd->tz);
    nd = elem->vrts[2];
    box.Expand(nd->tx, nd->ty, nd->tz);
    nd = elem->vrts[3];
    box.Expand(nd->tx, nd->ty, nd->tz);
}

void icy::BVHN::Update()
{
    // assume that leaves are already updated
    if(isLeaf) throw std::runtime_error("Updated called on BVHN leaf");
    if(!child1->isLeaf) child1->Update();
    if(!child2->isLeaf) child2->Update();
    box.Reset();
    box.Expand(child1->box);
    box.Expand(child2->box);
}

void icy::BVHN::SelfCollide()
{
    if (isLeaf) return;
    child1->SelfCollide();
    child2->SelfCollide();
    child1->Collide(child2);
}

void icy::BVHN::Collide(BVHN *b)
{
    if(!box.Overlaps(b->box)) return;
    if (this->isLeaf && b->isLeaf)
    {
        broad_list.push_back(elem);
        broad_list.push_back(b->elem);
    }
    else if (this->isLeaf)
    {
        Collide(b->child1);
        Collide(b->child2);
    }
    else
    {
        b->Collide(child1);
        b->Collide(child2);
    }
}

