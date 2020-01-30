#include <stdexcept>
#include <cfloat>
#include <algorithm>
#include "bvhn.h"

icy::SimpleObjectPool<std::vector<icy::BVHN*>> icy::BVHN::VectorFactory(50);
icy::SimpleObjectPool<icy::BVHN> icy::BVHN::BVHNFactory(50000);

icy::BVHN::BVHN() {}

void icy::BVHN::Initialize(std::vector<BVHN*> *bvs, int level_)
{
    if(level_ > 100) throw std::runtime_error("BVH level is over 100");
    level = level_;
    auto count = bvs->size();
    std::cout << count << std::endl;
    if(count == 0) throw new std::runtime_error("bvs->size==0 in BVHN::Initialize");
    else if(count == 1) throw new std::runtime_error("bvs->size==1 in BVHN::Initialize");

    isLeaf = false;

    box.Reset();
    for(auto const &bv : *bvs) box.Expand(bv->box); // expand box to the size of bvs collection

    // get length, width and height of the resulting box
    double dX, dY, dZ;
    box.Dimensions(dX, dY, dZ);

    std::vector<icy::BVHN*> *left = VectorFactory.take();
    std::vector<icy::BVHN*> *right = VectorFactory.take();
    left->clear();
    right->clear();
    left->reserve(bvs->size());
    right->reserve(bvs->size());

    if (dX >= dY && dX >= dZ)
    {
        std::cout << "branch 1" << std::endl;
        double center = box.centerX();
        for(auto const &bv : *bvs) {
            if(bv->box.centerX() < center) left->push_back(bv);
            else right->push_back(bv);
        }

        // make sure that there is at least one element on each side
        if(left->size() == 0)
        {
            BVHN *selected;
            double min1 = DBL_MAX;
            for(auto const &bv : *right)
                if(min1 >= bv->box.centerX()) {
                    min1 = bv->box.centerX();
                    selected = bv;
                }
            // move "selected" from left to right
            left->push_back(selected);
            size_t size_before = right->size();
            right->erase(std::remove(right->begin(), right->end(),selected),right->end());
            if(left->size() != 1) throw std::runtime_error("left->size() != 1");
            if(right->size()!= (size_before-1)) throw std::runtime_error("right->size()!= (size_before-1)");
        } else if(right->size() == 0)
        {
            BVHN *selected;
            double min1 = DBL_MAX;
            for(auto const &bv : *left)
                if(min1 >= bv->box.centerX()) {
                    min1 = bv->box.centerX();
                    selected = bv;
                }
            // move "selected" from right to left
            right->push_back(selected);
            size_t size_before = left->size();
            left->erase(std::remove(left->begin(), left->end(),selected),left->end());
            if(right->size() != 1) throw std::runtime_error("right->size() != 1");
            if(left->size()!= (size_before-1)) throw std::runtime_error("left->size()!= (size_before-1)");
        }
    }
    else if(dY >= dX && dY >= dZ)
    {
        std::cout << "branch 2" << std::endl;
        double center = box.centerY();
        for(auto const &bv : *bvs) {
            if(bv->box.centerY() < center) left->push_back(bv);
            else right->push_back(bv);
        }

        // make sure that there is at least one element on each side
        if(left->size() == 0)
        {
            BVHN *selected;
            double min1 = DBL_MAX;
            for(auto const &bv : *right)
                if(min1 >= bv->box.centerY()) {
                    min1 = bv->box.centerY();
                    selected = bv;
                }
            // move "selected" from left to right
            left->push_back(selected);
            size_t size_before = right->size();
            right->erase(std::remove(right->begin(), right->end(),selected),right->end());
            if(left->size() != 1) throw std::runtime_error("left->size() != 1");
            if(right->size()!= (size_before-1)) throw std::runtime_error("right->size()!= (size_before-1)");
        }
        else if(right->size() == 0)
        {
            BVHN *selected;
            double min1 = DBL_MAX;
            for(auto const &bv : *left)
                if(min1 >= bv->box.centerY()) {
                    min1 = bv->box.centerY();
                    selected = bv;
                }
            // move "selected" from right to left
            right->push_back(selected);
            size_t size_before = left->size();
            left->erase(std::remove(left->begin(), left->end(),selected),left->end());
            if(right->size() != 1) throw std::runtime_error("right->size() != 1");
            if(left->size()!= (size_before-1)) throw std::runtime_error("left->size()!= (size_before-1)");
        }
    }
    else
    {
        std::cout << "branch 3" << std::endl;
        double center = box.centerZ();
        for(auto const &bv : *bvs) {
            if(bv->box.centerZ() < center) left->push_back(bv);
            else right->push_back(bv);
        }

        // make sure that there is at least one element on each side
        if(left->size() == 0)
        {
            BVHN *selected;
            double min1 = DBL_MAX;
            for(auto const &bv : *right)
                if(min1 >= bv->box.centerZ()) {
                    min1 = bv->box.centerZ();
                    selected = bv;
                }
            // move "selected" from left to right
            left->push_back(selected);
            size_t size_before = right->size();
            right->erase(std::remove(right->begin(), right->end(),selected),right->end());
            if(left->size() != 1) throw std::runtime_error("left->size() != 1");
            if(right->size()!= (size_before-1)) throw std::runtime_error("right->size()!= (size_before-1)");
        }
        else if(right->size() == 0)
        {
            BVHN *selected;
            double min1 = DBL_MAX;
            for(auto const &bv : *left)
                if(min1 >= bv->box.centerZ()) {
                    min1 = bv->box.centerZ();
                    selected = bv;
                }
            // move "selected" from right to left
            right->push_back(selected);
            size_t size_before = left->size();
            left->erase(std::remove(left->begin(), left->end(),selected),left->end());
            if(right->size() != 1) throw std::runtime_error("right->size() != 1");
            if(left->size()!= (size_before-1)) throw std::runtime_error("left->size()!= (size_before-1)");
        }
    }

    if(left->size() == 1)
    {
        child1 = left->front();
        child1->level=level+1;
        if(!child1->isLeaf) throw std::runtime_error("lone child is not leaf");
    } else if(left->size() > 1) {
        child1 = BVHNFactory.take();
        child1->Initialize(left, level+1);
    } else throw std::runtime_error("left.size < 1");
    VectorFactory.release(left);

    if(right->size() == 1)
    {
        child2 = right->front();
        child2->level=level+1;
        if(!child2->isLeaf) throw std::runtime_error("lone child is not leaf");
    } else if(right->size() > 1) {
        child2 = BVHNFactory.take();
        child2->Initialize(right, level+1);
    } else throw std::runtime_error("right.size < 1");
    VectorFactory.release(right);
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
        BVHT::broad_list.push_back(elem);
        BVHT::broad_list.push_back(b->elem);
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

