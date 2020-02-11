#ifndef CZINSERTIONTOOL_H
#define CZINSERTIONTOOL_H

#include "geometry/mesh.h"
#include <list>
#include <vector>
#include <map>
#include <stdexcept>

namespace icy {
class CZInsertionTool;
class ExtendedNode;
class ExtendedFace;
class ExtendedElement;
class ExtendedCZ;
}


class icy::ExtendedNode {
public:
    double x0, y0, z0;
    int id;

    bool _belongs_to_cz = false;
    std::vector<int> grains;
    std::vector<int> elementsOfNode; // connected elems

    void SplitNode(std::list<int> &allNodesAsLL, std::vector<ExtendedNode> &allNodes)
    {
        // NI
    }
};


class icy::ExtendedFace {
public:
    int vrts[3];
    bool exposed = true;            // is this an outside surface
    bool created = false;           // got exposed at simulation time due to fracture

    int elementsOfFace[2] = {};
};

class icy::ExtendedElement {
public:
    int id;
    int vrts[4];
    int tag;

    std::vector<int> faces;

    void SubstituteNode(int oldNode, int newNode) {
        bool found = false;
        for(int i=0;i<4;i++)
            if(vrts[i]==oldNode) {
                found = true;
                vrts[i] = newNode;
                break;
            }
        if(!found) throw std::runtime_error("substituted node not found");
    }
};

class icy::ExtendedCZ {
public:
    int vrts[6];

    int e0, e1; // attached elements
    bool sameGrain = false; // czs should not be insreted between elements of the same grain

};



class icy::CZInsertionTool
{
public:
    std::vector<ExtendedNode> nodes;
    std::vector<ExtendedFace> faces;
    std::vector<ExtendedElement> elems;
    std::vector<ExtendedCZ> czs;
    std::list<int> allNodesAsLL;


    CZInsertionTool();
    void InsertCohesiveElements(Mesh &mg);
private:
    void Extend(Mesh &mg); // convert to the extended format stored in this class
    void ConvertBack(Mesh &mg); // clear and re-create from scratch
    void IdentifyParentsOfTriangles();
};




#endif // CZINSERTIONTOOL_H
