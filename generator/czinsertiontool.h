#ifndef CZINSERTIONTOOL_H
#define CZINSERTIONTOOL_H

#include "geometry/mesh.h"
#include <list>
#include <vector>
#include <map>
#include <stdexcept>
#include <unordered_set>

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
    std::unordered_set<int> grains;
    std::vector<int> elementsOfNode; // connected elems

    void SplitNode(std::list<int> &allNodesAsLL, std::vector<ExtendedNode> &allNodes)
    {
        // NI
    }
};


class icy::ExtendedFace {
public:
    int id = -1;
    int vrts[3] = {};
    bool exposed = true;            // is this an outside surface
    bool created = false;           // got exposed at simulation time due to fracture
    std::vector<int>elementsOfFace;
};

class icy::ExtendedElement {
public:
    int id;
    int vrts[4];
    int tag;

    static const int myConvention[4][3];

    std::vector<int> faces;
    void SubstituteNode(int oldNode, int newNode);
    int WhichFace(ExtendedFace &fc);

    bool ContainsFace(ExtendedFace &fc);
    void AddFaceIfContains(ExtendedFace &fc);
    ExtendedFace CreateFaceFromIdx(int idx);
    static int TetraIdxOfFaceVertex(int fidx, int faceVertIdx)
    { return myConvention[fidx][faceVertIdx]; }

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


    CZInsertionTool() {}
    void InsertCohesiveElements(Mesh &mg);
private:
    void Extend(Mesh &mg); // convert to the extended format stored in this class
    void ConvertBack(Mesh &mg); // clear and re-create from scratch
    void IdentifyParentsOfTriangles();
};




#endif // CZINSERTIONTOOL_H
