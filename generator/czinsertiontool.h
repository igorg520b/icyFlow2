#ifndef CZINSERTIONTOOL_H
#define CZINSERTIONTOOL_H

#include "geometry/mesh.h"
#include <list>
#include <vector>
#include <map>
#include <stdexcept>
#include <unordered_set>
#include <tuple>
#include <algorithm>

namespace icy {
class CZInsertionTool;
class ExtendedNode;
class ExtendedFace;
class ExtendedElement;
class ExtendedCZ;

}


namespace std
{
    template<> struct hash<tuple<int,int,int>>
    {
        std::size_t operator()(tuple<int,int,int> const& s) const noexcept
        {
            auto H3 = [](int a , int b, int c) {return .5*((.5*(a + b)*(a + b + 1) + b) + c)*((.5*(a + b)*(a + b + 1) + b) + c + 1) + c; };

            int v1 = std::get<0>(s);
            int v2 = std::get<1>(s);
            int v3 = std::get<2>(s);
            return H3(v1, v2,v3);
        }
    };
}

class icy::ExtendedNode {
public:
    double x0, y0, z0;
    bool isSurface = false;
    int id;

    bool belongs_to_cz = false;
    std::unordered_set<int> grains;
    std::vector<int> elementsOfNode; // connected elems

    void SplitNode(std::list<int> &allNodesAsLL,
                   std::vector<ExtendedNode> &allNodes,
                   std::vector<ExtendedElement> &elems);
};


class icy::ExtendedFace {
public:
    ExtendedFace(std::tuple<int,int,int> ftuple);
    ExtendedFace(){}
    int id = -1;
    int vrts[3] = {};
    std::tuple<int,int,int> ftuple;
    bool exposed = true;            // is this an outside surface
    bool created = false;           // got exposed at simulation time due to fracture
    std::vector<int>connectedElements;
};

class icy::ExtendedElement {
public:
    int id;
    int vrts[4];
    int tag;
    std::tuple<int,int,int> faces_as_tuples[4];
    std::vector<int> faces;

    static const int myConvention[4][3];

    void SubstituteNode(int oldNode, int newNode);
    int WhichFace(ExtendedFace &fc);

    bool ContainsFace(ExtendedFace &fc);
    void AddFaceIfContains(ExtendedFace &fc);
    ExtendedFace CreateFaceFromIdx(int idx);
    static int TetraIdxOfFaceVertex(int fidx, int faceVertIdx)
    { return myConvention[fidx][faceVertIdx]; }

    static std::tuple<int,int,int> make_sorted_tuple(int a, int b, int c);
    void GenerateFaces();

};

class icy::ExtendedCZ {
public:
    int vrts[6];

    int e0, e1; // attached elements
    int fidx0, fidx1;
    bool sameGrain = false; // czs should not be insreted between elements of the same grain
    int relativeOrientationIndex = -1; // used when creating CZ from Face
};



class icy::CZInsertionTool
{
public:
    std::vector<ExtendedNode> nodes;
    std::vector<ExtendedFace> faces;
    std::vector<ExtendedElement> elems;
    std::vector<ExtendedCZ> czs;
    std::list<int> allNodesAsLL;
    std::map<std::tuple<int,int,int>, ExtendedFace> face_map;


    CZInsertionTool() {}
    void InsertCohesiveElements(Mesh &mg);

private:
    void Extend(Mesh &mg); // convert to the extended format stored in this class
    void ConvertBack(Mesh &mg); // clear and re-create from scratch
    void IdentifyParentsOfTriangles();
};




#endif // CZINSERTIONTOOL_H
