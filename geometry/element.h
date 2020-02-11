#ifndef ELEMENT_H
#define ELEMENT_H

#include "node.h"
#include "face.h"
#include <unordered_set>
#include <vector>

namespace icy {
class Element;
class Node;
class Face;
}

class icy::Element
{
public:
    icy::Node* vrts[4];                 // References to element vertices
    int granule = 0;                    // Granule that the element belongs to (0 if not granular)
    int globalElementId = -1;           // Sequential numbering of exterior (surface) elements for collision detection
    int id = 0;                         // sequential numbering within same mesh for saving/loading faces
    bool isSurface = false;
    int tag = 0;                        // for cz insertion

    std::unordered_set<icy::Face*> adjFaces; // Collection of adjacent faces for collision detection (only on surface elements)
    double stress[6]={};
    double principal_stresses[3]={};

    // storage space for parallel computaiton/assembly of elastic forces and derivatives
    double rhs[12];
    double lhs[12][12]; // store left-hand side of the linearized equation of motion

    Element();
    void FindAdjFaces();
    double volume();
};

#endif // ELEMENT_H
