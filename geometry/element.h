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
    int granule;                        // Granule that the element belongs to (0 if not granular)
    int globalElementId = -1;           // Sequential numbering of exterior (surface) elements for collision detection
    int id;                             // sequential numbering within same mesh for saving/loading faces
    bool isSurface;

    std::unordered_set<icy::Face*> adjFaces; // Collection of adjacent faces for collision detection (only on surface elements)
    double stress[6];
    double principal_stresses[3];
    // public CPU_Linear_Tetrahedron.ElementExtension extension; // any additional data

    Element();
    void FindAdjFaces();
    double volume();
};

#endif // ELEMENT_H
