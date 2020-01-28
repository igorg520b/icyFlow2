#ifndef MESHCOLLECTION_H
#define MESHCOLLECTION_H

#include <vector>
#include "mesh.h"
#include "node.h"
#include "element.h"
#include "cz.h"
#include "face.h"

namespace icy {
class MeshCollection;
}

class icy::MeshCollection
{
public:
    MeshCollection();
    std::vector<icy::Mesh*> mgs;

    std::vector<icy::Node*> allNodes, activeNodes;
    std::vector<icy::Element*> surfaceElements, elasticElements;
    std::vector<icy::CZ*> nonFailedCZs, allCZs, failedCZs;
//    std::vector<icy::Face*> allFaces;
    std::vector<icy::Mesh*> deformables, nonDeformables, indenters;

    void Clear(); // return the collection to empty state

};

#endif // MESHCOLLECTION_H
