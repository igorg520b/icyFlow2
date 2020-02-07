#ifndef MESHCOLLECTION_H
#define MESHCOLLECTION_H

#include <vector>
#include "mesh.h"
#include "node.h"
#include "element.h"
#include "cz.h"
#include "face.h"
#include "simulation/bvh/bvht.h"
#include <vtkCellData.h>

namespace icy {
class MeshCollection;
}

class icy::MeshCollection
{
public:
    BVHT bvh;
    std::vector<icy::Mesh*> mgs;
    icy::Mesh* beam = nullptr;      // valid in the context of beam-indenter setup
    icy::Mesh* indenter = nullptr;  // valid in the context of beam-indenter setup

    std::vector<icy::Node*> allNodes, activeNodes;
    std::vector<icy::Element*> surfaceElements, elasticElements;
    std::vector<icy::CZ*> nonFailedCZs, allCZs, failedCZs;
//    std::vector<icy::Face*> allFaces;

    MeshCollection();
    void Clear(); // return the collection to empty state
    void Prepare();
    void ConstructBVH();
    void UpdateActors(); // transfer current nodal positions to vtk for rendering
};

#endif // MESHCOLLECTION_H
