#include "meshcollection.h"

icy::MeshCollection::MeshCollection()
{

}

void icy::MeshCollection::Clear()
{
    for(auto &mesh : mgs) delete mesh;
    mgs.clear();

    allNodes.clear();
    activeNodes.clear();

    nonFailedCZs.clear();
    allCZs.clear();
    failedCZs.clear();

    surfaceElements.clear();
    elasticElements.clear();
    beam=indenter=nullptr;
}

void icy::MeshCollection::ConstructBVH()
{
    // list surface elements
    surfaceElements.clear();
    for(auto &mesh : mgs) {
        mesh->IdentifySurfaceElements();
        for(auto &elem : mesh->surfaceElements) surfaceElements.push_back(elem);
    }

    bvh.Construct(surfaceElements);
}


void icy::MeshCollection::Prepare()
{
// populate allNodes, activeNodes;

    // populate elasticElements

    // populate nonFailedCZs, allCZs, failedCZs;
}
