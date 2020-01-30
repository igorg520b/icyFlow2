#include "meshcollection.h"

icy::MeshCollection::MeshCollection()
{

}

void icy::MeshCollection::Clear()
{
    for(auto &mesh : mgs) delete mesh;
    mgs.clear();
}

void icy::MeshCollection::IdentifySurfaceElements()
{

}

void icy::MeshCollection::ConstructBVH()
{
    // list surface elements
    surfaceElements.clear();
    for(auto &mesh : mgs) {
        mesh->IdentifySurfaceElements();
        for(auto &elem : mesh->elems) surfaceElements.push_back(&elem);
    }

    bvh.Construct(surfaceElements);
}
