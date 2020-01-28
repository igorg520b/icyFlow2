#include "meshcollection.h"

icy::MeshCollection::MeshCollection()
{

}

void icy::MeshCollection::Clear()
{
    for(auto &mesh : mgs) delete mesh;
    mgs.clear();
}
