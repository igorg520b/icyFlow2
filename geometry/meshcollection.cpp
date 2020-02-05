#include<algorithm>
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
    int count = 0;
    for(auto &mesh : mgs) {
        mesh->IdentifySurfaceElements();
        for(auto &elem : mesh->surfaceElements) {
            elem->globalElementId = count++;
            surfaceElements.push_back(elem);
        }
    }

    bvh.Construct(surfaceElements);
}


void icy::MeshCollection::Prepare()
{
    // populate allNodes, activeNodes; assign altId
    allNodes.clear();
    activeNodes.clear();

    allNodes.reserve(100000);
    activeNodes.reserve(100000);

    int count = 0;
    int globalCount = 0;
    for(auto const &mesh : mgs) {
        for(auto &node : mesh->nodes) {
            node.globalNodeId = globalCount++;
            allNodes.push_back(&node);
            if(!node.anchored) {
                node.altId = count++;
                activeNodes.push_back(&node);
            }
            else node.altId = -1;
        }
        mesh->ConnectFaces();
    }


    // populate elasticElements
    elasticElements.clear();
    for(auto &elem : beam->elems) elasticElements.push_back(&elem);

    // populate nonFailedCZs, allCZs, failedCZs;

    // surface elements
    /*
            // this is done after surface elements are marked
            List<Element> surfaceElementList = new List<Element>(mgs.Sum(mg => mg.surfaceElements.Count));
            foreach (Mesh mg in mgs) surfaceElementList.AddRange(mg.surfaceElements);
            surfaceElements = surfaceElementList.ToArray();
            for (int i = 0; i < surfaceElements.Length; i++) surfaceElements[i].globalElementId = i;
*/
}


void icy::MeshCollection::UpdateActors()
{
    for(auto &mesh : mgs) mesh->UpdateUGrid();
}
