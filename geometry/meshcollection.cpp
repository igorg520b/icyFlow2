#include<algorithm>
#include "meshcollection.h"
#include <vtkCellData.h>
#include <vtkPointData.h>


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
    std::cout << "active nodes: " << activeNodes.size() << std::endl;

    // populate elasticElements
    elasticElements.clear();
    for(auto &elem : beam->elems) elasticElements.push_back(&elem);

    // populate nonFailedCZs, allCZs, failedCZs;

    // beam data
    beam->principalStresses_cells->Reset();
    beam->principalStresses_cells->SetName("principal_stresses");
    beam->principalStresses_cells->SetNumberOfComponents(3);
//    beam->principalStresses_cells->Resize(beam->elems.size());
    beam->principalStresses_cells->SetNumberOfTuples(beam->elems.size());
    beam->ugrid->GetCellData()->SetScalars(beam->principalStresses_cells);

    beam->verticalDisplacements_nodes->Reset();
    beam->verticalDisplacements_nodes->SetName("vertical_displacements");
    beam->verticalDisplacements_nodes->SetNumberOfComponents(1);
    beam->verticalDisplacements_nodes->Resize(beam->nodes.size());
    beam->verticalDisplacements_nodes->SetNumberOfValues(beam->nodes.size());
    beam->ugrid->GetPointData()->SetScalars(beam->verticalDisplacements_nodes);

}


void icy::MeshCollection::UpdateActors()
{
    for(auto &mesh : mgs) mesh->UpdateUGrid();
    beam->UpdateGridData();
}
