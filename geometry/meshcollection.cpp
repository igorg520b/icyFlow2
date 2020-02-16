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

    activeCZs.clear();
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

    UpdateCZs();

    beam->verticalDisplacements_nodes->Reset();
    beam->verticalDisplacements_nodes->SetName("vertical_displacements");
    beam->verticalDisplacements_nodes->SetNumberOfComponents(1);
    beam->verticalDisplacements_nodes->SetNumberOfValues(beam->nodes.size());
    beam->ugrid->GetPointData()->SetScalars(beam->verticalDisplacements_nodes);

    // beam data

    beam->firstPrincipalStress_cells->Reset();
    beam->firstPrincipalStress_cells->SetName("first_principal_stress");
    beam->firstPrincipalStress_cells->SetNumberOfComponents(1);
    beam->firstPrincipalStress_cells->SetNumberOfTuples(beam->elems.size());
    beam->ugrid->GetCellData()->SetScalars(beam->firstPrincipalStress_cells);


    beam->principalStresses_cells->Reset();
    beam->principalStresses_cells->SetName("principal_stresses");
    beam->principalStresses_cells->SetNumberOfComponents(3);
    beam->principalStresses_cells->SetNumberOfTuples(beam->elems.size());
 //   beam->ugrid->GetCellData()->SetScalars(beam->principalStresses_cells);

    beam->tags_cells->Reset();
    beam->tags_cells->SetName("element_tags");
    beam->tags_cells->SetNumberOfComponents(1);
    beam->tags_cells->SetNumberOfTuples(beam->elems.size());
//    beam->ugrid->GetCellData()->SetScalars(beam->tags_cells);
//    beam->ugrid->GetCellData()->AddArray(beam->tags_cells);

//    beam->ugrid->GetPointData()->SetActiveScalars("principal_stresses");

    beam->mapper->SetInputData(beam->ugrid);
    beam->mapper->ScalarVisibilityOn();
    beam->mapper->SetScalarModeToUseCellData();
//    beam->mapper->SetScalarModeToUsePointData();
    beam->mapper->SetColorModeToMapScalars();


//    beam->hueLut->SetTableRange (-1, 1e8);
//    beam->hueLut->SetHueRange (-1000000, 100000000);
//    beam->hueLut->SetSaturationRange (-1000000, 100000000);
//    beam->hueLut->SetValueRange (-100, 1e3);

    beam->hueLut->Build();
    beam->mapper->SetLookupTable(beam->hueLut);
    beam->ugridActor->SetMapper(beam->mapper);
}

void icy::MeshCollection::UpdateCZs() {
    // populate nonFailedCZs, allCZs, failedCZs;
    allCZs.clear();
    failedCZs.clear();
    activeCZs.clear();
    for(auto &cz : beam->czs) {
        allCZs.push_back(&cz);
        if(cz.failed) failedCZs.push_back(&cz);
        else activeCZs.push_back(&cz);
    }

}


void icy::MeshCollection::UpdateActors()
{
    for(auto &mesh : mgs) mesh->UpdateUGrid();
    beam->UpdateGridData();
}
