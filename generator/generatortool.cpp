#include "generatortool.h"

//icy::GeneratorTool::GeneratorTool() {}

void icy::GeneratorTool::GenerateLBeamSetup(BeamParams *beamParams, MeshCollection *mc)
{
    mc->Clear();

    // generate beam
    Mesh *beam = new Mesh();

    // generate indenter
    Mesh *indenter = new Mesh();

    mc->mgs.push_back(beam);
    mc->mgs.push_back(indenter);

    // align the indenter

}

void icy::GeneratorTool::GenerateIndenter(BeamParams *beamParams, Mesh *outMesh)
{
    double CharacteristicLengthIndenter = beamParams->CharacteristicLengthIndenter;
//    double a = beamParams->beamA; // beamA
//    double b = beamParams->beamB; // beamB
//    double l1 = beamParams->beamL1;
    double l2 = beamParams->beamL2;
    double c = beamParams->beamGap; // beam gap
    double d = beamParams->beamMargin; // beam margin
    double h = beamParams->beamThickness; // thickness
    double indsize = beamParams->IndenterSize; // indentor size

    gmsh::initialize();
    gmsh::option::setNumber("General.Terminal", 1);
    model::add("indenter1");

//    double dy = l1 + c + d;
    double dx = d;
    double yoffset = -indsize/2;

    int point3 = factory::addPoint(dx + l2, yoffset, 0+2*h, 1.0);
    int point4 = factory::addPoint(dx + l2-indsize, yoffset, 0+2*h, 1.0);
    int point5 = factory::addPoint(dx + l2-indsize/2, yoffset, h+2*h, 1.0);
    int point6 = factory::addPoint(dx + l2-indsize, yoffset, c/2+2*h, 1.0);
    int point7 = factory::addPoint(dx + l2, yoffset, c/2+2*h, 1.0);

    int circle1 = factory::addCircleArc(point3, point5, point4);
    int line2 = factory::addLine(point4, point6);
    int line3 = factory::addLine(point6, point7);
    int line4 = factory::addLine(point7, point3);

    std::vector<int> curveTags;
    curveTags.push_back(circle1);
    curveTags.push_back(line2);
    curveTags.push_back(line3);
    curveTags.push_back(line4);
    int loopTag = factory::addCurveLoop(curveTags);

    std::vector<int> loops;
    loops.push_back(loopTag);
    factory::addPlaneSurface(loops);

    factory::synchronize();

    gmsh::vectorpair vp;
    gmsh::vectorpair vpOut;

    model::getEntities(vp, 2);
    factory::extrude(vp, 0, indsize, 0, vpOut);

    factory::synchronize();
    gmsh::option::setNumber("Mesh.CharacteristicLengthMax", CharacteristicLengthIndenter);
    model::mesh::generate(3);

    // retrieve nodes from gmsh and put them into outMesh one-by-one
    std::vector<std::size_t> nodeTags1;
    std::vector<double> nodeCoords, parametricCoords;
    model::mesh::getNodesByElementType(4, nodeTags1, nodeCoords, parametricCoords, -1, false);

    // reorder node tags
    std::map<int,int> nodeTagMap; // "gmsh tag" -> "sequential tag"

    outMesh->nodes.resize(nodeTags1.size());
    for(int i = 0; i<(int)nodeTags1.size(); i++) {
        outMesh->nodes[i].Initialize(nodeCoords[i*3+0], nodeCoords[i*3+1], nodeCoords[i*3+2], i);
        nodeTagMap[(int)nodeTags1[i]] = i;
    }

    // retrieve elements
    std::vector<std::size_t> elementTags, nodeTags2;
    model::mesh::getElementsByType(4, elementTags, nodeTags2);

    outMesh->elems.resize(elementTags.size());

    for(int i=0;i<(int)elementTags.size();i++)
    {
        Element *elem = &(outMesh->elems[i]);
        int idx0 = nodeTagMap[nodeTags2[i*4+0]];
        int idx1 = nodeTagMap[nodeTags2[i*4+1]];
        int idx2 = nodeTagMap[nodeTags2[i*4+2]];
        int idx3 = nodeTagMap[nodeTags2[i*4+3]];
        elem->vrts[0] = &(outMesh->nodes[idx0]);
        elem->vrts[1] = &(outMesh->nodes[idx1]);
        elem->vrts[2] = &(outMesh->nodes[idx2]);
        elem->vrts[3] = &(outMesh->nodes[idx3]);
    }




/*
//    vtkSmartPointer<vtkUnstructuredGrid> ugrid =
//            vtkSmartPointer<vtkUnstructuredGrid>::New();
    ugrid->Reset();

    ugrid->SetPoints(points);

    ugrid->Allocate(elementTags.size());

    vtkIdType pts2[4];
    for(int i=0;i<elementTags.size();i++)
    {
        pts2[0] = map[nodeTags2[i*4+0]];
        pts2[1] = map[nodeTags2[i*4+1]];
        pts2[2] = map[nodeTags2[i*4+2]];
        pts2[3] = map[nodeTags2[i*4+3]];
        ugrid->InsertNextCell(VTK_TETRA, 4,pts2);
    }




    vtkSmartPointer<vtkDataSetMapper> ugridMapper =
            vtkSmartPointer<vtkDataSetMapper>::New();
    ugridMapper->SetInputData(ugrid);

    vtkSmartPointer<vtkActor> ugridActor =
            vtkSmartPointer<vtkActor>::New();
    ugridActor->SetMapper(ugridMapper);
    ugridActor->GetProperty()->SetColor(colors->GetColor3d("Seashell").GetData());
    ugridActor->GetProperty()->EdgeVisibilityOn();


    renderer->AddActor(ugridActor);

    renderer->ResetCamera();
    renderWindow->Render();

    vtkActorCollection *vac = renderer->GetActors();
    qDebug() << "number of actors: " << vac->GetNumberOfItems();
*/
}

void icy::GeneratorTool::GenerateBeam(BeamParams *beamParams, Mesh *output)
{

}
