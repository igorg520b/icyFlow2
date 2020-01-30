#include "generatortool.h"

//icy::GeneratorTool::GeneratorTool() {}

void icy::GeneratorTool::GenerateLBeamSetup(BeamParams *beamParams, MeshCollection *mc)
{
    mc->Clear();

    // generate beam
    Mesh *beam = new Mesh();
    GeneratorTool::GenerateBeam(beamParams, beam);
    beam->setObjectName("beam");

    // generate indenter
    Mesh *indenter = new Mesh();
    GeneratorTool::GenerateIndenter(beamParams, indenter);
    indenter->setObjectName("indenter");

    mc->mgs.push_back(beam);
    mc->mgs.push_back(indenter);

    // align the indenter
    indenter->Translate(0,0,-indenter->zmin+beam->zmax + 1e-10);
    beam->CreateUGrid();
    indenter->CreateUGrid();
    indenter->ugridActor->GetProperty()->SetColor(indenter->colors->GetColor3d("Brown").GetData());

/*
            Translation t0 = new Translation(0, 0, 0, 0);
            Translation t1 = new Translation(0, 0, -0.6, 240);
            mgIndenter.translationCollection.Add(t0);
            mgIndenter.translationCollection.Add(t1);
*/
}

void icy::GeneratorTool::GenerateIndenter(BeamParams *beamParams, Mesh *outMesh)
{
    double CharacteristicLengthIndenter = beamParams->CharacteristicLengthIndenter;
    double l2 = beamParams->beamL2;
    double l1 = beamParams->beamL1;
    double a = beamParams->beamA;
    double c = beamParams->beamGap; // beam gap
    double d = beamParams->beamMargin; // beam margin
    double h = beamParams->beamThickness; // thickness
    double indsize = beamParams->IndenterSize; // indentor size

//    gmsh::initialize();
    gmsh::clear();
    gmsh::option::setNumber("General.Terminal", 1);
    model::add("indenter1");

    double dy = l1 + c + d;
    double dx = c+ d;
    double yoffset = dy-a/2-indsize/2;

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
    outMesh->ComputeBoundingBox();

    // get triangles
    elementTags.clear();
    nodeTags2.clear();
    model::mesh::getElementsByType(2, elementTags, nodeTags2);

    outMesh->faces.resize(elementTags.size());
    for(int i=0;i<(int)elementTags.size();i++)
    {
        Face &fc = outMesh->faces[i];
        int idx0 = nodeTagMap[nodeTags2[i*3+0]];
        int idx1 = nodeTagMap[nodeTags2[i*3+1]];
        int idx2 = nodeTagMap[nodeTags2[i*3+2]];
        fc.vrts[0]=&(outMesh->nodes[idx0]);
        fc.vrts[1]=&(outMesh->nodes[idx1]);
        fc.vrts[2]=&(outMesh->nodes[idx2]);
    }
}

void icy::GeneratorTool::GenerateBeam(BeamParams *beamParams, Mesh *outMesh)
{

    double CharacteristicLengthIndenter = beamParams->CharacteristicLengthIndenter;
    double CharacteristicLengthMax = beamParams->CharacteristicLengthMax;
    double rm = beamParams->RefinementMultiplier;
    double a = beamParams->beamA; // beamA
    double b = beamParams->beamB; // beamB
    double l1 = beamParams->beamL1;
    double l2 = beamParams->beamL2;
    double c = beamParams->beamGap; // beam gap
    double d = beamParams->beamMargin; // beam margin
    double h = beamParams->beamThickness; // thickness
    double indsize = beamParams->IndenterSize; // indentor size

    gmsh::clear();
    gmsh::option::setNumber("General.Terminal", 1);
    model::add("beam1");

    double dy = l1 + c + d;
    double dx = c+d;

    int point1 = factory::addPoint(dx + 0, dy + 0, 0, rm);
    int point2 = factory::addPoint(dx + l2, dy + 0, 0, rm);
    int point22 = factory::addPoint(dx + l2 / 2,dy + 0, 0, rm);
    int point3 = factory::addPoint(dx + l2, dy - a, 0, rm);
    int point4 = factory::addPoint(dx + 0,dy - l1 + c / 2, 0, rm);
    int point5 = factory::addPoint(dx - c,dy - l1 + c / 2, 0, rm);
    int point6 = factory::addPoint(dx - c / 2,dy - l1 + c / 2, 0, rm);
    int point7 = factory::addPoint(dx - c, dy + c, 0, 1.0);
    int point8 = factory::addPoint(dx + l2 + c, dy + c, 0, 1.0);
    int point9 = factory::addPoint(dx + l2 + c, dy - a - c, 0, 1.0);
    int point10 = factory::addPoint(dx + b + 2 * c,dy - a, 0, rm);
    int point11 = factory::addPoint(dx + b,dy - a - 2 * c, 0, rm);
    int point12 = factory::addPoint(dx + b + 2 * c,dy - a - 2 * c,0, rm);
    int point13 = factory::addPoint(dx + b + 2 * c,dy - a - c, 0, rm);
    int point14 = factory::addPoint(dx + b + c,dy - a - 2 * c, 0, rm);
    int point15 = factory::addPoint(dx + b,dy - l1 + c / 2, 0, rm);
    int point16 = factory::addPoint(dx + b + c,dy - l1 + c / 2, 0, rm);
    int point17 = factory::addPoint(dx + b + c / 2,dy - l1 + c / 2, 0, rm);
    int point18 = factory::addPoint(-d, dy + c + d, 0, 1.0);
    int point19 = factory::addPoint(dx + l2 + c + d, dy + c + d, 0, 1.0);
    int point20 = factory::addPoint(dx + l2 + c + d, dy - l1 - c - 2*d, 0, 1.0);
    int point21 = factory::addPoint(-d, -d, 0, rm);
    int point23 = factory::addPoint(dx/4 + l2/4 + c/4 + d/4, -d , 0, rm);
    int point24 = factory::addPoint(-d, dy/4 + c/4 + d/4, 0, rm);
    int point25 = factory::addPoint(-d , dy/2 + c/2 + d/2, 0, 1.0);
    int point26 = factory::addPoint(dx/2 + l2/2 + c/2 + d/2, -d, 0, 1.0);

    int circle1 = factory::addCircleArc(point4, point6, point5);
    int circle8 = factory::addCircleArc(point10, point12, point11);
    int circle11 = factory::addCircleArc(point13, point12, point14);
    int circle14 = factory::addCircleArc(point16, point17, point15);

    int line2 = factory::addLine(point1, point22);
    int line19 = factory::addLine(point22, point2);
    int line3 = factory::addLine(point2, point3);
    int line4 = factory::addLine(point4, point1);
    int line5 = factory::addLine(point5, point7);
    int line6 = factory::addLine(point7, point8);
    int line7 = factory::addLine(point8, point9);
    int line9 = factory::addLine(point3, point10);
    int line10 = factory::addLine(point9, point13);
    int line12 = factory::addLine(point11, point15);
    int line13 = factory::addLine(point14, point16);
    int line15 = factory::addLine(point18, point19);
    int line16 = factory::addLine(point19, point20);
    int line17 = factory::addLine(point20, point26);
    int line18 = factory::addLine(point26, point23);
    int line20 = factory::addLine(point23, point21);
    int line21 = factory::addLine(point21, point24);
    int line22 = factory::addLine(point24, point25);
    int line23 = factory::addLine(point25, point18);

    std::vector<int> curveTags;
    curveTags.push_back(line15);
    curveTags.push_back(line16);
    curveTags.push_back(line17);
    curveTags.push_back(line18);
    curveTags.push_back(line20);
    curveTags.push_back(line21);
    curveTags.push_back(line22);
    curveTags.push_back(line23);
    int loop2 = factory::addCurveLoop(curveTags);

    curveTags.clear();
    curveTags.push_back(line6);
    curveTags.push_back(line7);
    curveTags.push_back(line10);
    curveTags.push_back(circle11);
    curveTags.push_back(line13);
    curveTags.push_back(circle14);
    curveTags.push_back(-line12);
    curveTags.push_back(-circle8);
    curveTags.push_back(-line9);
    curveTags.push_back(-line3);
    curveTags.push_back(-line19);
    curveTags.push_back(-line2);
    curveTags.push_back(-line4);
    curveTags.push_back(circle1);
    curveTags.push_back(line5);
    int loop3 = factory::addCurveLoop(curveTags);

    std::vector<int> loops;
    loops.push_back(loop2);
    loops.push_back(loop3);
    int surfaceTag = factory::addPlaneSurface(loops);

    factory::synchronize();

    gmsh::vectorpair vp;
    gmsh::vectorpair vpOut;

    model::getEntities(vp, 2);
    factory::extrude(vp, 0, 0, h, vpOut);

    factory::synchronize();
    gmsh::option::setNumber("Mesh.CharacteristicLengthMax", CharacteristicLengthMax);
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
    outMesh->ComputeBoundingBox();

    // get triangles
    elementTags.clear();
    nodeTags2.clear();
    model::mesh::getElementsByType(2, elementTags, nodeTags2);

    outMesh->faces.resize(elementTags.size());
    for(int i=0;i<(int)elementTags.size();i++)
    {
        Face &fc = outMesh->faces[i];
        int idx0 = nodeTagMap[nodeTags2[i*3+0]];
        int idx1 = nodeTagMap[nodeTags2[i*3+1]];
        int idx2 = nodeTagMap[nodeTags2[i*3+2]];
        fc.vrts[0]=&(outMesh->nodes[idx0]);
        fc.vrts[1]=&(outMesh->nodes[idx1]);
        fc.vrts[2]=&(outMesh->nodes[idx2]);
    }
}
