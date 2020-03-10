#include "generatortool.h"
#include "czinsertiontool.h"
#include <unordered_set>
#include <utility>
#include <cmath>
#include <cfloat>

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
    mc->beam = beam;
    mc->indenter = indenter;
    for(auto &nd : indenter->nodes) nd.anchored = true;

    // align the indenter
    indenter->Translate(0,0,-indenter->zmin+beam->zmax + 1e-10);
    indenter->CreateUGrid();
    indenter->ugridActor->GetProperty()->SetColor(indenter->colors->GetColor3d("Brown").GetData());

    // insert czs
    CZInsertionTool cztool;
    cztool.InsertCohesiveElements(*beam);
    beam->AnchorSides();
    beam->CreateUGrid();
}

void icy::GeneratorTool::GenerateIndenter(BeamParams *beamParams, Mesh *outMesh)
{
    outMesh->isIndenter = true;
    double CharacteristicLengthIndenter = beamParams->CharacteristicLengthIndenter;
    double l2 = beamParams->beamL2;
    double l1 = beamParams->beamL1;
    double a = beamParams->beamA;
    double c = beamParams->beamGap; // beam gap
    double d = beamParams->beamMargin; // beam margin
    double h = beamParams->beamThickness; // thickness
    double indsize = beamParams->IndenterSize; // indentor size

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

    // process the result

    std::vector<std::size_t> nodeTags3;
    std::vector<double> nodeCoords3, parametricCoords3;
    model::mesh::getNodes(nodeTags3, nodeCoords3, parametricCoords3, -1, -1, false, false);

    // retrieve elements
    std::vector<std::size_t> elementTags, nodeTagsInElems;
    model::mesh::getElementsByType(4, elementTags, nodeTagsInElems);
    outMesh->elems.resize(elementTags.size());

    // get triangles
    std::vector<std::size_t> faceTags, nodeTagsInFaces;
    model::mesh::getElementsByType(2, faceTags, nodeTagsInFaces);
    outMesh->faces.resize(faceTags.size());

    // compile a set of node tags that are present in elements
    std::unordered_set<int> usedNodes;
    for(int i=0;i<(int)nodeTagsInElems.size();i++) usedNodes.insert(nodeTagsInElems[i]);
    for(int i=0;i<(int)nodeTagsInFaces.size();i++) usedNodes.insert(nodeTagsInFaces[i]);

    std::map<int,int> nodeTagMap; // "gmsh tag" -> "sequential tag"
    int count = 0;
    outMesh->nodes.resize(usedNodes.size());
    for(int i=0;i<(int)nodeTags3.size();i++) {
        int tag = nodeTags3[i];
        if(usedNodes.find(tag)!=usedNodes.end() && nodeTagMap.find(tag)==nodeTagMap.end())
        {
            if(count >= (int)usedNodes.size()) throw std::runtime_error("generator tool error");
            outMesh->nodes[count].Initialize(nodeCoords3[i*3+0], nodeCoords3[i*3+1], nodeCoords3[i*3+2], count);
            nodeTagMap[tag]=count;
            count++;
        }
    }

    // elements & faces
    for(int i=0;i<(int)elementTags.size();i++)
    {
        Element *elem = &(outMesh->elems[i]);
        for(int j=0;j<4;j++) {
            int ndidx = nodeTagsInElems[i*4+j];
            int newIdx = nodeTagMap[ndidx];
            elem->vrts[j] = &(outMesh->nodes[newIdx]);
        }
    }

    for(int i=0;i<(int)faceTags.size();i++)
    {
        Face &fc = outMesh->faces[i];
        for(int j=0;j<3;j++) {
            int ndidx = nodeTagsInFaces[i*3+j];
            int newIdx = nodeTagMap[ndidx];
            fc.vrts[j]=&(outMesh->nodes[newIdx]);
        }
    }
    outMesh->ComputeBoundingBox();
}

void icy::GeneratorTool::GenerateBeam(BeamParams *beamParams, Mesh *outMesh)
{
    outMesh->isDeformable = true;
    double CharacteristicLengthMax = beamParams->CharacteristicLengthMax;
    double rm = beamParams->RefinementMultiplier;
    double a = beamParams->beamA; // beamA
    double b = beamParams->beamB; // beamB
    double l1 = beamParams->beamL1;
    double l2 = beamParams->beamL2;
    double c = beamParams->beamGap; // beam gap
    double d = beamParams->beamMargin; // beam margin
    double h = beamParams->beamThickness; // thickness

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

    // process the result

    std::vector<std::size_t> nodeTags3;
    std::vector<double> nodeCoords3, parametricCoords3;
    model::mesh::getNodes(nodeTags3, nodeCoords3, parametricCoords3, -1, -1, false, false);

    // retrieve elements
    std::vector<std::size_t> elementTags, nodeTagsInElems;
    model::mesh::getElementsByType(4, elementTags, nodeTagsInElems);
    outMesh->elems.resize(elementTags.size());

    // get triangles
    std::vector<std::size_t> faceTags, nodeTagsInFaces;
    model::mesh::getElementsByType(2, faceTags, nodeTagsInFaces);
    outMesh->faces.resize(faceTags.size());

    // compile a set of node tags that are present in elements
    std::unordered_set<int> usedNodes;
    for(int i=0;i<(int)nodeTagsInElems.size();i++) usedNodes.insert(nodeTagsInElems[i]);
    for(int i=0;i<(int)nodeTagsInFaces.size();i++) usedNodes.insert(nodeTagsInFaces[i]);

    std::map<int,int> nodeTagMap; // "gmsh tag" -> "sequential tag"
    int count = 0;
    outMesh->nodes.resize(usedNodes.size());
    for(int i=0;i<(int)nodeTags3.size();i++) {
        int tag = nodeTags3[i];
        if(usedNodes.find(tag)!=usedNodes.end() && nodeTagMap.find(tag)==nodeTagMap.end())
        {
            if(count >= (int)usedNodes.size()) throw std::runtime_error("generator tool error");
            outMesh->nodes[count].Initialize(nodeCoords3[i*3+0], nodeCoords3[i*3+1], nodeCoords3[i*3+2], count);
            nodeTagMap[tag]=count;
            count++;
        }
    }

    // elements & faces
    for(int i=0;i<(int)elementTags.size();i++)
    {
        Element *elem = &(outMesh->elems[i]);
        for(int j=0;j<4;j++) {
            int ndidx = nodeTagsInElems[i*4+j];
            int newIdx = nodeTagMap[ndidx];
            elem->vrts[j] = &(outMesh->nodes[newIdx]);
        }
    }

    for(int i=0;i<(int)faceTags.size();i++)
    {
        Face &fc = outMesh->faces[i];
        for(int j=0;j<3;j++) {
            int ndidx = nodeTagsInFaces[i*3+j];
            int newIdx = nodeTagMap[ndidx];
            fc.vrts[j]=&(outMesh->nodes[newIdx]);
        }
    }
    outMesh->ComputeBoundingBox();

    // region of cz insertion
    std::vector<p2d> region;
    region.push_back(std::make_pair(dx + l2*0.6, dy+c/2));
    region.push_back(std::make_pair(dx + b + 2*c, dy -a-c/2));
    region.push_back(std::make_pair(dx + b + c/2, dy - a - 2*c));
    region.push_back(std::make_pair(dx + b + c/2, dy - l1 - c - d/2));
    region.push_back(std::make_pair(dx - c/2, dy - l1 - c - d/2));
    region.push_back(std::make_pair(dx - c/2, dy - l1 / 2));
    region.push_back(std::make_pair(dx - c / 2, dy + c / 2));


    // exterior point for the region
    double maxX = -DBL_MAX;
    double maxY = -DBL_MAX;
    for(auto &pt : region) {
        if(pt.first > maxX) maxX=pt.first;
        if(pt.second > maxY) maxY=pt.second;
    }
    p2d exteriorPt = std::make_pair(maxX+1, maxY+1);

    count = 2;
    for(auto &elem : outMesh->elems) {
        double centerX = 0, centerY = 0;
        for(int i=0;i<4;i++) {
            centerX+=elem.vrts[i]->x0/4;
            centerY+=elem.vrts[i]->y0/4;
        }
        p2d elemCenter = std::make_pair(centerX, centerY);
        if(PointInsideLoop(elemCenter, region, exteriorPt)) elem.tag = count++;
        else elem.tag = 1;
    }
}

bool icy::GeneratorTool::PointInsideLoop(
        p2d p, std::vector<p2d> &loop, p2d exteriorPt) {

    int count = 0;
    for (int i = 0; i < (int)loop.size(); i++)
    {
        p2d n0 = loop[i];
        p2d n1 = loop[(i + 1) % (int)loop.size()];
        if (SegmentSegmentIntersection(n0, n1, p, exteriorPt)) count++;
    }
    return (count % 2) == 1; // true if interior
}

bool icy::GeneratorTool::SegmentSegmentIntersection(p2d p, p2d p2, p2d q, p2d q2)
{
    auto cross = [](p2d v1, p2d v2) {return v1.first*v2.second-v1.second*v2.first;};
    auto sub = [](p2d v1, p2d v2) {return std::make_pair(v1.first-v2.first,v1.second-v2.second);};

    p2d r = sub(p2,p);
    p2d s = sub(q2,q);

    double rxs = cross(r, s);

    if (abs(rxs) < 1E-30) return false; // collinear (disjoint in our case)

    double p_zeta = cross(sub(q,p), s) / rxs;
    double q_zeta = cross(sub(q,p), r) / rxs;

    if (p_zeta > 0 && p_zeta < 1 && q_zeta > 0 && q_zeta < 1) return true;
    else return false;
}


//==================== generate beam for 4-point bending

void icy::GeneratorTool::GenerateCantileverBeamSetup(BeamParams *beamParams, MeshCollection *mc)
{
    mc->Clear();

    // generate beam
    Mesh *beam = new Mesh();
    GeneratorTool::GenerateCantileverBeam(beamParams, beam);
    beam->setObjectName("beam");
    mc->mgs.push_back(beam);
    mc->beam = beam;
    CZInsertionTool cztool;
    cztool.InsertCohesiveElements(*beam);
    beam->CreateUGrid();


    // generate indenters
    Mesh *indenter = new Mesh();
    GeneratorTool::GenerateCantileverIndenter(beamParams, indenter);
    indenter->setObjectName("indenter");

    mc->mgs.push_back(indenter);
    mc->indenter = indenter;
    for(auto &nd : indenter->nodes) nd.anchored = true;

    // align the indenter
    indenter->Translate(0,0,-indenter->zmin+beam->zmax + 1e-10);
    indenter->CreateUGrid();
    indenter->ugridActor->GetProperty()->SetColor(indenter->colors->GetColor3d("Brown").GetData());

    //support
    Mesh *support = new Mesh();
    GeneratorTool::GenerateCantileverSupport(beamParams, support);
    indenter->setObjectName("support");

    mc->mgs.push_back(support);
    for(auto &nd : support->nodes) nd.anchored = true;

    // align the support
    support->Translate(0,0,-indenter->zmin+beam->zmax + 1e-10);
    support->CreateUGrid();
    support->ugridActor->GetProperty()->SetColor(indenter->colors->GetColor3d("Blue").GetData());
}


void icy::GeneratorTool::GenerateCantileverBeam(BeamParams *beamParams, Mesh *outMesh)
{
    outMesh->isDeformable = true;
    double CharacteristicLengthMax = beamParams->CharacteristicLengthMax;
    double a = beamParams->beamA; // beamA
    double l1 = beamParams->beamL1;
    double c = beamParams->beamGap; // beam gap
    double h = beamParams->beamThickness; // thickness

    gmsh::clear();
    gmsh::option::setNumber("General.Terminal", 1);
    model::add("beam1");


    factory::addBox(0, 0, 0, l1, a, h);
    factory::synchronize();

    factory::synchronize();
    gmsh::option::setNumber("Mesh.CharacteristicLengthMax", CharacteristicLengthMax);
    model::mesh::generate(3);

    // process the result

    std::vector<std::size_t> nodeTags3;
    std::vector<double> nodeCoords3, parametricCoords3;
    model::mesh::getNodes(nodeTags3, nodeCoords3, parametricCoords3, -1, -1, false, false);

    // retrieve elements
    std::vector<std::size_t> elementTags, nodeTagsInElems;
    model::mesh::getElementsByType(4, elementTags, nodeTagsInElems);
    outMesh->elems.resize(elementTags.size());

    // get triangles
    std::vector<std::size_t> faceTags, nodeTagsInFaces;
    model::mesh::getElementsByType(2, faceTags, nodeTagsInFaces);
    outMesh->faces.resize(faceTags.size());

    // compile a set of node tags that are present in elements
    std::unordered_set<int> usedNodes;
    for(int i=0;i<(int)nodeTagsInElems.size();i++) usedNodes.insert(nodeTagsInElems[i]);
    for(int i=0;i<(int)nodeTagsInFaces.size();i++) usedNodes.insert(nodeTagsInFaces[i]);

    std::map<int,int> nodeTagMap; // "gmsh tag" -> "sequential tag"
    int count = 0;
    outMesh->nodes.resize(usedNodes.size());
    for(int i=0;i<(int)nodeTags3.size();i++) {
        int tag = nodeTags3[i];
        if(usedNodes.find(tag)!=usedNodes.end() && nodeTagMap.find(tag)==nodeTagMap.end())
        {
            if(count >= (int)usedNodes.size()) throw std::runtime_error("generator tool error");
            outMesh->nodes[count].Initialize(nodeCoords3[i*3+0], nodeCoords3[i*3+1], nodeCoords3[i*3+2], count);
            nodeTagMap[tag]=count;
            count++;
        }
    }

    // elements & faces
    for(int i=0;i<(int)elementTags.size();i++)
    {
        Element *elem = &(outMesh->elems[i]);
        for(int j=0;j<4;j++) {
            int ndidx = nodeTagsInElems[i*4+j];
            int newIdx = nodeTagMap[ndidx];
            elem->vrts[j] = &(outMesh->nodes[newIdx]);
        }
    }

    for(int i=0;i<(int)faceTags.size();i++)
    {
        Face &fc = outMesh->faces[i];
        for(int j=0;j<3;j++) {
            int ndidx = nodeTagsInFaces[i*3+j];
            int newIdx = nodeTagMap[ndidx];
            fc.vrts[j]=&(outMesh->nodes[newIdx]);
        }
    }
    outMesh->ComputeBoundingBox();

    count = 1;
    for(auto &elem : outMesh->elems) elem.tag = count++;
}


void icy::GeneratorTool::GenerateCantileverIndenter(BeamParams *beamParams, Mesh *outMesh)
{
    outMesh->isIndenter = true;
    double CharacteristicLengthIndenter = beamParams->CharacteristicLengthIndenter;
    double l1 = beamParams->beamL1;
    double a = beamParams->beamA;
    double d = beamParams->beamMargin; // beam margin
    double h = beamParams->beamThickness; // thickness
    double indsize = beamParams->IndenterSize; // indentor size
    double c = beamParams->beamGap; // beam gap


    gmsh::clear();
    gmsh::option::setNumber("General.Terminal", 1);
    model::add("indenter1");

    double dx1 = indsize/2 + (l1-indsize)/3;
    double dx2 = l1-dx1;

    int point3 = factory::addPoint(dx1+indsize/2, 0, 0, 1.0);
    int point4 = factory::addPoint(dx1-indsize/2, 0, 0, 1.0);
    int point5 = factory::addPoint(dx1, 0, h, 1.0);
    int point6 = factory::addPoint(dx1-indsize/2, 0, c/2, 1.0);
    int point7 = factory::addPoint(dx1+indsize/2, 0, c/2, 1.0);

    int point13 = factory::addPoint(dx2+indsize/2, 0, 0, 1.0);
    int point14 = factory::addPoint(dx2-indsize/2, 0, 0, 1.0);
    int point15 = factory::addPoint(dx2, 0, h, 1.0);
    int point16 = factory::addPoint(dx2-indsize/2, 0, c/2, 1.0);
    int point17 = factory::addPoint(dx2+indsize/2, 0, c/2, 1.0);

    int circle1 = factory::addCircleArc(point3, point5, point4);
    int line2 = factory::addLine(point4, point6);
    int line3 = factory::addLine(point6, point7);
    int line4 = factory::addLine(point7, point3);

    int circle11 = factory::addCircleArc(point13, point15, point14);
    int line12 = factory::addLine(point14, point16);
    int line13 = factory::addLine(point16, point17);
    int line14 = factory::addLine(point17, point13);

    std::vector<int> curveTags;
    curveTags.push_back(circle1);
    curveTags.push_back(line2);
    curveTags.push_back(line3);
    curveTags.push_back(line4);
    int loopTag1 = factory::addCurveLoop(curveTags);

    std::vector<int> curveTags2;
    curveTags2.push_back(circle11);
    curveTags2.push_back(line12);
    curveTags2.push_back(line13);
    curveTags2.push_back(line14);
    int loopTag2 = factory::addCurveLoop(curveTags2);

    std::vector<int> loops1;
    loops1.push_back(loopTag1);
    std::vector<int> loops2;
    loops2.push_back(loopTag2);
    factory::addPlaneSurface(loops1);
    factory::addPlaneSurface(loops2);

    factory::synchronize();

    gmsh::vectorpair vp;
    gmsh::vectorpair vpOut;

    model::getEntities(vp, 2);
    factory::extrude(vp, 0, a, 0, vpOut);

    factory::synchronize();
    gmsh::option::setNumber("Mesh.CharacteristicLengthMax", CharacteristicLengthIndenter);
    model::mesh::generate(3);

    // process the result

    std::vector<std::size_t> nodeTags3;
    std::vector<double> nodeCoords3, parametricCoords3;
    model::mesh::getNodes(nodeTags3, nodeCoords3, parametricCoords3, -1, -1, false, false);

    // retrieve elements
    std::vector<std::size_t> elementTags, nodeTagsInElems;
    model::mesh::getElementsByType(4, elementTags, nodeTagsInElems);
    outMesh->elems.resize(elementTags.size());

    // get triangles
    std::vector<std::size_t> faceTags, nodeTagsInFaces;
    model::mesh::getElementsByType(2, faceTags, nodeTagsInFaces);
    outMesh->faces.resize(faceTags.size());

    // compile a set of node tags that are present in elements
    std::unordered_set<int> usedNodes;
    for(int i=0;i<(int)nodeTagsInElems.size();i++) usedNodes.insert(nodeTagsInElems[i]);
    for(int i=0;i<(int)nodeTagsInFaces.size();i++) usedNodes.insert(nodeTagsInFaces[i]);

    std::map<int,int> nodeTagMap; // "gmsh tag" -> "sequential tag"
    int count = 0;
    outMesh->nodes.resize(usedNodes.size());
    for(int i=0;i<(int)nodeTags3.size();i++) {
        int tag = nodeTags3[i];
        if(usedNodes.find(tag)!=usedNodes.end() && nodeTagMap.find(tag)==nodeTagMap.end())
        {
            if(count >= (int)usedNodes.size()) throw std::runtime_error("generator tool error");
            outMesh->nodes[count].Initialize(nodeCoords3[i*3+0], nodeCoords3[i*3+1], nodeCoords3[i*3+2], count);
            nodeTagMap[tag]=count;
            count++;
        }
    }

    // elements & faces
    for(int i=0;i<(int)elementTags.size();i++)
    {
        Element *elem = &(outMesh->elems[i]);
        for(int j=0;j<4;j++) {
            int ndidx = nodeTagsInElems[i*4+j];
            int newIdx = nodeTagMap[ndidx];
            elem->vrts[j] = &(outMesh->nodes[newIdx]);
        }
    }

    for(int i=0;i<(int)faceTags.size();i++)
    {
        Face &fc = outMesh->faces[i];
        for(int j=0;j<3;j++) {
            int ndidx = nodeTagsInFaces[i*3+j];
            int newIdx = nodeTagMap[ndidx];
            fc.vrts[j]=&(outMesh->nodes[newIdx]);
        }
    }
    outMesh->ComputeBoundingBox();
}



void icy::GeneratorTool::GenerateCantileverSupport(BeamParams *beamParams, Mesh *outMesh)
{
    outMesh->isIndenter = true;
    double CharacteristicLengthIndenter = beamParams->CharacteristicLengthIndenter;
    double l1 = beamParams->beamL1;
    double a = beamParams->beamA;
    double d = beamParams->beamMargin; // beam margin
    double h = beamParams->beamThickness; // thickness
    double indsize = beamParams->IndenterSize; // indentor size
    double c = beamParams->beamGap; // beam gap


    gmsh::clear();
    gmsh::option::setNumber("General.Terminal", 1);
    model::add("indenter1");

    double dx1 = indsize/2;
    double dx2 = l1-dx1;

    int point3 = factory::addPoint(dx1+indsize/2, 0, 0, 1.0);
    int point4 = factory::addPoint(dx1-indsize/2, 0, 0, 1.0);
    int point5 = factory::addPoint(dx1, 0, -h, 1.0);
    int point6 = factory::addPoint(dx1-indsize/2, 0, -c/2, 1.0);
    int point7 = factory::addPoint(dx1+indsize/2, 0, -c/2, 1.0);

    int point13 = factory::addPoint(dx2+indsize/2, 0, 0, 1.0);
    int point14 = factory::addPoint(dx2-indsize/2, 0, 0, 1.0);
    int point15 = factory::addPoint(dx2, 0, -h, 1.0);
    int point16 = factory::addPoint(dx2-indsize/2, 0, -c/2, 1.0);
    int point17 = factory::addPoint(dx2+indsize/2, 0, -c/2, 1.0);

    int circle1 = factory::addCircleArc(point3, point5, point4);
    int line2 = factory::addLine(point4, point6);
    int line3 = factory::addLine(point6, point7);
    int line4 = factory::addLine(point7, point3);

    int circle11 = factory::addCircleArc(point13, point15, point14);
    int line12 = factory::addLine(point14, point16);
    int line13 = factory::addLine(point16, point17);
    int line14 = factory::addLine(point17, point13);

    std::vector<int> curveTags;
    curveTags.push_back(circle1);
    curveTags.push_back(line2);
    curveTags.push_back(line3);
    curveTags.push_back(line4);
    int loopTag1 = factory::addCurveLoop(curveTags);

    std::vector<int> curveTags2;
    curveTags2.push_back(circle11);
    curveTags2.push_back(line12);
    curveTags2.push_back(line13);
    curveTags2.push_back(line14);
    int loopTag2 = factory::addCurveLoop(curveTags2);

    std::vector<int> loops1;
    loops1.push_back(loopTag1);
    std::vector<int> loops2;
    loops2.push_back(loopTag2);
    factory::addPlaneSurface(loops1);
    factory::addPlaneSurface(loops2);

    factory::synchronize();

    gmsh::vectorpair vp;
    gmsh::vectorpair vpOut;

    model::getEntities(vp, 2);
    factory::extrude(vp, 0, a, 0, vpOut);

    factory::synchronize();
    gmsh::option::setNumber("Mesh.CharacteristicLengthMax", CharacteristicLengthIndenter);
    model::mesh::generate(3);

    // process the result

    std::vector<std::size_t> nodeTags3;
    std::vector<double> nodeCoords3, parametricCoords3;
    model::mesh::getNodes(nodeTags3, nodeCoords3, parametricCoords3, -1, -1, false, false);

    // retrieve elements
    std::vector<std::size_t> elementTags, nodeTagsInElems;
    model::mesh::getElementsByType(4, elementTags, nodeTagsInElems);
    outMesh->elems.resize(elementTags.size());

    // get triangles
    std::vector<std::size_t> faceTags, nodeTagsInFaces;
    model::mesh::getElementsByType(2, faceTags, nodeTagsInFaces);
    outMesh->faces.resize(faceTags.size());

    // compile a set of node tags that are present in elements
    std::unordered_set<int> usedNodes;
    for(int i=0;i<(int)nodeTagsInElems.size();i++) usedNodes.insert(nodeTagsInElems[i]);
    for(int i=0;i<(int)nodeTagsInFaces.size();i++) usedNodes.insert(nodeTagsInFaces[i]);

    std::map<int,int> nodeTagMap; // "gmsh tag" -> "sequential tag"
    int count = 0;
    outMesh->nodes.resize(usedNodes.size());
    for(int i=0;i<(int)nodeTags3.size();i++) {
        int tag = nodeTags3[i];
        if(usedNodes.find(tag)!=usedNodes.end() && nodeTagMap.find(tag)==nodeTagMap.end())
        {
            if(count >= (int)usedNodes.size()) throw std::runtime_error("generator tool error");
            outMesh->nodes[count].Initialize(nodeCoords3[i*3+0], nodeCoords3[i*3+1], nodeCoords3[i*3+2], count);
            nodeTagMap[tag]=count;
            count++;
        }
    }

    // elements & faces
    for(int i=0;i<(int)elementTags.size();i++)
    {
        Element *elem = &(outMesh->elems[i]);
        for(int j=0;j<4;j++) {
            int ndidx = nodeTagsInElems[i*4+j];
            int newIdx = nodeTagMap[ndidx];
            elem->vrts[j] = &(outMesh->nodes[newIdx]);
        }
    }

    for(int i=0;i<(int)faceTags.size();i++)
    {
        Face &fc = outMesh->faces[i];
        for(int j=0;j<3;j++) {
            int ndidx = nodeTagsInFaces[i*3+j];
            int newIdx = nodeTagMap[ndidx];
            fc.vrts[j]=&(outMesh->nodes[newIdx]);
        }
    }
    outMesh->ComputeBoundingBox();
}
