#include "mesh.h"
#include <stdexcept>
#include <cfloat>

icy::Mesh::~Mesh()
{
    Clear();

    hueLut->SetTableRange (-0.1, 0.01);
//    hueLut->SetHueRange (-0.0015, 0);
//    hueLut->SetSaturationRange (-0.0015, 0);
//    hueLut->SetValueRange (-0.0015, 0);
    hueLut->Build();
}

void icy::Mesh::Clear()
{
    // delete everything
    surfaceElements.clear();
    faces.clear();
    czs.clear();
    elems.clear();
    nodes.clear();
}


void icy::Mesh::ComputeBoundingBox()
{

    auto result = std::minmax_element( nodes.begin(), nodes.end(),
                                 []( const Node &a, const Node &b )
                                 { return a.x0 < b.x0; } );
    xmin = result.first->x0;
    xmax = result.second->x0;

    result = std::minmax_element( nodes.begin(), nodes.end(),
                                 []( const Node &a, const Node &b )
                                 { return a.y0 < b.y0; } );
    ymin = result.first->y0;
    ymax = result.second->y0;

    result = std::minmax_element( nodes.begin(), nodes.end(),
                                 []( const Node &a, const Node &b )
                                 { return a.z0 < b.z0; } );
    zmin = result.first->z0;
    zmax = result.second->z0;
}


void icy::Mesh::AnchorSides()
{
    auto almostEqual = [](double val1, double val2) {return ((val1-val2) < 1e-10 && (val2-val1) < 1e-10);};

    ComputeBoundingBox();
    int count = 0;
    for(auto &nd : nodes) {
        if(almostEqual(nd.x0,xmin) || almostEqual(nd.x0, xmax) ||
                almostEqual(nd.y0, ymin) || almostEqual(nd.y0, ymax))
        { nd.anchored = true; count++;}
    }
//    std::cout << "nodes anchored: " << count << std::endl;
//    std::cout << "bounding box (" << xmin << "," << ymin << "," << zmin << ")  (" << xmax << "," << ymax << "," << zmax << ");" << std::endl;
}

void icy::Mesh::IdentifySurfaceElements()
{
    // note: faces must be populated before this method is invoked

    // reset all nodes
    for(auto &nd : nodes) nd.isSurface = false;

    // nodes that belong to faces are marked as surface nodes
    for(auto const &fc : faces)
        for(int i=0;i<3;i++)
            fc.vrts[i]->isSurface = true;

    surfaceElements.clear();
    for(auto &elem : elems)
        if(elem.vrts[0]->isSurface ||
                elem.vrts[1]->isSurface ||
                elem.vrts[2]->isSurface ||
                elem.vrts[3]->isSurface) {
            elem.isSurface = true;
            surfaceElements.push_back(&elem);
        }
        else elem.isSurface = false;
}


void icy::Mesh::ConnectFaces()
{
    for(auto &nd : this->nodes) nd.faces.clear();
    for(auto &f : faces)
        for(int i=0;i<3;i++)
            f.vrts[i]->faces.push_back(&f);

}

void icy::Mesh::CenterSample(double &dx, double &dy)
{
    ComputeBoundingBox();
    double x_center = (xmax + xmin) / 2.0;
    double y_center = (ymax + ymin) / 2.0;

    for(auto &nd : nodes) {
        nd.x0 = (nd.x0 - x_center);
        nd.y0 = (nd.y0 - y_center);
        nd.cx = nd.x0;
        nd.cy = nd.y0;
        nd.cz = nd.z0;
    }
    ComputeBoundingBox();
    dx = x_center;
    dy = y_center;
}

void icy::Mesh::Translate(double dx, double dy, double dz)
{
    for(auto &nd : nodes)
    {
        nd.x0 += dx;
        nd.cx = nd.x0;

        nd.y0 += dy;
        nd.cy = nd.y0;

        nd.z0 += dz;
        nd.cz = nd.z0;
    }
    ComputeBoundingBox();
}

void icy::Mesh::CreateUGrid()
{
    points->Reset();
    points->Allocate(nodes.size());

    for(auto const &nd : nodes) points->InsertPoint(nd.id,nd.x0,nd.y0,nd.z0);
    ugrid->Reset();
    ugrid->Allocate(elems.size());
    ugrid->SetPoints(points);

    vtkIdType pts2[4];
    for(auto const &elem : elems)
    {
        pts2[0] = elem.vrts[0]->id;
        pts2[1] = elem.vrts[1]->id;
        pts2[2] = elem.vrts[2]->id;
        pts2[3] = elem.vrts[3]->id;
        ugrid->InsertNextCell(VTK_TETRA, 4,pts2);
    }
    /*
     // testing - draw some faces
    for(int i=0;i<(int)faces.size();i+=2)
    {
        Face &fc = faces[i];
        pts2[0] = fc.vrts[0]->id;
        pts2[1] = fc.vrts[1]->id;
        pts2[2] = fc.vrts[2]->id;
        ugrid->InsertNextCell(VTK_TRIANGLE, 3,pts2);
    }
    */

    dataSetMapper->SetInputData(ugrid);
    ugridActor->SetMapper(dataSetMapper);
    ugridActor->GetProperty()->SetColor(colors->GetColor3d("Seashell").GetData());
    ugridActor->GetProperty()->EdgeVisibilityOn();

    if(czs.size() == 0) return;

    ugrid_czs->Reset();
    ugrid_czs->Allocate(czs.size());
    ugrid_czs->SetPoints(points);

    for(auto &cz : czs) {
        for(int i=0;i<3;i++) pts2[i]=cz.vrts[i]->id;
        ugrid_czs->InsertNextCell(VTK_TRIANGLE, 3, pts2);
    }
    dataSetMapper_czs->SetInputData(ugrid_czs);
    ugridActor_czs->SetMapper(dataSetMapper_czs);
    ugridActor_czs->GetProperty()->SetColor(colors->GetColor3d("Blue").GetData());
    ugridActor_czs->GetProperty()->EdgeVisibilityOn();
}

void icy::Mesh::UpdateUGrid()
{
    for(auto const &nd : nodes) {
        points->SetPoint(nd.id,nd.cx,nd.cy,nd.cz);
    }
    points->Modified();
}

void icy::Mesh::UpdateGridData()
{
    // transfer computed variables onto the grid for vtk display and analysis in paraview
    double min_stress = DBL_MAX;
    double max_stress = -DBL_MAX;
    for(int i=0;i<(int)nodes.size();i++) {
        double val = nodes[i].uz;
        verticalDisplacements_nodes->SetValue(i, val);
        if(min_stress>val) min_stress = val;
        if(max_stress<val) max_stress = val;
    }
    verticalDisplacements_nodes->Modified();

    for(int i=0;i<(int)elems.size();i++) {
        Element *elem = &elems[i];
        principalStresses_cells->SetTuple3(i, elem->principal_stresses[0], elem->principal_stresses[1], elem->principal_stresses[2]);
        tags_cells->SetValue(i, elem->tag);
        firstPrincipalStress_cells->SetValue(i, elem->principal_stresses[0]);

//        double ps = elem->principal_stresses[0];
//        if(min_stress>ps) min_stress = ps;
//        if(max_stress<ps) max_stress = ps;
    }

    hueLut->SetTableRange (min_stress, max_stress);
    mapper->SetScalarRange(min_stress, max_stress);

    hueLut->Modified();
    principalStresses_cells->Modified();
    firstPrincipalStress_cells->Modified();
}

