#include "mesh.h"
#include <stdexcept>

icy::Mesh::Mesh()
{

}

icy::Mesh::~Mesh()
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
/*
            surfaceFragments.Clear();
            // sequential ids for faces
            for (int i = 0; i < faces.Count; i++) faces[i].id = i;

            for (int i = 0; i < 6; i++) surfaceFragments.Add(new SurfaceFragment() { id = i });
            surfaceFragments[0].name = "top";
            surfaceFragments[1].name = "bottom";
            surfaceFragments[2].name = "xmax";
            surfaceFragments[3].name = "xmin";
            surfaceFragments[4].name = "ymax";
            surfaceFragments[5].name = "ymin";

            foreach (Face f in faces)
            {
                if (f.vrts.All(nd => nd.z0 == zmax)) surfaceFragments[0].faces.Add(f.id);
                else if (f.vrts.All(nd => nd.z0 == zmin)) surfaceFragments[1].faces.Add(f.id);
                else if (f.vrts.All(nd => nd.x0 == xmax)) surfaceFragments[2].faces.Add(f.id);
                else if (f.vrts.All(nd => nd.x0 == xmin)) surfaceFragments[3].faces.Add(f.id);
                else if (f.vrts.All(nd => nd.y0 == ymax)) surfaceFragments[4].faces.Add(f.id);
                else if (f.vrts.All(nd => nd.y0 == ymin)) surfaceFragments[5].faces.Add(f.id);
            }
            foreach (SurfaceFragment sf in surfaceFragments) { sf.allFaces = faces; sf.ComputeArea(); }
            //            surfaceFragments[0].sensor = true;
            //            surfaceFragments[0].role = SurfaceFragment.SurfaceRole.Anchored;
            //            surfaceFragments[0].dz = -0.23;
            //           surfaceFragments[1].role = SurfaceFragment.SurfaceRole.Anchored;
            if (anchorsides)
            {
                surfaceFragments[2].role = SurfaceFragment.SurfaceRole.Anchored;
                surfaceFragments[3].role = SurfaceFragment.SurfaceRole.Anchored;
                surfaceFragments[4].role = SurfaceFragment.SurfaceRole.Anchored;
                surfaceFragments[5].role = SurfaceFragment.SurfaceRole.Anchored;
            }
*/
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

    for(auto &elem : elems)
        if(elem.vrts[0]->isSurface ||
                elem.vrts[1]->isSurface ||
                elem.vrts[2]->isSurface ||
                elem.vrts[3]->isSurface) elem.isSurface = true;
        else elem.isSurface = false;
}


void icy::Mesh::ConnectFaces()
{
throw std::runtime_error("not implemented");
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
}

void icy::Mesh::UpdateUGrid()
{
    throw std::runtime_error("not implemented");
}
