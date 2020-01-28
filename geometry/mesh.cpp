#include "mesh.h"

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


void icy::Mesh::DetectSurfaces(bool anchorsides)
{

}

void icy::Mesh::IdentifySurfaceElements()
{

}


void icy::Mesh::ConnectFaces()
{

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
