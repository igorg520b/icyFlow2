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

 //   for(auto &sf : surfaceFragments) delete sf;
 //   for(auto &fc : faces) delete fc;
 //   for(auto &cz : czs) delete cz;
 //   for(auto &elem : elems) delete elem;
 //   for(auto &nd : nodes) delete nd;
}

void icy::Mesh::ComputeBoundingBox()
{

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


