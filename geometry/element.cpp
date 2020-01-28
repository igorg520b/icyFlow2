#include "element.h"
#include <stdexcept>


icy::Element::Element()
{
}

void icy::Element::FindAdjFaces()
{
    adjFaces.clear();
    for(int i=0;i<4;i++)
        for(icy::Face* const& f : (vrts[i]->faces))
            adjFaces.insert(f);
}


double icy::Element::volume()
{
    double x1, x2, x3, x4, y1, y2, y3, y4, z1, z2, z3, z4;
    double x12, x13, x14, x23, x24, x34, x21, x31, x32, x42, x43, y12, y13, y14, y23, y24, y34;
    double y21, y31, y32, y42, y43, z12, z13, z14, z23, z24, z34, z21, z31, z32, z42, z43;
    double Jdet;
    x1 = vrts[0]->x0; y1 = vrts[0]->y0; z1 = vrts[0]->z0;
    x2 = vrts[1]->x0; y2 = vrts[1]->y0; z2 = vrts[1]->z0;
    x3 = vrts[2]->x0; y3 = vrts[2]->y0; z3 = vrts[2]->z0;
    x4 = vrts[3]->x0; y4 = vrts[3]->y0; z4 = vrts[3]->z0;

    x12 = x1 - x2; x13 = x1 - x3; x14 = x1 - x4; x23 = x2 - x3; x24 = x2 - x4; x34 = x3 - x4;
    x21 = -x12; x31 = -x13; x32 = -x23; x42 = -x24; x43 = -x34;
    y12 = y1 - y2; y13 = y1 - y3; y14 = y1 - y4; y23 = y2 - y3; y24 = y2 - y4; y34 = y3 - y4;
    y21 = -y12; y31 = -y13; y32 = -y23; y42 = -y24; y43 = -y34;
    z12 = z1 - z2; z13 = z1 - z3; z14 = z1 - z4; z23 = z2 - z3; z24 = z2 - z4; z34 = z3 - z4;
    z21 = -z12; z31 = -z13; z32 = -z23; z42 = -z24; z43 = -z34;
    Jdet = x21 * (y23 * z34 - y34 * z23) + x32 * (y34 * z12 - y12 * z34) + x43 * (y12 * z23 - y23 * z12);
    double V = Jdet / 6.0;

    if(V<=0) throw std::runtime_error("element volume: V<=0");
    return V; // supposed to be positive
}

