
#include <cmath>
#include "face.h"

icy::Face::Face()
{

}


double icy::Face::area()
{
    auto cross = [](double e0x, double e0y, double e0z,
            double e1x, double e1y, double e1z,
            double &x, double &y, double &z) {
                x = -(e0z * e1y) + e0y * e1z;
                y = e0z * e1x - e0x * e1z;
                z = -(e0y * e1x) + e0x * e1y;
    };

    double tx0 = vrts[0]->x0;
    double ty0 = vrts[0]->y0;
    double tz0 = vrts[0]->z0;
    double tx1 = vrts[1]->x0;
    double ty1 = vrts[1]->y0;
    double tz1 = vrts[1]->z0;
    double tx2 = vrts[2]->x0;
    double ty2 = vrts[2]->y0;
    double tz2 = vrts[2]->z0;

    double x, y, z;
    cross(tx1 - tx0, ty1 - ty0, tz1 - tz0, tx2 - tx0, ty2 - ty0, tz2 - tz0,x,y,z);
    double halfMag = sqrt(x * x + y * y + z * z)/2;
    return halfMag;
}
