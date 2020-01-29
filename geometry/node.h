#ifndef NODE_H
#define NODE_H

#include <vector>
#include "face.h"
namespace icy { class Node; class Face;}

class icy::Node
{
public:
    int id, altId = -1, globalNodeId;
    double x0, y0, z0;       // undisplaced position
    double cx, cy, cz;       // current position

    double ux, uy, uz;       // displacement
    double vx, vy, vz;       // velocity
    double ax, ay, az;       // acceleration
    double fx, fy, fz;       // force acting on node

    // tentative
    double unx, uny, unz;       // displacement at timestep t+h, at iteration n
    double vnx, vny, vnz;       // velocity at timestep t+h, at iteration n
    double anx, any, anz;       // acceleration at timestep t+h, at iteration n
    double tx, ty, tz;          // new position
    double dux, duy, duz;       // un - u

    bool anchored;           // anchored nodes do not contribute to the stiffness matrix
    bool isSurface;         // used to build a list of surface elements

    std::vector<icy::Face*> faces;     // Adjacent faces, if any, for collisions


    Node();
    void Initialize(double x, double y, double z, int id);
    void AcceptTentativeValues(double h);
    void InferTentativeValues(double h, double beta, double gamma);

};

#endif // NODE_H
