#ifndef CZ_H
#define CZ_H

#include "node.h"
#include "face.h"

namespace icy {
class CZ;
class Element;
class Node;
class Face;
}

class icy::CZ
{
public:
    icy::Node* vrts[6];                 // References to element vertices
    icy::Face* faces[6];                 // References to element vertices
    double pmax[3]={};
    double tmax[3]={};
    bool failed = false;
//    double avgDn, avgDt, avgTn, avgTt; // average traction-separations for subsequent analysis
//    double maxAvgDn, maxAvgDt;

    CZ();
};


#endif // CZ_H
