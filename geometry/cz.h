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
    double avgDn, avgDt, avgTn, avgTt; // average traction-separations for subsequent analysis
    double maxAvgDn, maxAvgDt;

    // these need to be reset before assembly (tentative values)
    bool _contact, _failed, _damaged;
    double _avgDn, _avgDt, _avgTn, _avgTt;
    double _pmax, _tmax;
    double pmax_[3], tmax_[3];
    double rhs[18];
    double lhs[18][18];

    CZ();
    void AcceptTentativeValues();
};


#endif // CZ_H
