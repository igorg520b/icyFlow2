#ifndef FACE_H
#define FACE_H

#include "node.h"
#include "element.h"

namespace icy {
class Face;
class Node;
class Element;
}

class icy::Face
{
public:
    icy::Node *vrts[3];
    icy::Element *elem = nullptr;   // element to which the face belongs (null unless CZs are inserted)
    int id;                         // sequential id
    int globalFaceId;               // in global array
    int granule;                    // to which granule the element belongs
    int tag;                        // surface partition id
    bool exposed = true;            // is this an outside surface
    bool created = false;           // got exposed at simulation time due to fracture
    double pnorm;

    Face();
    double area();
};

#endif // FACE_H
