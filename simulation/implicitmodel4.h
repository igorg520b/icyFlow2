#ifndef IMPLICITMODEL4_H
#define IMPLICITMODEL4_H

#include <vector>
#include "geometry/meshcollection.h"
#include "modelprms.h"
#include "frameinfo.h"
#include "linearsystem.h"
#include "bvh/bvht.h"

namespace icy {
class ImplicitModel4;
}

class icy::ImplicitModel4
{
public:
    MeshCollection mc;
    ModelPrms prms;
    FrameInfo cf, tcf0;
//    std::vector<FrameInfo> allFrames;
    LinearSystem linearSystem;

    bool isReady = false;

    ImplicitModel4();
    void Clear();
    bool Step();    // return true if aborted

private:
    bool explodes, diverges;

    void _prepare();    // this is called once
    void _beginStep();
    void _addCollidingNodesToStructure();
    bool _checkDamage();
    bool _checkDivergence();
    void _XtoDU();

    void _adjustTimeStep();
    void _acceptFrame();
    void _assemble();
    void _transferUpdatedState();

};

#endif // IMPLICITMODEL4_H
