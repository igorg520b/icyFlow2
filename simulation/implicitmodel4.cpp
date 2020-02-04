#include "implicitmodel4.h"
#include "numbercrunching.h"
#include "bvh/bvht.h"
#include <thread>         // std::this_thread::sleep_for
#include <chrono>         // std::chrono::seconds

icy::ImplicitModel4::ImplicitModel4()
{
    Clear();
}

void icy::ImplicitModel4::Clear()
{
//    allFrames.clear();
    cf.Reset();
    tcf0.Reset();
    isReady = false;
    mc.Clear();
}


void icy::ImplicitModel4::_prepare()
{
//    allFrames.push_back(cf);
    // re-create static contents of the linear system
    mc.Prepare();   // populate activeNodes
    linearSystem.csrd.ClearStatic();
//    linearSystem.csrd.ClearDynamic();

    // add entries for elastic elements to mesh collection
    // in general, it would be better have a collection for deformables
    for(auto const &elem : mc.beam->elems)
    {
        for(int i=0;i<4;i++) {
            Node *nd1 = elem.vrts[i];
            for(int j=0;j<4;j++) {
                Node *nd2 = elem.vrts[j];
                if(!nd1->anchored && !nd2->anchored && nd2->altId >= nd1->altId)
                    linearSystem.csrd.AddStatic(nd1->altId, nd2->altId);
            }
        }
    }
    isReady = true;
}

void icy::ImplicitModel4::_beginStep()
{
    // prepare tentative entry about the simulation step
    tcf0.CopyFrom(cf); // copy current frame data
    tcf0.IncrementTime(prms.InitialTimeStep);
    tcf0.ActiveNodes = (int)mc.activeNodes.size();

    // position the indenter
    for(auto &nd : mc.indenter->nodes) {
        nd.unz = nd.uz = - tcf0.SimulationTime*prms.IndentationVelocity;
        nd.cz = nd.tz = nd.z0 + nd.uz;
    }
}


void icy::ImplicitModel4::_addCollidingNodesToStructure(){}

bool icy::ImplicitModel4::_checkDamage(){}

bool icy::ImplicitModel4::_checkDivergence(){}

void icy::ImplicitModel4::_XtoDU(){}

void icy::ImplicitModel4::_adjustTimeStep(){}

void icy::ImplicitModel4::_acceptFrame()
{
    cf.CopyFrom(tcf0);
}

void icy::ImplicitModel4::_assemble()
{
    for(auto &nd : mc.allNodes) nd->fx = nd->fz = nd->fy = 0;

    // assemble elements

    // assemble cohesive zones
}

void icy::ImplicitModel4::_narrowPhase()
{
}

void icy::ImplicitModel4::_collisionResponse()
{

}

void icy::ImplicitModel4::_transferUpdatedState()
{

}


bool icy::ImplicitModel4::Step()
{
    std::this_thread::sleep_for(std::chrono::milliseconds(500)); // for testing

    if(!isReady) _prepare();

    do {
        _beginStep();

        for(auto &nd : mc.activeNodes) nd->InferTentativeValues(tcf0.TimeStep, prms.NewmarkBeta, prms.NewmarkGamma);

        mc.ConstructBVH();                                              // this should _not_ be called on every step
        mc.bvh.Traverse();                                              // traverse BVH
        NumberCrunching::NarrowPhase(icy::BVHT::broad_list, mc);        // narrow phase

        linearSystem.csrd.ClearDynamic();
        _addCollidingNodesToStructure();  // add colliding nodes to structure
        linearSystem.CreateStructure();

        _assemble(); // prepare and compute forces
        explodes = _checkDamage();
        if(!explodes) {
            _collisionResponse();
//            linearSystem.Solve();
            diverges = _checkDivergence();
            _XtoDU();
            tcf0.IterationsPerformed++;
        }
    } while(false);

    cf.TimeScaleFactorThisStep = tcf0.TimeScaleFactor; // record what TSF was used for this step
    _adjustTimeStep();
    _transferUpdatedState();
    _acceptFrame();


    return false; // step not aborted

}
