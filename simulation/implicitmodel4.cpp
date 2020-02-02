#include "implicitmodel4.h"
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
    linearSystem.csrd.ClearDynamic();

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
    tcf0.CopyFrom(cf); // copy current frame data

    // position the indenter
    for(auto &nd : mc.indenter->nodes) {nd.cz = nd.tz = nd.z0 - tcf0.SimulationTime*prms.IndentationVelocity;}
}


void icy::ImplicitModel4::_addCollidingNodesToStructure(){}

bool icy::ImplicitModel4::_checkDamage(){}

bool icy::ImplicitModel4::_checkDivergence(){}

void icy::ImplicitModel4::_XtoDU(){}

void icy::ImplicitModel4::_adjustTimeStep(){}

void icy::ImplicitModel4::_acceptFrame(){}


bool icy::ImplicitModel4::Step()
{
    if(!isReady) _prepare();
    _beginStep();

//    for(auto &nd : mc.activeNodes) nd->InferTentativeValues(tcf0.TimeStep, prms.NewmarkBeta, prms.NewmarkGamma);

    linearSystem.csrd.ClearDynamic();
    mc.ConstructBVH();                  // this should _not_ be called on every step

    mc.bvh.Traverse();  // traverse BVH
//    std::cout << "broad list size " << icy::BVHT::broad_list.size() << std::endl;

    std::this_thread::sleep_for (std::chrono::milliseconds(500));
    cf.SimulationTime += prms.InitialTimeStep;
    cf.StepNumber++;

    return false; // step not aborted

}
