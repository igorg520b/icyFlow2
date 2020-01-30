#include "implicitmodel4.h"

icy::ImplicitModel4::ImplicitModel4()
{
    Clear();
}

void icy::ImplicitModel4::Clear()
{
    allFrames.clear();
    cf.Reset();
    tcf0.Reset();
    isReady = false;
    mc.Clear();
}


void icy::ImplicitModel4::_prepare()
{

}

void icy::ImplicitModel4::_beginStep()
{
    linearSystem.csrd.ClearStatic();    // this should be called in _prepare
    linearSystem.csrd.ClearDynamic();

    mc.ConstructBVH();                  // this should _not_ be called on every step

    // position the indenter
    for(auto &nd : mc.indenter->nodes) nd.tz = nd.z0 - tcf0.SimulationTime*prms.IndentationVelocity;

}


void icy::ImplicitModel4::_addCollidingNodesToStructure(){}

bool icy::ImplicitModel4::_checkDamage(){}

bool icy::ImplicitModel4::_checkDivergence(){}

void icy::ImplicitModel4::_XtoDU(){}

void icy::ImplicitModel4::_adjustTimeStep(){}

void icy::ImplicitModel4::_acceptFrame(){}


void icy::ImplicitModel4::Step()
{

}
