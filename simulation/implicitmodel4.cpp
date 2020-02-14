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


void icy::ImplicitModel4::_updateStaticStructure()
{
    linearSystem.csrd.ClearStatic();
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

    for(auto const &cz : mc.beam->czs)
    {
        for(int i=0;i<6;i++) {
            Node *nd1 = cz.vrts[i];
            for(int j=0;j<6;j++) {
                Node *nd2 = cz.vrts[j];
                if(!nd1->anchored && !nd2->anchored && nd2->altId >= nd1->altId)
                    linearSystem.csrd.AddStatic(nd1->altId, nd2->altId);
            }
        }
    }
}

void icy::ImplicitModel4::_prepare()
{
    icy::NumberCrunching::InitializeConstants();

//    allFrames.push_back(cf);
    // re-create static contents of the linear system
    mc.Prepare();   // populate activeNodes
    _updateStaticStructure();
//    linearSystem.csrd.ClearDynamic();

    if(cf.StepNumber == 0) cf.nCZ_Initial = (int)mc.allCZs.size();
//    std::cout << "static csrd size " << linearSystem.csrd.staticCount() << std::endl;
    isReady = true;
}

void icy::ImplicitModel4::_beginStep()
{
    // prepare tentative entry about the simulation step
    tcf0.CopyFrom(cf); // copy current frame data
    tcf0.IncrementTime(prms.InitialTimeStep);
    tcf0.nActiveNodes = (int)mc.activeNodes.size();
    tcf0.ConvergenceReached = false;
    tcf0.IterationsPerformed = 0;

    // position the indenter
    for(auto &nd : mc.indenter->nodes) {
        nd.unz = nd.uz = - tcf0.SimulationTime*prms.IndentationVelocity;
        nd.cz = nd.tz = nd.z0 + nd.uz;
    }

    // initial displacement guess for active nodes

    double h = tcf0.TimeStep;
    for(auto &nd : mc.activeNodes) {
        nd->dux = nd->vx*h;
        nd->duy = nd->vy*h;
        nd->duz = nd->vz*h;
    }

    // for testing
    for(auto &nd : mc.beam->nodes) {
        if(nd.anchored) {
            nd.unz = nd.uz = tcf0.SimulationTime*prms.IndentationVelocity;
            nd.cz = nd.tz = nd.z0 + nd.uz;
       }
   }
}


void icy::ImplicitModel4::_addCollidingNodesToStructure()
{
    linearSystem.csrd.ClearDynamic();
    int count = (int)NumberCrunching::cprList.size();
    tcf0.nCollisions = count;
    for(int k=0;k<count;k++) {
        CPResult &cpr = NumberCrunching::cprList[k];
        Node *nd = cpr.nd;
        Face *fc = cpr.fc;

        Node *nds[4] = {nd, fc->vrts[0], fc->vrts[1], fc->vrts[2]};

        for(int i=0;i<4;i++) {
            Node *n1 = nds[i];
            for(int j=0;j<4;j++) {
                Node *n2 = nds[j];
                if(!n1->anchored && !n2->anchored && (n2->altId > n1->altId))
                    linearSystem.csrd.AddDynamic(n1->altId,n2->altId);
            }
        }
    }
}

bool icy::ImplicitModel4::_checkDamage()
{

    if (tcf0.TimeScaleFactor == tcf0.Parts) return false; // can't reduce time step anyway

    double dn_damaged = (double)(tcf0.nCZDamagedThisStep) / cf.nCZ_Initial;
    double dn_failed = (double)(tcf0.nCZFailedThisStep) / tcf0.nCZ_Initial;

    bool result = (dn_damaged > prms.maxDamagePerStep || dn_failed > prms.maxFailPerStep);

//    std::cout << "_checkDamage: " << result << std::endl;
    return result;
}

bool icy::ImplicitModel4::_checkDivergence()
{
    bool divergence;
    double cutoff = prms.ConvergenceCutoff; // somewhat arbitrary constant
    double norm = linearSystem.NormOfDx();
    if (tcf0.IterationsPerformed == 0)
    {
        tcf0.Error0 = norm;
        if (tcf0.Error0 == 0) tcf0.Error0 = prms.ConvergenceEpsilon;
        tcf0.ConvergenceReached = false;
        divergence = false;
    }
    else if (norm < cutoff)
    {
        tcf0.ConvergenceReached = true;
        tcf0.RelativeError = -1;
        divergence = false;
    }
    else
    {
        if(tcf0.Error0 == 0) throw std::runtime_error("tcf0.Error0 == 0");
        tcf0.RelativeError = sqrt(norm / tcf0.Error0);
        if (tcf0.RelativeError <= prms.ConvergenceEpsilon) tcf0.ConvergenceReached = true;
        divergence = (tcf0.RelativeError > 1.01); // return true diverges
    }
    return divergence;
}

void icy::ImplicitModel4::_XtoDU()
{
    for(auto &nd : mc.activeNodes) {
        int idx3 = nd->altId*3;
        nd->dux += linearSystem.dx[idx3+0];
        nd->duy += linearSystem.dx[idx3+1];
        nd->duz += linearSystem.dx[idx3+2];
    }
}

void icy::ImplicitModel4::_adjustTimeStep()
{
    std::cout << "_adjustTimeStep: ";
    const int holdFactorDelay = 4;
    cf.AttemptsTaken++;
    if (explodes)
    {
        cf.TimeScaleFactor *= 2;
        if (cf.TimeScaleFactor > cf.Parts) cf.TimeScaleFactor = cf.Parts;
        std::cout << "explodes, new sf: " << cf.TimeScaleFactor << std::endl;
        cf.StepsWithCurrentFactor = holdFactorDelay;
    }
    else if ((diverges || !tcf0.ConvergenceReached) && tcf0.TimeScaleFactor < cf.Parts && prms.maxIterations > 1)
    {
        // does not converge
        cf.TimeScaleFactor *= 2;
        cf.StepsWithCurrentFactor = holdFactorDelay;
        std::cout << "double ts: diverges " << diverges << "; convergenceReached " << tcf0.ConvergenceReached << "; factor " << cf.TimeScaleFactor << std::endl;
    }
    else if (tcf0.ConvergenceReached && cf.TimeScaleFactor > 1)
    {
        // converges
        if (tcf0.StepsWithCurrentFactor > 0) tcf0.StepsWithCurrentFactor--;
        else
        {
            tcf0.TimeScaleFactor /= 2;
            if (tcf0.TimeScaleFactor < 1) tcf0.TimeScaleFactor = 1;
            tcf0.StepsWithCurrentFactor = 2;
        }
        std::cout << "converges " << std::endl;
    }
}

void icy::ImplicitModel4::_acceptFrame()
{
    cf.CopyFrom(tcf0);
    cf.nCZFailedTotal+=tcf0.nCZFailedThisStep;
    std::cout << "accepting frame with ts " << tcf0.TimeStep;
    for(auto &nd : mc.allNodes)
        nd->AcceptTentativeValues(tcf0.TimeStep);

    // accept new state variables in CZs
    for(auto &cz : mc.nonFailedCZs) {
        cz->AcceptTentativeValues();
        if(cz->failed) {
            cz->faces[0]->created = cz->faces[0]->exposed = true;
            cz->faces[1]->created = cz->faces[1]->exposed = true;
        }
    }

    // remove CZs that failed and update matrix structure
    if(tcf0.nCZFailedThisStep > 0) {
        mc.UpdateCZs();
        _updateStaticStructure();
    }
}

void icy::ImplicitModel4::_assemble()
{
    linearSystem.CreateStructure();
    for(auto &nd : mc.allNodes) nd->fx = nd->fz = nd->fy = 0;

    // assemble elements
    if(tcf0.TimeStep == 0) throw std::runtime_error("time step is zero");
    NumberCrunching::AssembleElems(linearSystem, mc.elasticElements, prms, tcf0.TimeStep);

    // assemble cohesive zones
//    NumberCrunching::AssembleCZs(linearSystem, mc.allCZs, prms, tcf0.nCZFailedThisStep, tcf0.nCZDamagedThisStep);

    // assemble collisions
//    NumberCrunching::CollisionResponse(linearSystem, prms.DistanceEpsilon, prms.penaltyK);
}



bool icy::ImplicitModel4::Step()
{
    kill = false;
    if(!isReady) _prepare();
    _beginStep();

    do {

        for(auto &nd : mc.activeNodes) nd->InferTentativeValues(tcf0.TimeStep, prms.NewmarkBeta, prms.NewmarkGamma);

        mc.ConstructBVH();                                              // this should _not_ be called on every step
        mc.bvh.Traverse();                                              // traverse BVH
        NumberCrunching::NarrowPhase(icy::BVHT::broad_list, mc);        // narrow phase
        _addCollidingNodesToStructure();  // add colliding nodes to structure

        _assemble(); // prepare and compute forces
        explodes = _checkDamage();
        if(kill) { std::cout << "killing "; return true; }
        linearSystem.Solve();
        if(kill) { std::cout << "killing "; return true; }
        diverges = _checkDivergence();
        _XtoDU();
        tcf0.IterationsPerformed++;
        std::cout << "st " << tcf0.StepNumber;
        std::cout << ";it " << tcf0.IterationsPerformed;
        std::cout << "; div " << diverges;
        std::cout << "; expl " << explodes;
        std::cout << "; cvg " << tcf0.ConvergenceReached;
        std::cout << "; trf " << tcf0.TimeScaleFactor;
        std::cout << "; ts " << tcf0.TimeStep << std:: endl;
    } while(tcf0.IterationsPerformed < prms.minIterations ||
      (!explodes && !diverges && !tcf0.ConvergenceReached && tcf0.IterationsPerformed < prms.maxIterations));

//    cf.TimeScaleFactorThisStep = tcf0.TimeScaleFactor; // record what TSF was used for this step
    _adjustTimeStep();

    if (prms.maxIterations == 1 ||
        tcf0.TimeScaleFactor == tcf0.Parts ||
        (!tcf0.ConvergenceReached && tcf0.TimeScaleFactor >= tcf0.Parts) ||
        (!explodes && !diverges && tcf0.ConvergenceReached)) _acceptFrame();

    return false; // step not aborted
}
