#include "implicitmodel4.h"
#include "numbercrunching.h"
#include "bvh/bvht.h"
#include <thread>         // std::this_thread::sleep_for
#include <chrono>         // std::chrono::seconds
//#include <omp.h>


icy::ImplicitModel4::ImplicitModel4()
{
    Clear();
    allFrames.reserve(1000);
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
    icy::NumberCrunching::InitializeConstants(*prms);

    allFrames.push_back(cf);
    // re-create static contents of the linear system
    mc.Prepare();   // populate activeNodes
    _updateStaticStructure();

    if(cf.StepNumber == 0) cf.nCZ_Initial = (int)mc.allCZs.size();
    isReady = true;

    //extensometers

    if(beamParams->beamType == 0)
    {
        // obb tree
        filter1->SetInputDataObject(mc.beam->ugrid);
        filter1->Update();
        obbTree->SetDataSet(filter1->GetOutput());
    }
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
    //    std::cout << "static csrd size " << linearSystem.csrd.staticCount() << std::endl;
}

void icy::ImplicitModel4::_beginStep()
{
    if(cf.StepNumber == 0) cf.nCZ_Initial = mc.activeCZs.size();

        // prepare tentative entry about the simulation step
    tcf0.CopyFrom(cf); // copy current frame data
    tcf0.IncrementTime(prms->InitialTimeStep);
    tcf0.nActiveNodes = (int)mc.activeNodes.size();
    tcf0.ConvergenceReached = false;
    tcf0.IterationsPerformed = 0;

    // position the indenter
    for(auto &nd : mc.indenter->nodes) {
        nd.unz = nd.uz = - tcf0.SimulationTime*prms->IndentationVelocity;
        nd.cz = nd.tz = nd.z0 + nd.uz;
    }

    // initial displacement guess for active nodes
    double h = tcf0.TimeStep;
    for(auto &nd : mc.activeNodes) {
        nd->dux = nd->vx*h;
        nd->duy = nd->vy*h;
        nd->duz = nd->vz*h;
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
    double dn_damaged = (double)(tcf0.nCZDamagedThisStep-cf.nCZDamagedTotal) / tcf0.nCZ_Initial;
    double dn_failed = (double)(tcf0.nCZFailedThisStep) / tcf0.nCZ_Initial;
    bool result = (dn_damaged > prms->maxDamagePerStep || dn_failed > prms->maxFailPerStep);
    if(result)
        std::cout << "damaged " << tcf0.nCZDamagedThisStep << "; f " << tcf0.nCZFailedThisStep << ";i "  << tcf0.nCZ_Initial << std::endl;
    return result;
}

bool icy::ImplicitModel4::_checkDivergence()
{
    bool divergence;
    double cutoff = prms->ConvergenceCutoff; // somewhat arbitrary constant
    double norm = linearSystem.NormOfDx();
    if (tcf0.IterationsPerformed == 0)
    {
        tcf0.Error0 = norm;
        if (tcf0.Error0 == 0) tcf0.Error0 = prms->ConvergenceEpsilon;
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
        if (tcf0.RelativeError <= prms->ConvergenceEpsilon) tcf0.ConvergenceReached = true;
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
    const int holdFactorDelay = 4;
    cf.AttemptsTaken++;
    if (tcf0.explodes)
    {
        cf.TimeScaleFactor *= 2;
        if (cf.TimeScaleFactor > cf.Parts) cf.TimeScaleFactor = cf.Parts;
        std::cout << "_adjustTimeStep: explodes, new sf: " << cf.TimeScaleFactor << std::endl;
        cf.StepsWithCurrentFactor = holdFactorDelay;
    }
    else if ((tcf0.diverges || !tcf0.ConvergenceReached) &&
             tcf0.TimeScaleFactor < cf.Parts && prms->maxIterations > 1)
    {
        // does not converge
        cf.TimeScaleFactor *= 2;
        cf.StepsWithCurrentFactor = holdFactorDelay;
        std::cout << "_adjustTimeStep: diverges " << tcf0.diverges << "; convergenceReached " << tcf0.ConvergenceReached << "; factor " << cf.TimeScaleFactor << std::endl;
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
    }
}

void icy::ImplicitModel4::_acceptFrame()
{
    cf.CopyFrom(tcf0);
    cf.nCZFailedTotal += tcf0.nCZFailedThisStep;
    cf.nCZDamagedTotal = tcf0.nCZDamagedThisStep;

    size_t N = mc.activeNodes.size();
#pragma omp parallel for
    for(size_t i=0;i<N;i++) {
        Node *nd = mc.activeNodes[i];
        nd->InferTentativeValues(tcf0.TimeStep, prms->NewmarkBeta, prms->NewmarkGamma);
        nd->AcceptTentativeValues(tcf0.TimeStep);
    }

    // accept new state variables in CZs
    for(auto &cz : mc.activeCZs) {
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

    // read indenter force
    double fx, fy, fz;
    fx = fy = fz = 0;
    for(auto &nd : mc.indenter->nodes) {
        fx += nd.fx;
        fy += nd.fy;
        fz += nd.fz;
    }
    cf.IndenterForce = sqrt(fx*fx+fy*fy+fz*fz);

    if(beamParams->beamType == 0)
    {
        // "measure" the vertical deflection at extensometer locations
        // extensometers
        filter1->Update();
        obbTree->BuildLocator();
        for(int i=0;i<5;i++) {
            int result = obbTree->IntersectWithLine(&extPointsS[i*3], &extPointsE[i*3],extPoints, extIdList);
            if(result != 0) {
                double pt[3] = {};
                extPoints->GetPoint(0, pt);
                //            std::cout << i << "; " << pt[0] << ", " << pt[1] << ", " << pt[2] << std::endl;
                cf.extensometerDisplacements[i]=pt[2]-beamParams->beamThickness;
            } else cf.extensometerDisplacements[i] = 0;
        }

        myfile.open ("results.csv", std::fstream::in | std::fstream::out | std::fstream::app);
        myfile << cf.StepNumber << "," << cf.SimulationTime << "," << cf.IndenterForce;
        for(int i=0;i<5;i++) myfile << "," << cf.extensometerDisplacements[i];
        myfile << std::endl;
        myfile.close();
    } else if(beamParams->beamType == 1) {
        // cantilever beam

    }

    allFrames.push_back(cf);
}

void icy::ImplicitModel4::_assemble()
{
    linearSystem.CreateStructure();
    for(auto &nd : mc.allNodes) nd->fx = nd->fz = nd->fy = 0;

    // assemble elements
    if(tcf0.TimeStep == 0) throw std::runtime_error("time step is zero");
    NumberCrunching::AssembleElems(linearSystem, mc.elasticElements, tcf0.TimeStep);

    // assemble cohesive zones
    NumberCrunching::AssembleCZs(linearSystem, mc.activeCZs, tcf0.nCZFailedThisStep, tcf0.nCZDamagedThisStep);

    // assemble collisions
    NumberCrunching::CollisionResponse(linearSystem, prms->DistanceEpsilon, prms->penaltyK);

    if(beamParams->beamType == 1) {
        double m = (mc.beam->ymax - mc.beam->ymin)/100;
        NumberCrunching::ForceBox(linearSystem, mc.activeNodes, prms->penaltyK,
          mc.beam->xmin-m, mc.beam->xmax+m, mc.beam->ymin-m, mc.beam->ymax+m);
    }
}

bool icy::ImplicitModel4::Step()
{
    kill = false;
    if(!isReady) _prepare();
    _beginStep();

    do {

        for(auto &nd : mc.activeNodes)
            nd->InferTentativeValues(tcf0.TimeStep, prms->NewmarkBeta, prms->NewmarkGamma);

        mc.ConstructBVH();                                              // this should _not_ be called on every step
        mc.bvh.Traverse();                                              // traverse BVH
        NumberCrunching::NarrowPhase(icy::BVHT::broad_list, mc);        // narrow phase
        _addCollidingNodesToStructure();  // add colliding nodes to structure

        _assemble(); // prepare and compute forces
        tcf0.explodes = _checkDamage();

        linearSystem.Solve();
        cf.TotalSolves++;
        if(kill) return true;
        tcf0.diverges = _checkDivergence();
        _XtoDU();
        tcf0.IterationsPerformed++;
//        tcf0.Print();
    } while(tcf0.IterationsPerformed < prms->minIterations ||
      (!tcf0.explodes && !tcf0.diverges && !tcf0.ConvergenceReached && tcf0.IterationsPerformed < prms->maxIterations));

    _adjustTimeStep();

    if (prms->maxIterations == 1 ||
        tcf0.TimeScaleFactor == tcf0.Parts ||
        (!tcf0.ConvergenceReached && tcf0.TimeScaleFactor >= tcf0.Parts) ||
        (!tcf0.explodes && !tcf0.diverges && tcf0.ConvergenceReached)) _acceptFrame();

    return false; // step not aborted
}


void icy::ImplicitModel4::writeCSV(int termination_reason)
{
    if(beamParams->beamType != 1) return;
    std::cout << "writing file" << std::endl;

    double beamLength = beamParams->beamL1;
    double beamWidth = beamParams->beamA;
    double beamThickness = beamParams->beamThickness;
    double volume = beamLength * beamWidth * beamThickness;
    double indenterSize = beamParams->IndenterSize;

    char fileName[20];
    sprintf(fileName, "%d.csv", prms->InstanceNumber);
    myfile.open (fileName, std::fstream::out | std::fstream::trunc);
    myfile << termination_reason << std::endl;
    myfile << "length, " << beamLength << endl;
    myfile << "width, " << beamWidth << endl;
    myfile << "thickness, " << beamThickness << endl;
    myfile << "volume, " << volume << endl;
    myfile << "elems, " << mc.beam->elems.size() << endl;
    myfile << "czs, " << mc.beam->czs.size() << endl;
    myfile << "nodes, " << mc.beam->nodes.size() << endl;
    myfile << "CharacteristicLengthMax, " << beamParams->CharacteristicLengthMax << endl;
    myfile << "cz_sigma_max, " << prms->sigma_max << endl;
    myfile << "cz_tau_max, " << prms->tau_max << endl;

    // determine max force
    auto result = std::max_element(allFrames.begin(), allFrames.end(),
                                  [](const FrameInfo &a, const FrameInfo &b)
                                    {return a.IndenterForce < b.IndenterForce;});
    double maxForce = result->IndenterForce;
    myfile << "maxForce, " << maxForce << endl;

    // calculate "s" distance
    double s = (beamLength-indenterSize)/3.0;
    myfile << "s, " << s << endl;

    // calculate flexural strength
    double sigma_f = 3*s*maxForce/(beamWidth*beamThickness*beamThickness);
    myfile << "sigma_f, " << sigma_f << endl;

    myfile << "step_number, time, indenter_force" << endl;
    int N = allFrames.size();
    for(int i=0;i<N;i++) {
        FrameInfo &fi = allFrames[i];
        myfile << fi.StepNumber << "," << fi.SimulationTime << "," << fi.IndenterForce << endl;;
    }
    myfile.close();
}
