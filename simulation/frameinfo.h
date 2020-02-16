#ifndef FRAMEINFO_H
#define FRAMEINFO_H

#include <iostream>

namespace icy { class FrameInfo; }

class icy::FrameInfo
{
public:
    // recorded values
    double IndenterForce;
    double extensometerDisplacements[5]={};

    bool explodes, diverges; // for time step adjustment algorithm

    // geometry
    int nCZFailedThisStep;
    int nCZDamagedThisStep;

    int nCZFailedTotal;
    int nCZDamagedTotal;

    int nCollisions;
    int nCZ_Initial;
    int nActiveNodes;


    // time
    int StepNumber;
    double SimulationTime;
    int TimeScaleFactor;
    const int Parts = 1024*32;
    unsigned long SimulationIntegerTime;  // // measures time in 1/Parts intervals of InitialTimeStep
    double TimeStep;// time difference between current frame and last frame
    int StepsWithCurrentFactor; // time steps left with current factor (if TSF > 1)
//    int TimeScaleFactorThisStep = 1;    // Time scale used for this step


    // solution analysis
    int IterationsPerformed;
    int AttemptsTaken;
    double RelativeError;
    double Error0;
    bool ConvergenceReached;

    void Print() {
        std::cout << "st " << StepNumber;
        std::cout << ";it " << IterationsPerformed;
        std::cout << "; div " << diverges;
        std::cout << "; expl " << explodes;
        std::cout << "; cvg " << ConvergenceReached;
        std::cout << "; trf " << TimeScaleFactor;
        std::cout << "; ts " << TimeStep;
        std::cout << std:: endl;
    }

    void Reset() {
        nActiveNodes = nCZ_Initial = nCZFailedThisStep = nCZDamagedThisStep = nCollisions = 0;
        nCZFailedTotal = nCZDamagedTotal = 0;
        TimeStep = SimulationTime = 0;
        SimulationIntegerTime = StepsWithCurrentFactor = 0;
        StepNumber = -1;
        TimeScaleFactor = 1;
        ConvergenceReached = false;
        IndenterForce = 0;
    }

    void CopyFrom(FrameInfo &other) {
        Reset();
        StepNumber = other.StepNumber;
        SimulationTime = other.SimulationTime;
        SimulationIntegerTime = other.SimulationIntegerTime;
        StepsWithCurrentFactor = other.StepsWithCurrentFactor;
        TimeScaleFactor = other.TimeScaleFactor;
        nCZ_Initial = other.nCZ_Initial;
        nCZFailedTotal = other.nCZFailedTotal;
        nCZDamagedTotal = other.nCZDamagedTotal;
    }

    void IncrementTime(double initialTimestep)
    {
        // we use integer increments to avoid dealing with doubles
        int ticks = Parts / TimeScaleFactor;
        SimulationIntegerTime += ticks;
        TimeStep = initialTimestep / (double) TimeScaleFactor;
        SimulationTime = initialTimestep * (double)SimulationIntegerTime / Parts;
        StepNumber++;
    }

};

#endif // FRAMEINFO_H
