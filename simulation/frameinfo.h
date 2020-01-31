#ifndef FRAMEINFO_H
#define FRAMEINFO_H

#include <QObject>

namespace icy { class FrameInfo; }

class icy::FrameInfo : public QObject
{
    Q_OBJECT
    Q_PROPERTY(int StepNumber MEMBER StepNumber NOTIFY propertyChanged)


public:

    // geometry
    int nCZFailed;
    int nCZDamaged;
    int nCollisions;
    int nCZ_Initial;


    // time
    int StepNumber;
    double SimulationTime;
    int TimeScaleFactor;
    const int Parts = 1024*32;
    int SimulationIntegerTime;  // // measures time in 1/Parts intervals of InitialTimeStep
    double TimeStep;// time difference between current frame and last frame
    int StepsWithCurrentFactor; // time steps left with current factor (if TSF > 1)
    int TimeScaleFactorThisStep = 1;    // Time scale used for this step

    double IndenterForce;

    // solution analysis
    int IterationsPerformed;
    int AttemptsTaken;
    double RelativeError;
    double Error0;

    void Reset() {
        nCZ_Initial = nCZFailed = nCZDamaged = nCollisions = 0;
        TimeStep = SimulationTime = 0;
        SimulationIntegerTime = StepsWithCurrentFactor = 0;
        StepNumber = -1;
        TimeScaleFactor = 1;
        TimeScaleFactorThisStep = 1;

    }

    int MaxIntTimestep() { return Parts - SimulationIntegerTime % Parts; }

    void CopyFrom(FrameInfo &other) {
        Reset();
        StepNumber = other.StepNumber;
        SimulationTime = other.SimulationTime;
        SimulationIntegerTime = other.SimulationIntegerTime;
        StepsWithCurrentFactor = other.StepsWithCurrentFactor;
        TimeScaleFactor = other.TimeScaleFactor;
        nCZ_Initial = other.nCZ_Initial;
        nCZFailed = other.nCZFailed;
    }


signals:
    void propertyChanged();
};

#endif // FRAMEINFO_H
