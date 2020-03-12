#include "backgroundworker.h"
#include <algorithm>

icy::BackgroundWorker::BackgroundWorker(ImplicitModel4 *model, QObject *)
{
    this->model = model;
    this->start();
}

icy::BackgroundWorker::~BackgroundWorker()
{
//    terminate();
    kill = true;
    model->kill = true;
    condition.wakeOne();
    wait();
}


void icy::BackgroundWorker::toggle()
{
    if(running && timeToPause) return; // already trying to pause
    else if(running && !timeToPause) {
        timeToPause = true;
        // try to abort MKL computation
    }
    else {
        // not running
        condition.wakeOne(); // resume
    }
}

void icy::BackgroundWorker::run()
{
    char fileName[20];
    forever {
        mutex.lock();
        if (!running || timeToPause) {
            running = false;
            timeToPause = false;
            condition.wait(&mutex);
        }
        running = true;
        mutex.unlock();

        if(kill) break;
        bool aborted = model->Step();
        if(model->prms->SaveVTU) {
            // save
            sprintf(fileName, "beam_%05d.vtu", model->cf.StepNumber);
            writer->SetFileName(fileName);
            writer->SetInputData(model->mc.beam->ugrid);
            writer->Write();
        }

        if(kill) break;
        if(aborted || model->cf.StepNumber >= model->prms->MaxSteps) timeToPause = true;
        if(model->prms->MaxSolves > 0 && model->cf.TotalSolves > model->prms->MaxSolves) timeToPause = true;
        if(model->prms->StepAfterWhichFractionDetectionIsTriggered > 0 &&
                model->cf.IndenterForce < model->prms->FractionDetectionForceThreshold) {
            timeToPause = true;
            saveCSV(1);
        }
        emit stepCompleted(aborted);    // notify GUI that a step was completed
    }
}

void icy::BackgroundWorker::saveCSV(int termination_reason)
{
    if(model->beamParams->beamType != 1) return;
    std::cout << "writing file" << std::endl;

    double beamLength = model->beamParams->beamL1;
    double beamWidth = model->beamParams->beamA;
    double beamThickness = model->beamParams->beamThickness;
    double volume = beamLength * beamWidth * beamThickness;
    double indenterSize = model->beamParams->IndenterSize;

    char fileName[20];
    sprintf(fileName, "%05d.csv", model->prms->InstanceNumber);
    myfile.open (fileName, std::fstream::out | std::fstream::trunc);
    myfile << termination_reason << std::endl;
    myfile << "length, " << beamLength << endl;
    myfile << "width, " << beamWidth << endl;
    myfile << "thickness, " << beamThickness << endl;
    myfile << "volume, " << volume << endl;

    // determine max force
    auto result = std::max_element(model->allFrames.begin(), model->allFrames.end(),
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
    int N = model->allFrames.size();
    for(int i=0;i<N;i++) {
        FrameInfo &fi = model->allFrames[i];
        myfile << fi.StepNumber << "," << fi.SimulationTime << "," << fi.IndenterForce << endl;;
    }
    myfile.close();
}
