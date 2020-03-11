#include "backgroundworker.h"

icy::BackgroundWorker::BackgroundWorker(ImplicitModel4 *model, QObject *parent)
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
                model->cf.IndenterForce < model->prms->FractionDetectionForceThreshold) timeToPause = true;
        emit stepCompleted(aborted);    // notify GUI that a step was completed
    }
}

void icy::BackgroundWorker::saveCSV(int termination_reason)
{
    if(model->beamParams->beamType != 1) return;

    char fileName[20];
    sprintf(fileName, "%05d.csv", model->prms->InstanceNumber);
    myfile.open ("results.csv", std::fstream::out | std::fstream::trunc);
    myfile << termination_reason << std::endl;
    myfile << "length, " << model->beamParams->beamL1 << endl;
    myfile << "width, " << model->beamParams->beamA << endl;
    myfile << "thickness, " << model->beamParams->beamThickness << endl;

    // determine max force

    // calculate "s" distance

    // calculate flexural strength

    int N = model->allFrames.size();
    for(int i=0;i<N;i++) {
        FrameInfo &fi = model->allFrames[i];
        //    myfile << cf.StepNumber << "," << cf.SimulationTime << "," << cf.IndenterForce;
    }



    myfile.close();


}
