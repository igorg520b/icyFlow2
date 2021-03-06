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
            model->writeCSV(1);
        }
        emit stepCompleted(aborted);    // notify GUI that a step was completed
    }
}


