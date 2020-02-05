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
        if(kill) break;
        if(aborted || model->cf.StepNumber >= model->prms.MaxSteps) timeToPause = true;
        emit stepCompleted(aborted);    // notify GUI that a step was completed
    }
}


