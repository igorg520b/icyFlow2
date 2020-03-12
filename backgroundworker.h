#ifndef BACKGROUNDWORKER_H
#define BACKGROUNDWORKER_H

#include <QObject>
#include <QThread>
#include <QMutex>
#include <QWaitCondition>
#include <QString>
#include <vtkXMLUnstructuredGridWriter.h>
#include "simulation/implicitmodel4.h"
#include <fstream>

namespace icy { class BackgroundWorker; }

class icy::BackgroundWorker : public QThread
{
    Q_OBJECT
public:
    BackgroundWorker(ImplicitModel4 *model, QObject *parent = nullptr);
    ~BackgroundWorker();
    bool timeToPause = false;
    bool running = false;
    void toggle();  // pause if running; resume if paused

private:
    vtkNew<vtkXMLUnstructuredGridWriter> writer;

signals:
    void stepCompleted(bool aborted);

protected:
    void run() override;

private:
    QMutex mutex;
    QWaitCondition condition;
    ImplicitModel4 *model;
    bool kill = false;
    void saveCSV(int termination_reason);
    fstream myfile;

};

#endif // BACKGROUNDWORKER_H
