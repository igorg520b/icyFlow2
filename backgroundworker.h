#ifndef BACKGROUNDWORKER_H
#define BACKGROUNDWORKER_H

#include <QObject>
#include <QThread>
#include <QMutex>
#include <QWaitCondition>
#include <QString>
#include "simulation/implicitmodel4.h"

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

signals:
    void stepCompleted(bool aborted);

protected:
    void run() override;

private:
    QMutex mutex;
    QWaitCondition condition;
    ImplicitModel4 *model;
};

#endif // BACKGROUNDWORKER_H
