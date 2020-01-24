#ifndef MODELPRMS_H
#define MODELPRMS_H

#include <QObject>

namespace icy {
class ModelPrms;
}

class icy::ModelPrms : public QObject
{
    Q_OBJECT
public:
    ModelPrms();
};

#endif // MODELPRMS_H
