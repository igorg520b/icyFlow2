#ifndef FRAMEINFO_H
#define FRAMEINFO_H

#include <QObject>

namespace icy {
class FrameInfo;
}

class icy::FrameInfo : public QObject
{
    Q_OBJECT
public:
    explicit FrameInfo(QObject *parent = nullptr);

signals:
    void propertyChanged();
};

#endif // FRAMEINFO_H
