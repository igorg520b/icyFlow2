#ifndef FRAMEINFO_H
#define FRAMEINFO_H

#include <QObject>

class FrameInfo : public QObject
{
    Q_OBJECT
public:
    explicit FrameInfo(QObject *parent = nullptr);

signals:

};

#endif // FRAMEINFO_H
