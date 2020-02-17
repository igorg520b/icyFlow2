#ifndef BEAMPARAMS_H
#define BEAMPARAMS_H

#include <QObject>

class BeamParams : public QObject
{
    Q_OBJECT

    Q_PROPERTY(double beamA MEMBER beamA NOTIFY propertyChanged)
    Q_PROPERTY(double beamB MEMBER beamB NOTIFY propertyChanged)
    Q_PROPERTY(double beamL1 MEMBER beamL1 NOTIFY propertyChanged)
    Q_PROPERTY(double beamL2 MEMBER beamL2 NOTIFY propertyChanged)
    Q_PROPERTY(double beamGap MEMBER beamGap NOTIFY propertyChanged)
    Q_PROPERTY(double beamMargin MEMBER beamMargin NOTIFY propertyChanged)
    Q_PROPERTY(double beamThickness MEMBER beamThickness NOTIFY propertyChanged)
    Q_PROPERTY(double CharacteristicLengthMax MEMBER CharacteristicLengthMax NOTIFY propertyChanged)
    Q_PROPERTY(double CharacteristicLengthIndenter MEMBER CharacteristicLengthIndenter NOTIFY propertyChanged)
    Q_PROPERTY(double RefinementMultiplier MEMBER RefinementMultiplier NOTIFY propertyChanged)
    Q_PROPERTY(double IndenterSize MEMBER IndenterSize NOTIFY propertyChanged)

public:
    explicit BeamParams(QObject *parent = nullptr);

    double beamA = 0.4;     // a
    double beamB = 0.4;     // b
    double beamL1 = 0.83;
    double beamL2 = 1.3;
    double beamGap = 0.1;   // c
    double beamMargin = 0.35; // d
    double beamThickness = 0.5; // h
    double CharacteristicLengthMax = 0.08;//0.2;
    double CharacteristicLengthIndenter = 0.02;
    double RefinementMultiplier = 0.04;//0.07;
    double IndenterSize = 0.15;

signals:
    void propertyChanged();
};

#endif // BEAMPARAMS_H
