#include "generatortool.h"
#include <QtCore>
#include <QCommandLineParser>
#include <QDebug>
#include "simulation/modelprms.h"
#include <cmath>
#include "simulation/implicitmodel4.h"
#include <fstream>
#include <gmsh.h>


int main(int argc, char *argv[])
{
    gmsh::initialize();

    QCoreApplication a(argc, argv);
    QCoreApplication::setApplicationName("icFlow");
    QCoreApplication::setApplicationVersion("2.1");


    // parse command line options
    QCommandLineParser parser;
    parser.addHelpOption();
    parser.addVersionOption();
    // An option with a value
    QCommandLineOption idxOption(QStringList() << "i" << "idx",
            QCoreApplication::translate("main", "Integer process index starting from zero <idx>"),
            QCoreApplication::translate("main", "idx"));

    QCommandLineOption idxMaxOption(QStringList() << "m" << "max",
            QCoreApplication::translate("main", "Maximum value of idx parameter <max>"),
            QCoreApplication::translate("main", "max"),QString("5"));

    QCommandLineOption lowOption(QStringList() << "l" << "from",
            QCoreApplication::translate("main", "Vary from lowest value <val>"),
            QCoreApplication::translate("main", "val"), QString("1"));

    QCommandLineOption highOption(QStringList() << "t" << "to",
            QCoreApplication::translate("main", "Vary to highest value <val>"),
            QCoreApplication::translate("main", "val"), QString("5"));

    QCommandLineOption characteristicSizeOption(QStringList() << "characteristic-size",
            QCoreApplication::translate("main", "Vary to highest value <val>"),
            QCoreApplication::translate("main", "val"), QString("5"));

    QCommandLineOption cantileverBeamOption(QStringList() << "c" << "cantilever", QCoreApplication::translate("main", "Cantilever beam"));
    parser.addOption(idxOption);
    parser.addOption(idxMaxOption);
    parser.addOption(lowOption);
    parser.addOption(highOption);
    parser.addOption(characteristicSizeOption);
    parser.addOption(cantileverBeamOption);
    parser.process(a);

    // prepare params based on the command line options
    BeamParams beamParams;
    icy::ModelPrms modelPrms;


    if(parser.isSet(cantileverBeamOption)) {
        modelPrms.IndentationVelocity = 0.003;
        modelPrms.InitialTimeStep = 0.05;
        beamParams.IndenterSize = 0.25;
        beamParams.beamType = 1;
        beamParams.beamL1 = 5.80;
        beamParams.beamA = 0.70;
        beamParams.beamThickness = 0.64;
        beamParams.CharacteristicLengthMax = 0.2;
        modelPrms.penaltyK = 10000;
        modelPrms.gravity = -0.05;
        modelPrms.SaveVTU = false;
        modelPrms.NewmarkBeta = 0.5;
        modelPrms.NewmarkGamma = 1;
//        modelPrms.Y = 3.7e6;
        if(parser.isSet(characteristicSizeOption))
            beamParams.CharacteristicLengthMax = parser.value(characteristicSizeOption).toDouble();

        if(parser.isSet(idxOption)) {
            int idx = parser.value(idxOption).toInt();
            modelPrms.InstanceNumber = idx;
            double idx_max = parser.value(idxOption).toDouble()-1;
            double low = parser.value(lowOption).toDouble();
            double high = parser.value(highOption).toDouble();
            double idxNormalized = (double)idx/idx_max;
            double val = low+(high-low)*idxNormalized; // interpolated value

            // set volume of the beam equal to val
            double refVolume = beamParams.beamA*beamParams.beamL1*beamParams.beamThickness;
            double normalizedVolume = val/refVolume;
            double coeff = std::cbrt(normalizedVolume);
            beamParams.beamA*=coeff;
            beamParams.beamL1*=coeff;
            beamParams.beamThickness*=coeff;
        }
    }

    // run either using GUI or CLI
    modelPrms.NoGUI = true;
    qDebug() << "no-GUI mode";

    icy::ImplicitModel4 model;
    model.prms = &modelPrms;
    model.beamParams = &beamParams;
    icy::GeneratorTool::GenerateCantileverBeamSetup(&beamParams, &model.mc);
    bool aborted = false;
    do {
        aborted = model.Step();
        cout << "step " << model.cf.StepNumber << "; force " << model.cf.IndenterForce << endl;
    } while(!aborted &&
            (model.cf.StepNumber < model.prms->MaxSteps) &&
            !(model.prms->MaxSolves > 0 && model.cf.TotalSolves > model.prms->MaxSolves) &&
            !(model.prms->StepAfterWhichFractionDetectionIsTriggered > 0 &&
              model.cf.IndenterForce < model.prms->FractionDetectionForceThreshold));

    int reason = 0;
    if(model.prms->StepAfterWhichFractionDetectionIsTriggered > 0 &&
            model.cf.IndenterForce < model.prms->FractionDetectionForceThreshold) reason = 1;

    model.writeCSV(reason);
    return 0;
}

