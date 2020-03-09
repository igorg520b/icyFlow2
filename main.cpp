#include "mainwindow.h"

#include <QApplication>
#include <QSurfaceFormat>
#include <QCommandLineParser>
#include <QDebug>
#include "simulation/modelprms.h"
#include <cmath>

int main(int argc, char *argv[])
{
    QApplication a(argc, argv);
    QApplication::setApplicationName("icFlow");
    QApplication::setApplicationVersion("2.1");


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

    QCommandLineOption highOption(QStringList() << "h" << "to",
            QCoreApplication::translate("main", "Vary to highest value <val>"),
            QCoreApplication::translate("main", "val"), QString("5"));

    QCommandLineOption cantileverBeamOption(QStringList() << "c" << "cantilever", QCoreApplication::translate("main", "Cantilever beam"));
    QCommandLineOption noGUIOption(QStringList() << "n" << "no-gui", QCoreApplication::translate("main", "Run without GUI"));
    parser.addOption(idxOption);
    parser.addOption(idxMaxOption);
    parser.addOption(lowOption);
    parser.addOption(highOption);
    parser.addOption(cantileverBeamOption);
    parser.addOption(noGUIOption);
    parser.process(a);

    // prepare params based on the command line options
    BeamParams beamParams;
    icy::ModelPrms modelPrms;


    if(parser.isSet(cantileverBeamOption)) {
        beamParams.beamType = 1;
        if(parser.isSet(idxOption)) {
            double idx = parser.value(idxOption).toDouble();
            double idx_max = parser.value(idxOption).toDouble()-1;
            double low = parser.value(lowOption).toInt();
            double high = parser.value(highOption).toInt();
            double idxNormalized = idx/idx_max;
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
    if(parser.isSet(noGUIOption)) {
        qDebug() << "no-GUI mode";
        return 0;
    } else {
        QSurfaceFormat fmt = QVTKOpenGLNativeWidget::defaultFormat();
        fmt.setAlphaBufferSize(0);
        QSurfaceFormat::setDefaultFormat(fmt);
        MainWindow w;
        w.beamParams = &beamParams;
        w.model.prms = &modelPrms;
        w.showMaximized();
        return a.exec();
    }
}
