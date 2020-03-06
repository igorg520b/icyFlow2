#include "mainwindow.h"

#include <QApplication>
#include <QSurfaceFormat>
#include <QCommandLineParser>
#include <QDebug>

int main(int argc, char *argv[])
{
    QApplication a(argc, argv);
    QApplication::setApplicationName("icFlow");
    QApplication::setApplicationVersion("2.1");

    QSurfaceFormat fmt = QVTKOpenGLNativeWidget::defaultFormat();
    fmt.setAlphaBufferSize(0);
    QSurfaceFormat::setDefaultFormat(fmt);
    MainWindow w;

    // parse command line options
    QCommandLineParser parser;
    parser.addHelpOption();
    parser.addVersionOption();
    // An option with a value
    QCommandLineOption fpbOption(QStringList() << "b" << "four-point-bending",
            QCoreApplication::translate("main", "Model four-point beam bending with index parameter <idx>"),
            QCoreApplication::translate("main", "idx"));
    parser.addOption(fpbOption);
    parser.process(a);

    if(parser.isSet(fpbOption)) {
        qDebug() << "no-GUI mode";
        QString fpbIndexValue = parser.value(fpbOption);
        qDebug() << fpbIndexValue;
    } else {
        w.showMaximized();
    }
    return a.exec();
}
