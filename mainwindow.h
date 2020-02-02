#ifndef MAINWINDOW_H
#define MAINWINDOW_H

#include <QMainWindow>
#include <QHBoxLayout>
#include <QSizePolicy>
#include <QPushButton>
#include <QSplitter>
#include <QLabel>

#include <vtkGenericOpenGLRenderWindow.h>
#include <vtkNamedColors.h>
#include <vtkNew.h>
#include <vtkPolyDataMapper.h>
#include <vtkProperty.h>
#include <vtkRenderWindow.h>
#include <vtkRenderer.h>
#include <vtkSphereSource.h>
#include <vtkVersion.h>
#include <vtkUnstructuredGrid.h>
#include <vtkDataSetMapper.h>
#include <vtkCamera.h>
#include <vtkXMLUnstructuredGridWriter.h>
#include <vtkWindowToImageFilter.h>
#include <vtkPNGWriter.h>

#include <QVTKOpenGLNativeWidget.h>
#include <QVTKOpenGLWidget.h>

#include "objectpropertybrowser.h"
#include <QDebug>
#include <gmsh.h>
#include "beamparams.h"
#include "simulation/implicitmodel4.h"
#include "backgroundworker.h"

namespace model = gmsh::model;
namespace factory = gmsh::model::occ;

QT_BEGIN_NAMESPACE
namespace Ui { class MainWindow; }
QT_END_NAMESPACE

class MainWindow : public QMainWindow
{
    Q_OBJECT

public:
    MainWindow(QWidget *parent = nullptr);
    ~MainWindow();
    void showEvent( QShowEvent* event ) override;

public slots:
    void updateGUI(bool aborted);

private slots:

    void on_actionWrite_VTU_triggered();

    void on_actionScreenshot_triggered();

    void on_actionGenerator_Tool_triggered();

    void on_actionPause_Resume_triggered();

private:
    Ui::MainWindow *ui;
    QVTKOpenGLNativeWidget *qtw;
    QSplitter *sp;
    ObjectPropertyBrowser *pb;
    QLabel *statusPausedOrRunning;
    QLabel *statusFrameNumber;


    // VTK objects
    vtkNew<vtkGenericOpenGLRenderWindow> renderWindow;
    vtkNew<vtkRenderer> renderer;
    vtkNew<vtkNamedColors> colors;
    vtkNew<vtkXMLUnstructuredGridWriter> writer;
    vtkNew<vtkWindowToImageFilter> windowToImageFilter;
    vtkNew<vtkPNGWriter> writerPNG;

    icy::ImplicitModel4 model;
    BeamParams beamParams;
    icy::BackgroundWorker* worker;

};
#endif // MAINWINDOW_H
