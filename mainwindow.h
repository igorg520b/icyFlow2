#ifndef MAINWINDOW_H
#define MAINWINDOW_H

#include <QMainWindow>
#include <QHBoxLayout>
#include <QSizePolicy>
#include <QPushButton>
#include <QSplitter>
#include <QLabel>
#include <QVBoxLayout>

#include <QtCharts>
#include <QLineSeries>
#include <QChart>
#include <QChartView>

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
#include <vtkPointSource.h>
#include <vtkLineSource.h>
#include <vtkOBBTree.h>
#include <vtkPolyDataMapper.h>
#include <vtkMultiBlockDataSet.h>
#include <vtkCompositePolyDataMapper2.h>
#include <vtkPoints.h>
#include <vtkIdList.h>
#include <vtkDataSetSurfaceFilter.h>

#include <QVTKOpenGLNativeWidget.h>
#include <QVTKOpenGLWidget.h>

#include "objectpropertybrowser.h"
#include <QDebug>
#include <gmsh.h>
#include "beamparams.h"
#include "simulation/implicitmodel4.h"
#include "backgroundworker.h"

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

    enum ViewMode { Elements, AllCZs, DamagedCZs };
    ViewMode viewMode = ViewMode::Elements;
    BeamParams *beamParams = nullptr;
    icy::ImplicitModel4 model;

public slots:
    void updateGUI(bool aborted);

private slots:

    void on_actionScreenshot_triggered();

    void on_actionGenerator_Tool_triggered();

    void on_actionPause_Resume_triggered();

    void on_actionStep_triggered();

    void on_actionElements_triggered();

    void on_actionAllCZs_triggered();

private:
    Ui::MainWindow *ui;
    QVTKOpenGLNativeWidget *qtw;
    QSplitter *sp;
    QVBoxLayout *layout;
    ObjectPropertyBrowser *pb;
    QLabel *statusPausedOrRunning;
    QLabel *statusFrameNumber;
    QLabel *tcf;
    QLabel *statusDamagedCZs;
    QLabel *statusTotalSolves;

    // chart
    QLineSeries *series;
    QChart *chart;
    QChartView *chartView;

    // extensometer chart
    QLineSeries *series2[5];
    QChart *chart2;
    QChartView *chartView2;

    // VTK objects
    vtkNew<vtkGenericOpenGLRenderWindow> renderWindow;
    vtkNew<vtkRenderer> renderer;
    vtkNew<vtkNamedColors> colors;
    vtkNew<vtkXMLUnstructuredGridWriter> writer;
    vtkNew<vtkWindowToImageFilter> windowToImageFilter;
    vtkNew<vtkPNGWriter> writerPNG;

    double extPointsS[5][3];
    double extPointsE[5][3];
    vtkNew<vtkLineSource> extLines[5];
    vtkNew<vtkMultiBlockDataSet> mbds;
    vtkNew<vtkCompositePolyDataMapper2> extMapper;
    vtkNew<vtkActor> extActor;


    // model
    icy::BackgroundWorker* worker;

};
#endif // MAINWINDOW_H
