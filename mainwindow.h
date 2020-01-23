#ifndef MAINWINDOW_H
#define MAINWINDOW_H

#include <QMainWindow>
#include <QHBoxLayout>
#include <QSizePolicy>
#include <QPushButton>
#include <QSplitter>

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


#include <QVTKOpenGLNativeWidget.h>
#include <QVTKOpenGLWidget.h>

#include "objectpropertybrowser.h"
#include <QDebug>
#include <gmsh.h>
#include "beamparams.h"

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

private slots:

    void on_actionGenerate_triggered();

private:
    Ui::MainWindow *ui;
    QVTKOpenGLNativeWidget *qtw;
    QSplitter *sp;
    ObjectPropertyBrowser *pb;

    BeamParams *beamParams;

    // VTK objects
    vtkNew<vtkGenericOpenGLRenderWindow> renderWindow;
    vtkNew<vtkRenderer> renderer;
    vtkNew<vtkNamedColors> colors;

};
#endif // MAINWINDOW_H
