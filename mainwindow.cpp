#include "mainwindow.h"
#include "./ui_mainwindow.h"
#include "generatortool.h"

MainWindow::MainWindow(QWidget *parent)
    : QMainWindow(parent)
    , ui(new Ui::MainWindow)
{
    ui->setupUi(this);

    gmsh::initialize();


    sp=new QSplitter(Qt::Orientation::Horizontal);
    qtw = new QVTKOpenGLNativeWidget();

    renderWindow->SetAlphaBitPlanes(1); //?

    qtw->SetRenderWindow(renderWindow);
    vtkNew<vtkSphereSource> sphereSource;
    vtkNew<vtkPolyDataMapper> sphereMapper;
    sphereMapper->SetInputConnection(sphereSource->GetOutputPort());
    vtkNew<vtkActor> sphereActor;
    sphereActor->SetMapper(sphereMapper);
    sphereActor->GetProperty()->SetColor(colors->GetColor4d("Yellow").GetData());
    sphereActor->GetProperty()->EdgeVisibilityOn();

    renderer->AddActor(sphereActor);
    renderer->SetBackground(colors->GetColor3d("LightGrey").GetData());

    qtw->GetRenderWindow()->AddRenderer(renderer);
    qtw->GetRenderWindow()->SetWindowName("RenderWindowNoUIFile");

    pb = new ObjectPropertyBrowser(nullptr);
    pb->setActiveObject(&beamParams);

    sp->addWidget(pb);
    sp->addWidget(qtw);
    setCentralWidget(sp);
}

MainWindow::~MainWindow()
{
    delete ui;
}


void MainWindow::on_actionWrite_VTU_triggered()
{
    writer->SetFileName("beam.vtu");
    writer->SetInputData(ugrid);
    writer->Write();
}

void MainWindow::on_actionScreenshot_triggered()
{
    windowToImageFilter->SetInput(renderWindow);
    windowToImageFilter->SetScale(3);
    writerPNG->SetFileName("screenshot2.png");
    writerPNG->SetInputConnection(windowToImageFilter->GetOutputPort());
    writerPNG->Write();
}

void MainWindow::on_actionGenerator_Tool_triggered()
{
    // user generator tool to set up the scene
    icy::GeneratorTool::GenerateLBeamSetup(&beamParams, &model.mc);
    renderer->RemoveAllViewProps();
    renderer->AddActor(model.mc.mgs[0]->ugridActor);
    renderer->AddActor(model.mc.mgs[1]->ugridActor);

    renderer->ResetCamera();
    renderWindow->Render();
}
