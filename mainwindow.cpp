#include "mainwindow.h"
#include "./ui_mainwindow.h"
#include "generatortool.h"

MainWindow::MainWindow(QWidget *parent)
    : QMainWindow(parent)
    , ui(new Ui::MainWindow)
{
    ui->setupUi(this);

    statusPausedOrRunning = new QLabel("paused");
    statusFrameNumber = new QLabel("-");

    ui->statusbar->addWidget(statusPausedOrRunning);
    ui->statusbar->addWidget(statusFrameNumber);

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

    worker = new icy::BackgroundWorker(&model);
    connect(worker, SIGNAL(stepCompleted(bool)), this, SLOT(updateGUI(bool)));
}

MainWindow::~MainWindow()
{
    delete worker;
    delete ui;
}


void MainWindow::showEvent( QShowEvent*)
{
    // for testing
    ui->actionGenerator_Tool->trigger();
}

void MainWindow::updateGUI(bool aborted)
{
    if(aborted) statusPausedOrRunning->setText("aborted");
    if(!worker->running) statusPausedOrRunning->setText("paused");
    else statusPausedOrRunning->setText("running");

    statusFrameNumber->setText(QString::number(model.cf.StepNumber));
    model.mc.UpdateActors();
    renderWindow->Render();
}


void MainWindow::on_actionWrite_VTU_triggered()
{
    writer->SetFileName("beam.vtu");
    writer->SetInputData(model.mc.beam->ugrid);
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

void MainWindow::on_actionPause_Resume_triggered()
{
    worker->toggle();
    if(worker->timeToPause)statusPausedOrRunning->setText("stopping");
    else statusPausedOrRunning->setText("resuming");
}

void MainWindow::on_actionStep_triggered()
{
    model.Step();
    updateGUI(false);
}
