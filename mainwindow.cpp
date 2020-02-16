#include "mainwindow.h"
#include "./ui_mainwindow.h"
#include "generatortool.h"
#include "simulation/frameinfo.h"

MainWindow::MainWindow(QWidget *parent)
    : QMainWindow(parent)
    , ui(new Ui::MainWindow)
{
    ui->setupUi(this);

    series = new QLineSeries();
    series->append(0, 0);
    series->append(1, 5000);
    QChart *chart = new QChart();
    chart->legend()->hide();
    chart->addSeries(series);

/*
    QValueAxis *axisX = new QValueAxis();
    axisX->setTitleText("time");
    axisX->setLabelFormat("%g");
    axisX->setTickCount(100);
    chart->addAxis(axisX, Qt::AlignBottom);
    series->attachAxis(axisX);

    QValueAxis *axisY = new QValueAxis();
    axisX->setTitleText("force");
    axisX->setLabelFormat("%g");
    axisX->setTickCount(100);
    chart->addAxis(axisX, Qt::AlignBottom);
    series->attachAxis(axisY);
*/
    chart->createDefaultAxes();
//    chart->setTitle("Simple line chart example");
    QChartView *chartView = new QChartView(chart);
    chartView->setRenderHint(QPainter::Antialiasing);

    QVBoxLayout *layout = new QVBoxLayout;
    QWidget *leftContainer = new QWidget;
    leftContainer->setLayout(layout);


    statusPausedOrRunning = new QLabel("paused");
    statusFrameNumber = new QLabel("-");

    ui->statusbar->addWidget(statusPausedOrRunning);
    ui->statusbar->addWidget(statusFrameNumber);

    gmsh::initialize();

    sp=new QSplitter(Qt::Orientation::Horizontal);
    qtw = new QVTKOpenGLNativeWidget();

    renderWindow->SetAlphaBitPlanes(1); //?

    qtw->SetRenderWindow(renderWindow);
    renderer->SetBackground(colors->GetColor3d("LightGrey").GetData());
    qtw->GetRenderWindow()->AddRenderer(renderer);

    pb = new ObjectPropertyBrowser(nullptr);
    pb->setActiveObject(&model.prms);

    layout->addWidget(pb);
    layout->addWidget(chartView);
    sp->addWidget(leftContainer);
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

    // update chart

    series->clear();
    for(size_t i=0;i<model.allFrames.size();i++) {
        icy::FrameInfo &f = model.allFrames[i];
        series->append(f.SimulationTime,f.IndenterForce);
    }

//    chartView->update();

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

void MainWindow::on_actionElements_triggered()
{
    viewMode = ViewMode::Elements;
    renderer->RemoveAllViewProps();
    renderer->AddActor(model.mc.mgs[0]->ugridActor);
    renderer->AddActor(model.mc.mgs[1]->ugridActor);
    renderer->ResetCamera();
    renderWindow->Render();
}

void MainWindow::on_actionAllCZs_triggered()
{
    viewMode = ViewMode::AllCZs;
    renderer->RemoveAllViewProps();
    renderer->AddActor(model.mc.beam->ugridActor_czs);
    renderer->ResetCamera();
    renderWindow->Render();
}
