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
    chart = new QChart();
    chart->legend()->hide();
    chart->addSeries(series);

    chart->createDefaultAxes();
    chartView = new QChartView(chart);
    chartView->setRenderHint(QPainter::Antialiasing);

    QVBoxLayout *layout = new QVBoxLayout;
    QWidget *leftContainer = new QWidget;
    leftContainer->setLayout(layout);

    // chart for extensometers
    chart2 = new QChart();
    chart2->legend()->hide();
    for(int i=0;i<5;i++) {
        series2[i] = new QLineSeries();
        series2[i]->append(0, -0.001);
        series2[i]->append(1, 0);
        chart2->addSeries(series2[i]);
    }
    chart2->createDefaultAxes();
    chartView2 = new QChartView(chart2);
    chartView2->setRenderHint(QPainter::Antialiasing);


    statusPausedOrRunning = new QLabel("paused");
    statusFrameNumber = new QLabel("-");
    tcf = new QLabel("-");

    ui->statusbar->addWidget(statusPausedOrRunning);
    ui->statusbar->addWidget(statusFrameNumber);
    ui->statusbar->addWidget(tcf);

    gmsh::initialize();

    sp=new QSplitter(Qt::Orientation::Horizontal);
    qtw = new QVTKOpenGLNativeWidget();

    renderWindow->SetAlphaBitPlanes(1); //?

    qtw->SetRenderWindow(renderWindow);
    renderer->SetBackground(colors->GetColor3d("White").GetData());
    qtw->GetRenderWindow()->AddRenderer(renderer);

    pb = new ObjectPropertyBrowser(nullptr);
    pb->setActiveObject(&model.prms);

    layout->addWidget(pb);
    layout->addWidget(chartView2);
    layout->addWidget(chartView);
    sp->addWidget(leftContainer);
    sp->addWidget(qtw);
    setCentralWidget(sp);

    worker = new icy::BackgroundWorker(&model);
    connect(worker, SIGNAL(stepCompleted(bool)), this, SLOT(updateGUI(bool)));
    model.beamParams = &beamParams;
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

    double extOffsets[5][2]={{0.07,0.07},
                             {0.855, 0.075},
                             {1.16, 0.07},
                             {0.93, 0.61},
                             {1.23, 0.70}};


    mbds->SetNumberOfBlocks(5);
    for(int i=0;i<5;i++) {
        extPointsS[i][0] = extPointsE[i][0] = beamParams.beamGap+beamParams.beamMargin+beamParams.beamL2-extOffsets[i][0];
        extPointsS[i][1] = extPointsE[i][1] = beamParams.beamGap+beamParams.beamMargin+beamParams.beamL1-extOffsets[i][1];
        extPointsS[i][2] = beamParams.beamThickness+0.1;
        extPointsE[i][2] = beamParams.beamThickness-0.1;

        extLines[i]->SetPoint1(&extPointsS[i][0]);
        extLines[i]->SetPoint2(&extPointsE[i][0]);

        extLines[i]->Update();
        mbds->SetBlock(i, extLines[i]->GetOutput());
    }

    extMapper->SetInputDataObject(mbds.GetPointer());
    extActor->SetMapper(extMapper);
    extActor->GetProperty()->SetLineWidth(6);
    extActor->GetProperty()->SetColor(1,1,1);

    renderer->AddActor(extActor);
    model.extPointsE = &extPointsE[0][0];
    model.extPointsS = &extPointsS[0][0];
}

void MainWindow::updateGUI(bool aborted)
{
    if(aborted) statusPausedOrRunning->setText("aborted");
    if(!worker->running) statusPausedOrRunning->setText("paused");
    else statusPausedOrRunning->setText("running");

    statusFrameNumber->setText(QString::number(model.cf.StepNumber));
    tcf->setText(QString::number(model.cf.TimeScaleFactor));
    model.mc.UpdateActors();
    renderWindow->Render();

    // update chart

    for(int i=0;i<5;i++) series2[i]->clear();
    series->clear();
    for(size_t i=0;i<model.allFrames.size();i++) {
        icy::FrameInfo &f = model.allFrames[i];
        series->append(f.SimulationTime,f.IndenterForce);

        for(int i=0;i<5;i++) series2[i]->append(f.SimulationTime, f.extensometerDisplacements[i]);
    }


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
