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

    chart->createDefaultAxes();
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
    double extOffsets[5][2]={{0.07,0.07},
                             {0.855, 0.075},
                             {1.16, 0.07},
                             {0.93, 0.61},
                             {1.23, 0.70}};
    double extPointsS[5][3];
    double extPointsE[5][3];

    mbds->SetNumberOfBlocks(5);
    for(int i=0;i<5;i++) {
        // set up extensometers
        extLines[i]->SetPoint1(beamParams.beamGap+beamParams.beamMargin+beamParams.beamL2-extOffsets[i][0],
                beamParams.beamGap+beamParams.beamMargin+beamParams.beamL1-extOffsets[i][1],
                beamParams.beamThickness+0.1);
        extLines[i]->SetPoint2(beamParams.beamGap+beamParams.beamMargin+beamParams.beamL2-extOffsets[i][0],
                beamParams.beamGap+beamParams.beamMargin+beamParams.beamL1-extOffsets[i][1],
                beamParams.beamThickness-0.1);

        extPointsS[i][0] = extPointsE[i][0] = beamParams.beamGap+beamParams.beamMargin+beamParams.beamL2-extOffsets[i][0];
        extPointsS[i][1] = extPointsE[i][1] = beamParams.beamGap+beamParams.beamMargin+beamParams.beamL1-extOffsets[i][1];
        extPointsS[i][2] = beamParams.beamThickness+0.1;
        extPointsE[i][2] = beamParams.beamThickness-0.1;

        extLines[i]->Update();
//        extMapper->SetInputConnection(extLines[i]->GetOutputPort());
        mbds->SetBlock(i, extLines[i]->GetOutput());
    }

    extMapper->SetInputDataObject(mbds.GetPointer());
    extActor->SetMapper(extMapper);
    extActor->GetProperty()->SetLineWidth(4);
    extActor->GetProperty()->SetColor(0,0,1);
    renderer->AddActor(extActor);

    // obb tree
//    filter1->SetInputConnection(model.mc.beam->ugrid);
    filter1->SetInputDataObject(model.mc.beam->ugrid);
    filter1->Update();

    obbTree->SetDataSet(filter1->GetOutput());
    obbTree->BuildLocator();

    int result = obbTree->IntersectWithLine(&extPointsS[0][0], &extPointsE[0][0],extPoints, extIdList);
    std::cout << "result " << result << std::endl;
    double pt[3] = {};
    extPoints->GetPoint(0, pt);
    std::cout << pt[0] << ", " << pt[1] << ", " << pt[2] << std::endl;


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
