#include "mainwindow.h"
#include "./ui_mainwindow.h"


MainWindow::MainWindow(QWidget *parent)
    : QMainWindow(parent)
    , ui(new Ui::MainWindow)
{
    ui->setupUi(this);

    beamParams = new BeamParams;

    sp=new QSplitter(Qt::Orientation::Horizontal);
    qtw = new QVTKOpenGLNativeWidget();

    vtkNew<vtkNamedColors> colors;

    vtkNew<vtkGenericOpenGLRenderWindow> renderWindow;
    qtw->SetRenderWindow(renderWindow);
    vtkNew<vtkSphereSource> sphereSource;
    vtkNew<vtkPolyDataMapper> sphereMapper;
    sphereMapper->SetInputConnection(sphereSource->GetOutputPort());
    vtkNew<vtkActor> sphereActor;
    sphereActor->SetMapper(sphereMapper);
    sphereActor->GetProperty()->SetColor(colors->GetColor4d("Tomato").GetData());

    vtkNew<vtkRenderer> renderer;
    renderer->AddActor(sphereActor);
    renderer->SetBackground(colors->GetColor3d("SteelBlue").GetData());

    qtw->GetRenderWindow()->AddRenderer(renderer);
    qtw->GetRenderWindow()->SetWindowName("RenderWindowNoUIFile");


    pb = new ObjectPropertyBrowser(nullptr);
    pb->setActiveObject(beamParams);

    sp->addWidget(pb);
    sp->addWidget(qtw);
    setCentralWidget(sp);

}

MainWindow::~MainWindow()
{
    delete ui;
}


void MainWindow::on_actionGenerate_triggered()
{
    // generate simulation setup

    /*
    gmsh::initialize(argc, argv);
    gmsh::option::setNumber("General.Terminal", 1);

    model::add("indenter1");
    double CharacteristicLengthIndenter = 0.02;
    double a = 0.4; // beamA
    double b = 0.4; // beamB
    double l1 = 0.83;
    double l2 = 1.3;
    double c = 0.1; // beam gap
    double d = 0.35; // beam margin
    double h = 0.5; // thickness
    double indsize = 0.15; // indentor size

    double dy = l1 + c + d;
    double dx = d;
    double yoffset = -indsize/2;

    int point3 = factory::addPoint(dx + l2, yoffset, 0+2*h, 1.0);
    int point4 = factory::addPoint(dx + l2-indsize, yoffset, 0+2*h, 1.0);
    int point5 = factory::addPoint(dx + l2-indsize/2, yoffset, h+2*h, 1.0);
    int point6 = factory::addPoint(dx + l2-indsize, yoffset, c/2+2*h, 1.0);
    int point7 = factory::addPoint(dx + l2, yoffset, c/2+2*h, 1.0);

    int circle1 = factory::addCircleArc(point3, point5, point4);
    int line2 = factory::addLine(point4, point6);
    int line3 = factory::addLine(point6, point7);
    int line4 = factory::addLine(point7, point3);

    std::vector<int> curveTags;
    curveTags.push_back(circle1);
    curveTags.push_back(line2);
    curveTags.push_back(line3);
    curveTags.push_back(line4);
    int loopTag = factory::addCurveLoop(curveTags);

    std::vector<int> loops;
    loops.push_back(loopTag);
    int surfaceTag = factory::addPlaneSurface(loops);

    factory::synchronize();

    gmsh::vectorpair vp;
    gmsh::vectorpair vpOut;

    model::getEntities(vp, 2);
    factory::extrude(vp, 0, indsize, 0, vpOut);

    factory::synchronize();
    gmsh::option::setNumber("Mesh.CharacteristicLengthMax", CharacteristicLengthIndenter);
    model::mesh::generate(3);

    gmsh::write("indenter.msh");


    gmsh::finalize();

*/
}

