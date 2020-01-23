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
    double CharacteristicLengthIndenter = 0.02;
    double a = beamParams->beamA; // beamA
    double b = beamParams->beamB; // beamB
    double l1 = beamParams->beamL1;
    double l2 = beamParams->beamL2;
    double c = beamParams->beamGap; // beam gap
    double d = beamParams->beamMargin; // beam margin
    double h = beamParams->beamThickness; // thickness
    double indsize = beamParams->IndenterSize; // indentor size

    // generate simulation setup
    qDebug() << "generating indenter";

    gmsh::initialize();
    gmsh::option::setNumber("General.Terminal", 1);

    model::add("indenter1");

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

    //    gmsh::write("indenter.msh");

    std::vector<std::size_t> nodeTags1;
    std::vector<double> nodeCoords, parametricCoords;
    model::mesh::getNodesByElementType(4, nodeTags1, nodeCoords, parametricCoords, -1, false);

    qDebug() << "size node tags " << nodeTags1.size();
    qDebug() << "size node coords " << nodeCoords.size();
    qDebug() << "size parametric coords " << parametricCoords.size();

    std::vector<std::size_t> elementTags, nodeTags2;

    model::mesh::getElementsByType(4, elementTags, nodeTags2);
    qDebug() << "size tags: " << elementTags.size();
    qDebug() << "node tags: " << nodeTags2.size();




    qDebug() << "done generating";
    renderer->RemoveAllViewProps();
    // update qtw

    vtkSmartPointer<vtkPoints> points =
            vtkSmartPointer<vtkPoints>::New();

    points->Allocate(nodeTags1.size());
    QMap<vtkIdType,int>map;

    for(int i = 0; i<nodeTags1.size(); i++) {
        vtkIdType nodeTag = (vtkIdType)nodeTags1[i];
        points->InsertPoint((vtkIdType)i,
                nodeCoords[i*3+0],
                nodeCoords[i*3+1],
                nodeCoords[i*3+2]);
        map.insert(nodeTag,i);
    }


    vtkSmartPointer<vtkUnstructuredGrid> ugrid =
            vtkSmartPointer<vtkUnstructuredGrid>::New();
    ugrid->SetPoints(points);

    ugrid->Allocate(elementTags.size());

    vtkIdType pts2[4];
    for(int i=0;i<elementTags.size();i++)
    {
        pts2[0] = map[nodeTags2[i*4+0]];
        pts2[1] = map[nodeTags2[i*4+1]];
        pts2[2] = map[nodeTags2[i*4+2]];
        pts2[3] = map[nodeTags2[i*4+3]];
        ugrid->InsertNextCell(VTK_TETRA, 4,pts2);
    }




    vtkSmartPointer<vtkDataSetMapper> ugridMapper =
            vtkSmartPointer<vtkDataSetMapper>::New();
    ugridMapper->SetInputData(ugrid);

    vtkSmartPointer<vtkActor> ugridActor =
            vtkSmartPointer<vtkActor>::New();
    ugridActor->SetMapper(ugridMapper);
    ugridActor->GetProperty()->SetColor(colors->GetColor3d("Seashell").GetData());
    ugridActor->GetProperty()->EdgeVisibilityOn();


    renderer->AddActor(ugridActor);

    renderer->ResetCamera();
    renderWindow->Render();

    vtkActorCollection *vac = renderer->GetActors();
    qDebug() << "number of actors: " << vac->GetNumberOfItems();
    gmsh::finalize();

}

