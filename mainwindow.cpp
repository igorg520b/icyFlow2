#include "mainwindow.h"
#include "./ui_mainwindow.h"


MainWindow::MainWindow(QWidget *parent)
    : QMainWindow(parent)
    , ui(new Ui::MainWindow)
{
    ui->setupUi(this);

    gmsh::initialize();

    beamParams = new BeamParams;

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
    double CharacteristicLengthIndenter = beamParams->CharacteristicLengthIndenter;
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


//    vtkSmartPointer<vtkUnstructuredGrid> ugrid =
//            vtkSmartPointer<vtkUnstructuredGrid>::New();
    ugrid->Reset();

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
//    gmsh::finalize();

}


void MainWindow::on_actionGenerate_Beam_triggered()
{
    double CharacteristicLengthIndenter = beamParams->CharacteristicLengthIndenter;
    double CharacteristicLengthMax = beamParams->CharacteristicLengthMax;
    double rm = beamParams->RefinementMultiplier;
    double a = beamParams->beamA; // beamA
    double b = beamParams->beamB; // beamB
    double l1 = beamParams->beamL1;
    double l2 = beamParams->beamL2;
    double c = beamParams->beamGap; // beam gap
    double d = beamParams->beamMargin; // beam margin
    double h = beamParams->beamThickness; // thickness
    double indsize = beamParams->IndenterSize; // indentor size

//    gmsh::clear();
    gmsh::option::setNumber("General.Terminal", 1);
    model::add("beam1");

    double dy = l1 + c + d;
    double dx = c+d;

    int point1 = factory::addPoint(dx + 0, dy + 0, 0, rm);
    int point2 = factory::addPoint(dx + l2, dy + 0, 0, rm);
    int point22 = factory::addPoint(dx + l2 / 2,dy + 0, 0, rm);
    int point3 = factory::addPoint(dx + l2, dy - a, 0, rm);
    int point4 = factory::addPoint(dx + 0,dy - l1 + c / 2, 0, rm);
    int point5 = factory::addPoint(dx - c,dy - l1 + c / 2, 0, rm);
    int point6 = factory::addPoint(dx - c / 2,dy - l1 + c / 2, 0, rm);
    int point7 = factory::addPoint(dx - c, dy + c, 0, 1.0);
    int point8 = factory::addPoint(dx + l2 + c, dy + c, 0, 1.0);
    int point9 = factory::addPoint(dx + l2 + c, dy - a - c, 0, 1.0);
    int point10 = factory::addPoint(dx + b + 2 * c,dy - a, 0, rm);
    int point11 = factory::addPoint(dx + b,dy - a - 2 * c, 0, rm);
    int point12 = factory::addPoint(dx + b + 2 * c,dy - a - 2 * c,0, rm);
    int point13 = factory::addPoint(dx + b + 2 * c,dy - a - c, 0, rm);
    int point14 = factory::addPoint(dx + b + c,dy - a - 2 * c, 0, rm);
    int point15 = factory::addPoint(dx + b,dy - l1 + c / 2, 0, rm);
    int point16 = factory::addPoint(dx + b + c,dy - l1 + c / 2, 0, rm);
    int point17 = factory::addPoint(dx + b + c / 2,dy - l1 + c / 2, 0, rm);
    int point18 = factory::addPoint(-d, dy + c + d, 0, 1.0);
    int point19 = factory::addPoint(dx + l2 + c + d, dy + c + d, 0, 1.0);
    int point20 = factory::addPoint(dx + l2 + c + d, dy - l1 - c - 2*d, 0, 1.0);
    int point21 = factory::addPoint(-d, -d, 0, rm);
    int point23 = factory::addPoint(dx/4 + l2/4 + c/4 + d/4, -d , 0, rm);
    int point24 = factory::addPoint(-d, dy/4 + c/4 + d/4, 0, rm);
    int point25 = factory::addPoint(-d , dy/2 + c/2 + d/2, 0, 1.0);
    int point26 = factory::addPoint(dx/2 + l2/2 + c/2 + d/2, -d, 0, 1.0);

    int circle1 = factory::addCircleArc(point4, point6, point5);
    int circle8 = factory::addCircleArc(point10, point12, point11);
    int circle11 = factory::addCircleArc(point13, point12, point14);
    int circle14 = factory::addCircleArc(point16, point17, point15);

    int line2 = factory::addLine(point1, point22);
    int line19 = factory::addLine(point22, point2);
    int line3 = factory::addLine(point2, point3);
    int line4 = factory::addLine(point4, point1);
    int line5 = factory::addLine(point5, point7);
    int line6 = factory::addLine(point7, point8);
    int line7 = factory::addLine(point8, point9);
    int line9 = factory::addLine(point3, point10);
    int line10 = factory::addLine(point9, point13);
    int line12 = factory::addLine(point11, point15);
    int line13 = factory::addLine(point14, point16);
    int line15 = factory::addLine(point18, point19);
    int line16 = factory::addLine(point19, point20);
    int line17 = factory::addLine(point20, point26);
    int line18 = factory::addLine(point26, point23);
    int line20 = factory::addLine(point23, point21);
    int line21 = factory::addLine(point21, point24);
    int line22 = factory::addLine(point24, point25);
    int line23 = factory::addLine(point25, point18);


    std::vector<int> curveTags;
    curveTags.push_back(line15);
    curveTags.push_back(line16);
    curveTags.push_back(line17);
    curveTags.push_back(line18);
    curveTags.push_back(line20);
    curveTags.push_back(line21);
    curveTags.push_back(line22);
    curveTags.push_back(line23);
    int loop2 = factory::addCurveLoop(curveTags);

    curveTags.clear();
    curveTags.push_back(line6);
    curveTags.push_back(line7);
    curveTags.push_back(line10);
    curveTags.push_back(circle11);
    curveTags.push_back(line13);
    curveTags.push_back(circle14);
    curveTags.push_back(-line12);
    curveTags.push_back(-circle8);
    curveTags.push_back(-line9);
    curveTags.push_back(-line3);
    curveTags.push_back(-line19);
    curveTags.push_back(-line2);
    curveTags.push_back(-line4);
    curveTags.push_back(circle1);
    curveTags.push_back(line5);
    int loop3 = factory::addCurveLoop(curveTags);

    std::vector<int> loops;
    loops.push_back(loop2);
    loops.push_back(loop3);
    int surfaceTag = factory::addPlaneSurface(loops);

    factory::synchronize();

    gmsh::vectorpair vp;
    gmsh::vectorpair vpOut;

    model::getEntities(vp, 2);
    factory::extrude(vp, 0, 0, h, vpOut);

    factory::synchronize();
    gmsh::option::setNumber("Mesh.CharacteristicLengthMax", CharacteristicLengthMax);
    model::mesh::generate(3);

//    gmsh::write("lbeam.msh");

    std::vector<std::size_t> nodeTags1;
    std::vector<double> nodeCoords, parametricCoords;
    model::mesh::getNodesByElementType(4, nodeTags1, nodeCoords, parametricCoords, -1, false);

    std::vector<std::size_t> elementTags, nodeTags2;

    model::mesh::getElementsByType(4, elementTags, nodeTags2);


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


//    vtkSmartPointer<vtkUnstructuredGrid> ugrid =
//            vtkSmartPointer<vtkUnstructuredGrid>::New();
    ugrid->Reset();
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
