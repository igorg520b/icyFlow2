#ifndef MESH_H
#define MESH_H

#include <QObject>
#include <vector>
#include <algorithm>
#include <vtkNew.h>
#include <vtkUnstructuredGrid.h>
#include <vtkCellType.h>
#include <vtkDataSetMapper.h>
#include <vtkActor.h>
#include <vtkProperty.h>
#include <vtkNamedColors.h>
#include "node.h"
#include "element.h"
#include "cz.h"
#include "face.h"
#include "surfacefragment.h"
#include "translation.h"


namespace icy { class Mesh; }

class icy::Mesh : public QObject
{
    Q_OBJECT
    Q_PROPERTY(bool isDeformable MEMBER isDeformable NOTIFY propertyChanged)
    Q_PROPERTY(bool isIndenter MEMBER isIndenter NOTIFY propertyChanged)
    Q_PROPERTY(int nNodes READ getNumberOfNodes)
    Q_PROPERTY(int nElems READ getNumberOfElems)
    Q_PROPERTY(int nCZs READ getNumberOfCZs)
    Q_PROPERTY(int nFaces READ getNumberOfFaces)

public:
    std::vector<Node> nodes;
    std::vector<Element> elems;
    std::vector<CZ> czs;
    std::vector<Face> faces;

    // the is is a subset of elems
    std::vector<Element*> surfaceElements; // elements that can potentially come in contact

    //  public TranslationCollection translationCollection = new TranslationCollection();

    bool isDeformable = false;
    bool isIndenter = false;
    // bounding box
    double xmin, xmax, ymin, ymax, zmin, zmax;
    vtkNew<vtkUnstructuredGrid> ugrid;
    vtkNew<vtkDataSetMapper> dataSetMapper;
    vtkNew<vtkActor> ugridActor;
    vtkNew<vtkPoints> points;
    vtkNew<vtkNamedColors> colors;

    Mesh();
    ~Mesh();
    void ComputeBoundingBox();
    void AnchorSides();
    void IdentifySurfaceElements();
    void ConnectFaces();
    void CenterSample(double &dx, double &dy);
    void Translate(double dx, double dy, double dz);

    void CreateUGrid();
    void UpdateUGrid();

private:
    int getNumberOfNodes() {return (int)nodes.size();}
    int getNumberOfElems() {return (int)elems.size();}
    int getNumberOfCZs() {return (int)czs.size();}
    int getNumberOfFaces() {return (int)faces.size();}

signals:
    void propertyChanged();
};

#endif // MESH_H
