#ifndef IMPLICITMODEL4_H
#define IMPLICITMODEL4_H

#include <fstream>
#include <vector>
#include "geometry/meshcollection.h"
#include "modelprms.h"
#include "beamparams.h"
#include "frameinfo.h"
#include "linearsystem.h"
#include "bvh/bvht.h"

#include <vtkCellData.h>
#include <vtkPointData.h>
#include <vtkDoubleArray.h>
#include <vtkPointSource.h>
#include <vtkLineSource.h>
#include <vtkOBBTree.h>
#include <vtkPolyDataMapper.h>
#include <vtkMultiBlockDataSet.h>
#include <vtkCompositePolyDataMapper2.h>
#include <vtkPoints.h>
#include <vtkIdList.h>
#include <vtkDataSetSurfaceFilter.h>

namespace icy {
class ImplicitModel4;
}

class icy::ImplicitModel4
{
public:
    MeshCollection mc;
    ModelPrms *prms = nullptr;
    BeamParams *beamParams;
    FrameInfo cf, tcf0;
    std::vector<FrameInfo> allFrames;
    LinearSystem linearSystem;

    bool isReady = false;
    bool kill = false;

    ImplicitModel4();
    void Clear();
    bool Step();    // return true if aborted
    double *extPointsS; // [5][3];
    double *extPointsE;//[5][3];

private:

    void _prepare();    // this is called once
    void _updateStaticStructure();    // this is called once
    void _beginStep();
    void _addCollidingNodesToStructure();
    bool _checkDamage();
    bool _checkDivergence();
    void _XtoDU();

    void _adjustTimeStep();
    void _acceptFrame();
    void _assemble();

    // extensometers
    fstream myfile;
    vtkNew<vtkOBBTree> obbTree;
    vtkNew<vtkDataSetSurfaceFilter> filter1;
    vtkNew<vtkPoints> extPoints;
    vtkNew<vtkIdList> extIdList;
};

#endif // IMPLICITMODEL4_H
