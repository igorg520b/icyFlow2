#ifndef GENERATORTOOL_H
#define GENERATORTOOL_H

#include "geometry/meshcollection.h"
#include "beamparams.h"
#include <gmsh.h>
#include <vector>
#include <map>

namespace model = gmsh::model;
namespace factory = gmsh::model::occ;

namespace icy {
class GeneratorTool;
typedef std::pair<double,double> p2d;
}

class icy::GeneratorTool
{
public:
    static void GenerateLBeamSetup(BeamParams *beamParams, MeshCollection *output);
    static void GenerateCantileverBeamSetup(BeamParams *beamParams, MeshCollection *output);
private:
    GeneratorTool() {}
    static void GenerateIndenter(BeamParams *beamParams, Mesh *output);
    static void GenerateCantileverIndenter(BeamParams *beamParams, Mesh *output);
    static void GenerateCantileverSupport(BeamParams *beamParams, Mesh *output);
    static void GenerateBeam(BeamParams *beamParams, Mesh *output);
    static void GenerateCantileverBeam(BeamParams *beamParams, Mesh *output);

//    static double cross(p2d v1, p2d v2) { return v1.first*v2.second-v1.second*v2.first; }

    static bool SegmentSegmentIntersection(p2d p, p2d p2, p2d q, p2d q2);
    static bool PointInsideLoop(p2d p, std::vector<p2d> &loop, p2d exteriorPt);

};

#endif // GENERATORTOOL_H
