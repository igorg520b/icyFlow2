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
}

class icy::GeneratorTool
{
public:
    static void GenerateLBeamSetup(BeamParams *beamParams, MeshCollection *output);
private:
    GeneratorTool() {}
    static void GenerateIndenter(BeamParams *beamParams, Mesh *output);
    static void GenerateBeam(BeamParams *beamParams, Mesh *output);
    static void GenerateTest(BeamParams *beamParams, Mesh *output);

};

#endif // GENERATORTOOL_H
