#ifndef RODDATA_HH_
#define RODDATA_HH_

#include <BASim/BASim>

using namespace BASim;
using namespace std;

// This class holds all the info needed to simulate and render a single rod.
// It feels like stepper and RodRenderer should perhaps be members of the
// ElasticRod class.
class RodData
{
public:
    RodData();
    RodData( ElasticRod* i_rod, RodTimeStepper* i_stepper, RodRenderer* i_rodRenderer );
    ~RodData();
    
    ElasticRod* rod;
    RodTimeStepper* stepper;
    RodRenderer* rodRenderer;
    
    // These variables are updated directly by the WmBunsenRodNode each frame.
    vector<Vec3d> undeformedVertexPositions;
    vector<Vec3d> initialVertexPositions;
    RodOptions rodOptions;
};

#endif
