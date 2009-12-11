#ifndef RODDATA_HH_
#define RODDATA_HH_

#include <BASim/BASim>
#include "SplineAttrEval.hh"

using namespace BASim;
using namespace std;

// This class holds all the info needed to simulate and render a single rod.
// It feels like stepper and RodRenderer should perhaps be members of the
// ElasticRod class.
class RodData
{
public:
    RodData();
    RodData( ElasticRod* i_rod, RodCollisionTimeStepper* i_stepper, RodRenderer* i_rodRenderer );
    ~RodData();
    
    
    // If for some reason this rod shouldn't be simulated then set this flag to false. This
    // usually happens when the user has set the rod node to playback from cache.
    bool shouldSimulate;
    
    
    // Should we keep the ObjectHandle returned by World rather than the actual rod?
    ElasticRod* rod;
    RodCollisionTimeStepper* stepper;
    RodRenderer* rodRenderer;
    
    // These variables are updated directly by the WmBunsenRodNode each frame.
    vector<Vec3d> undeformedVertexPositions;
    vector<Vec3d> initialVertexPositions;
    RodOptions rodOptions;

    // These variables are for interpolating substeps between frames.
    // The position of the input curve vertices at the last frame
    vector<Vec3d> prevVertexPositions;
    // The position of the input curve vertices at the current substep
    vector<Vec3d> currVertexPositions;
    // The position of the input curve vertices at the next frame
    vector<Vec3d> nextVertexPositions;

    // whether this vertex is locked in place or not
   // vector<Bool> vertexFixed;

    double hairSprayScaleFactor;
    cvDataMap forceWeightMap;
};

#endif
