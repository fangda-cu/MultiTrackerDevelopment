#ifndef RODDATA_HH_
#define RODDATA_HH_

#ifdef WETA
#include <weta/Wfigaro/Physics/ElasticRods/ElasticRod.hh>
#include <weta/Wfigaro/Physics/ElasticRods/AnisotropicRod.hh>
#include <weta/Wfigaro/Physics/ElasticRods/RodCollisionTimeStepper.hh>
#include <weta/Wfigaro/Physics/ElasticRods/RodUtils.hh>
#include <weta/Wfigaro/Render/RodRenderer.hh>
#else
#include <BASim/src/Physics/ElasticRods/ElasticRod.hh>
#include <BASim/src/Physics/ElasticRods/AnisotropicRod.hh>
#include <BASim/src/Physics/ElasticRods/RodCollisionTimeStepper.hh>
#include <BASim/src/Physics/ElasticRods/RodUtils.hh>
#include <BASim/src/Render/RodRenderer.hh>
#endif

#include "SplineAttrEval.hh"
#include <tr1/unordered_map>

using namespace BASim;

class MaterialFrame
{
public:
    Vec3d m1;
    Vec3d m2;
    Vec3d m3;
};


class EdgeTransform
{
public:
    MaterialFrame materialFrame;
    Vec3d position;
};


// FIXME:
// This class has too much in it. The functionality should be moved up to RodData and just the data
// kept in here.

class KinematicEdgeData
{
public:
    KinematicEdgeData( unsigned int i_edgeNumber, ElasticRod* i_rod = NULL, MaterialFrame* i_materialframe = NULL ) 
    { 
        if ( i_materialframe != NULL )
        {
            materialFrame = *i_materialframe;
            rootFrameDefined = true;
        }
        else
            rootFrameDefined = false;
            
        edgeNumber = i_edgeNumber;
        rod = i_rod;
        
        offsetFromRodRefFrame = getAngleBetweenReferenceAndStrand();
    }
    
    // Since we were probably initislised before the rods were create we need to let the user
    // add the rod ptr here.
    void resetMaterialFrame( ElasticRod* i_rod, MaterialFrame& i_materialframe )
    {
        rod = i_rod;
        materialFrame = i_materialframe;
        offsetFromRodRefFrame = 0;
        offsetFromRodRefFrame = getAngleBetweenReferenceAndStrand();
    }

    void updateMaterialFrame( MaterialFrame& i_materialframe )
    {
        materialFrame = i_materialframe;
    }
    // As the rod material frame and input material frame both require one of the vectors
    // to be tangent to the edge then the other two vectors are offset only by some rotation
    // around that edge.
    double getAngleFromRodmaterialFrame();
    double getAngleBetweenReferenceAndStrand();
       
    ElasticRod *rod;
    bool rootFrameDefined;
    double offsetFromRodRefFrame;
    unsigned int edgeNumber; // This defines which edge on the rod the frame refers to    
    MaterialFrame materialFrame;
};

typedef std::tr1::unordered_map< unsigned int, KinematicEdgeData* > KinematicEdgeDataMap;
    

// This class holds all the info from Maya needed to simulate and render a single rod.
class RodData
{
public:
    RodData();
    RodData( ElasticRod* i_rod, RodCollisionTimeStepper* i_stepper, RodRenderer* i_rodRenderer );
    ~RodData();
    
    void removeKinematicEdge( unsigned int i_edgeNumber );
    void addKinematicEdge( unsigned int i_edgeNumber, ElasticRod* i_rod = NULL, MaterialFrame* i_materialframe = NULL );
    void resetKinematicEdge( unsigned int i_edgeNumber, ElasticRod* i_rod, MaterialFrame& i_materialframe );
    void updateKinematicEdge( unsigned int i_edgeNumber, MaterialFrame& i_materialframe );
    void allocateStorage( size_t i_numCVs );
    void resetVertexPositions( vector< Vec3d >& i_vertexPositions );
    
    // If for some reason this rod shouldn't be simulated then set this flag to false. This
    // usually happens when the user has set the rod node to playback from cache.
    bool shouldSimulate;
    
    // Should we keep the ObjectHandle returned by World rather than the actual rod?
    ElasticRod* rod;
    
    RodCollisionTimeStepper* stepper;
    RodRenderer* rodRenderer;
    
    // These variables are updated directly by the WmBunsenRodNode each frame.
    std::vector<Vec3d> undeformedVertexPositions;
    std::vector<Vec3d> initialVertexPositions;
    RodOptions rodOptions;

    // These variables are for interpolating substeps between frames.
    // The position of the input curve vertices at the last frame
    std::vector<Vec3d> prevVertexPositions;
    // The position of the input curve vertices at the current substep
    std::vector<Vec3d> currVertexPositions;
    // The position of the input curve vertices at the next frame
    std::vector<Vec3d> nextVertexPositions;
    // whether this vertex is locked in place or not
   // vector<Bool> vertexFixed;

    double hairSprayScaleFactor;
    cvDataMap forceWeightMap;
    double massDamping;

    // debug info
    std::vector<Vec3d> ALLprevVertexPositions;
    std::vector<Vec3d> ALLnextVertexPositions;
    std::vector<Vec3d> ALLcurrVertexPositions;
    
    // The undeformed material frames for this rod at startTime
    std::vector<MaterialFrame> undeformedMaterialFrame;
    
    // A list of frames for this rod. Each entry corresponds to an edge and a frame for that
    // edge. This is used for the first edge when Barbershop is driving the rods and for
    // other edges when the user is controlling the rod edge by a Maya transform.
    KinematicEdgeDataMap kinematicEdgeDataMap;
};

#endif
