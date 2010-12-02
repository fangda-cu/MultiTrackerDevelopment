#ifndef RODDATA_HH_
#define RODDATA_HH_

#ifdef WETA
#include <weta/Wfigaro/Physics/ElasticRods/ElasticRod.hh>
#include <weta/Wfigaro/Physics/ElasticRods/AnisotropicRod.hh>
#include <weta/Wfigaro/Physics/ElasticRods/RodCollisionTimeStepper.hh>
#include <weta/Wfigaro/Physics/ElasticRods/AdaptiveBinaryStepper.hh>
#include <weta/Wfigaro/Physics/ElasticRods/RodMassDamping.hh>
#include <weta/Wfigaro/Physics/ElasticRods/RodGravity.hh>
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
    BASim::Vec3d m1;
    BASim::Vec3d m2;
    BASim::Vec3d m3;
};


class EdgeTransform
{
public:
    MaterialFrame materialFrame;
    BASim::Vec3d position;
};


// FIXME:
// This class has too much in it. The functionality should be moved up to RodData and just the data
// kept in here.

class KinematicEdgeData
{
public:
    explicit KinematicEdgeData( unsigned int i_edgeNumber, ElasticRod* i_rod = NULL, MaterialFrame* i_materialframe = NULL ) 
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
    RodData( RodOptions& i_rodOptions, std::vector<BASim::Vec3d>& i_rodVertexPositions,
             double i_massDamping, BASim::Vec3d& i_gravity, RodTimeStepper::Method i_solverType,
             bool i_isReadingFromCache = false, bool i_doReverseHairdo = false );
    //RodData( ElasticRod* i_rod, RodCollisionTimeStepper* i_stepper, RodRenderer* i_rodRenderer );
    ~RodData();
    
    void removeKinematicEdge( unsigned int i_edgeNumber );

    // if no frame passed in then Null means just lock the edge and ignore rotation.
    void addKinematicEdge( unsigned int i_edgeNumber, MaterialFrame* i_materialframe = NULL );
    void resetKinematicEdge( unsigned int i_edgeNumber, MaterialFrame& i_materialframe );
    void updateKinematicEdge( unsigned int i_edgeNumber, MaterialFrame& i_materialframe );
    void allocateStorage( int i_numCVs );
    void resetVertexPositions( vector< BASim::Vec3d >& i_vertexPositions );
    void updateNextRodVertexPositions( vector< BASim::Vec3d >& i_vertexPositions );

    void updateBoundaryConditions();

    bool isPlaceHolderRod()
    {
        return m_isPlaceHolderRod;
    }

    int numberOfVerticesInRod()
    {
        if ( m_rod != NULL )
        {
            return m_rod->nv();
        }
        else
        {
            return 0;
        }
    }
	
	int numberOfEdgesInRod()
    {
        if ( m_rod != NULL )
        {
            return m_rod->ne();
        }
        else
        {
            return 0;
        }
    }
    
    void disableRod()
    {
        if ( m_stepper != NULL )
        {
            m_stepper->setEnabled( false );
        }
    }


    // We store the number of vertices in the fake rod because this is used to skip over the input
    // when coming from Barbershop as the input vertices are one big array so we need to keep
    // track of which verts come from which strand.
    /*void initialiseFakeRod( int i_numberOfVertices )
    {
        m_isFakeRod = true;
        m_numVerticesInFakeRod = i_numberOfVertices;
    }

    int verticesInFakeRod()
    {
        return m_numVerticesInFakeRod;
    }

    bool isFakeRod()
    {
        return m_isFakeRod;
    }

    void setRodComingFromCache( bool i_isRodComingFromCache )
    {
        m_isRodComingFromCache = i_isRodComingFromCache;
    }

    bool isRodComingFromCache()
    {
        return m_isRodComingFromCache;
    }

    bool shouldSimulate()
    {
        if ( m_isFakeRod || m_isRodComingFromCache )
        {
            return false;
        }
        
        return true;
    }

    void setSimulationNeedsReset( bool i_simulationNeedsReset )
    {
        m_simulationNeedsReset = i_simulationNeedsReset;
    }

    /*void setStepperEnabled( bool i_isEnabled )
    {
        if ( !m_isFakeRod )
            stepper->setEnabled( i_isEnabled );
    }

    bool isStepperEnabled()
    {
        if ( m_isFakeRod )
        {
            cerr << "rod is fake\n";
            return false;
        }
        else
        {
            return stepper->isEnabled();
        }
    }*/

    void setRodParameters( double i_radiusA, double i_radiusB, double i_youngsModulus,
                           double i_shearModulus, double i_viscosity, double i_density )
    {
        if ( !m_isPlaceHolderRod )
        {
            m_rod->setRadius( i_radiusA, i_radiusB );
            m_rod->setYoungsModulus( i_youngsModulus );
            m_rod->setShearModulus( i_shearModulus );
            m_rod->setViscosity( i_viscosity );
            m_rod->setDensity( i_density );
        }
    }

    void setDrawScale( double i_drawScale )
    {
        if ( !m_isPlaceHolderRod ) {
            //m_rodRenderer->setDrawScale( i_drawScale );
            m_rod->setRadiusScale( i_drawScale );
        }

    }

    void setDrawMode( RodRenderer::DrawMode i_drawMode )
    {
        if ( !m_isPlaceHolderRod )
            m_rodRenderer->setMode( i_drawMode );
    }

    void render()
    {
        if ( !m_isPlaceHolderRod )
            m_rodRenderer->render();
    }
       
    RodRenderer* rodRenderer()
    {
        return m_rodRenderer;
    }

 /*   static int numberOfRealRods( std::vector< RodData* >* i_pRodData )
    {
        if ( i_pRodData == NULL )
        {
            return 0;
        }

        int numberOfRods = 0;
        for ( int r = 0; r < (int)i_pRodData->size(); ++r )
        {
            if ( !( (*i_pRodData)[ r ]->isFakeRod() ) )
            {
                numberOfRods++;
            }            
        }

        return numberOfRods;
    }*/

    ElasticRod* elasticRod()
    {
        return m_rod;
    }

     RodCollisionTimeStepper* collisionStepper()
    {
        return m_stepper;
    }

    // FIXME: commented out as it does nothing useful.
    /*void setMassDamping( double i_massDamping )
    {
        m_massDamping = i_massDamping;
    }*/

    BASim::Vec3d nextVertexPosition( int i_vertexIndex )
    {
        return nextVertexPositions[ i_vertexIndex ];
    }

    void setUndeformedMaterialFrame( int i_frameIndex, BASim::Vec3d i_v1, BASim::Vec3d i_v2, BASim::Vec3d i_v3 )
    {        
        undeformedMaterialFrame[ i_frameIndex ].m1 = i_v1;
        undeformedMaterialFrame[ i_frameIndex ].m2 = i_v1;
        undeformedMaterialFrame[ i_frameIndex ].m3 = i_v3;
    }

    MaterialFrame getUndeformedMaterialFrame( int i_frameIndex )
    {
        return undeformedMaterialFrame[ i_frameIndex ];
    }


    void setUndeformedVertexPosition( int i_vertexIndex, BASim::Vec3d& i_newPosition )
    {
        undeformedVertexPositions[ i_vertexIndex ] = i_newPosition;
    }

    void setNextVertexPosition( int i_vertexIndex, BASim::Vec3d& i_newPosition )
    {
        nextVertexPositions[ i_vertexIndex ] = i_newPosition;
    }

     double massDamping()
    {
        return m_massDamping;
    }

//private:
    
    // Should we keep the ObjectHandle returned by World rather than the actual rod?
    ElasticRod* m_rod;
    
    RodCollisionTimeStepper* m_stepper;
    RodRenderer* m_rodRenderer;
    
    // These variables are updated directly by the WmBunsenRodNode each frame.
    std::vector<BASim::Vec3d> undeformedVertexPositions;
    std::vector<BASim::Vec3d> initialVertexPositions;
    RodOptions rodOptions;

    // These variables are for interpolating substeps between frames.
    // The position of the input curve vertices at the last frame
    std::vector<BASim::Vec3d> prevVertexPositions;
    // The position of the input curve vertices at the current substep
    std::vector<BASim::Vec3d> currVertexPositions;
    // The position of the input curve vertices at the next frame
    std::vector<BASim::Vec3d> nextVertexPositions;
    // whether this vertex is locked in place or not
   // vector<Bool> vertexFixed;

    double hairSprayScaleFactor;
    cvDataMap forceWeightMap;
    
    // debug info
    //std::vector<BASim::Vec3d> ALLprevVertexPositions;
    //std::vector<BASim::Vec3d> ALLnextVertexPositions;
    //std::vector<BASim::Vec3d> ALLcurrVertexPositions;
    
    // The undeformed material frames for this rod at startTime
    std::vector<MaterialFrame> undeformedMaterialFrame;
    
    // A list of frames for this rod. Each entry corresponds to an edge and a frame for that
    // edge. This is used for the first edge when Barbershop is driving the rods and for
    // other edges when the user is controlling the rod edge by a Maya transform.
    KinematicEdgeDataMap kinematicEdgeDataMap;

    bool m_isFakeRod;

    // As we don't create a rod for fake rods we need to store the number of vertices the rod would
    // have had. So that when updating from input we know how many input vertices to skip past
    // in the big input vertex array.
    int m_numVerticesInFakeRod;

    bool m_isRodComingFromCache;
    bool m_simulationNeedsReset;

    bool m_isPlaceHolderRod;
    double m_massDamping;
    BASim::Vec3d m_gravity;
    
    std::vector< float > m_targetDensityPerEdge;
};

#endif
