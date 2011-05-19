#ifndef WMFIGRODGROUP_HH_
#define WMFIGRODGROUP_HH_

#include "RodData.hh"

class WmFigRodGroup
{
public:
    WmFigRodGroup();
    ~WmFigRodGroup();

  //  void setRodOptions( int i_rodNumber, RodOptions i_rodOptions );

    void setUndeformedVertexPosition( int i_rodIndex, int i_vertexIndex, BASim::Vec3d& i_newPosition )
    {
        m_rodData[ i_rodIndex ]->setUndeformedVertexPosition( i_vertexIndex, i_newPosition );
    }

    bool simulationNeedsReset()
    {
        return m_simulationNeedsReset;
    }

    void setSimulationNeedsReset( bool i_simulationNeedsReset )
    {
        m_simulationNeedsReset = i_simulationNeedsReset ;
    }

    int addRod();
    int addRod( std::vector<BASim::Vec3d>& i_rodVertices, RodOptions& i_rodOptions, 
                double i_massDamping, BASim::Vec3d& i_gravity,  RodTimeStepper::Method i_solverType, 
                const bool i_isFromCache = false, const bool i_doReverseHairdo = false );

    void addRodsFromCache( std::vector<std::vector< BASim::Vec3d > >& i_rodVertices, RodOptions& i_rodOptions, 
                             double i_massDamping )
    {
        for ( int r=0; r < i_rodVertices.size(); ++r )
        {
            RodOptions rodOptions = i_rodOptions;
            rodOptions.numVertices = (int)i_rodVertices[ r ].size();

            BASim::Vec3d gravity( 0, 0, 0 );
            addRod( i_rodVertices[ r ], rodOptions, i_massDamping, gravity, RodTimeStepper::NONE, true );
        }
    }

    RodOptions& getRodOptions( int i_rodIndex )
    {
        return m_rodData[ i_rodIndex ]->rodOptions;
    }

    BASim::Vec3d& getGravity( int i_rodIndex )
    {
        return m_rodData[ i_rodIndex ]->m_gravity;
    }

    double getMassDamping( int i_rodIndex )
    {
        return m_rodData[ i_rodIndex ]->massDamping();
    }

    void removeAllRods()
    {
        for ( int r = 0; r < m_rodData.size(); ++r )
        {
            delete m_rodData[ r ];
        }

        m_rodData.clear();
    }

    int numberOfRods()
    {
        return (int)m_rodData.size();
    }

    int numberOfRealRods()
    {
        int numberOfRods = 0;
        for( int r = 0; r < m_rodData.size(); ++r )
        {
            if ( !isPlaceHolderRod( r ) )
            {
               ++numberOfRods;
            }
        }
    
        return numberOfRods;
    }

    int numberOfVerticesInRod( int i_rodIndex )
    {
        return (int) m_rodData[ i_rodIndex ]->numberOfVerticesInRod();
    }

    int numberOfEdgesInRod( int i_rodIndex )
    {
        return (int) m_rodData[ i_rodIndex ]->numberOfEdgesInRod();
    }

    bool shouldSimulateRod( int i_rodIndex )
    {
        if ( !m_isReadingFromCache && !m_simulationNeedsReset && !isPlaceHolderRod( i_rodIndex ) && isRodSimulationEnabled(i_rodIndex))
        {
            return true;
        }
        else
        {
            return false;
        }
    }

    void setIsReadingFromCache( bool i_isReadingFromcache )
    {
        m_isReadingFromCache = i_isReadingFromcache;
    }

    bool isPlaceHolderRod( int rodIndex )
    {
        return m_rodData[ rodIndex ]->isPlaceHolderRod();
    }

    bool isRodSimulationEnabled( int rodIndex )
    {
        return m_rodData[ rodIndex ]->isSimulationEnabled();
    }

    void disableSimulationAllRods()
    {
        for( int r = 0; r < m_rodData.size(); ++r )
        {
            m_rodData[ r ]->disableSimulation();
        }
    }

    void enableSimulationAllRods()
    {
        for( int r = 0; r < m_rodData.size(); ++r )
        {
            m_rodData[ r ]->enableSimulation();
        }
    }

    void disableSimulation( int rodIndex )
    {
        m_rodData[ rodIndex ]->disableSimulation();
    }

    void enableSimulation( int rodIndex )
    {
        m_rodData[ rodIndex ]->enableSimulation();
    }


    BASim::Vec3d nextVertexPosition( int i_rodIndex, int i_vertexIndex )
    {
        return m_rodData[ i_rodIndex ]->nextVertexPosition( i_vertexIndex );
    }
    
    void updateRodNextVertexPositions( int i_rodIndex, std::vector<BASim::Vec3d>& inputStrandVertices )
    {
        m_rodData[ i_rodIndex ]->updateNextRodVertexPositions( inputStrandVertices );
    }

    void setUndeformedMaterialFrame( int i_rodIndex, int i_frameIndex, BASim::Vec3d i_v1, BASim::Vec3d i_v2, BASim::Vec3d i_v3 )
    {        
        m_rodData[ i_rodIndex ]->setUndeformedMaterialFrame( i_frameIndex, i_v1, i_v2, i_v3 );
    }

    MaterialFrame undeformedMaterialFrame( int i_rodIndex, int i_frameIndex )
    {
        return m_rodData[ i_rodIndex ]->getUndeformedMaterialFrame( i_frameIndex );
    }

    void updateCurrentVertexPositions( double i_interpolationFactor )
    {
        for( int r = 0; r < m_rodData.size(); ++r )
        {
            if ( !isPlaceHolderRod( r ) )
            {
                int numberOfVertices = m_rodData[ r ]->elasticRod()->nv();

                for ( int v = 0; v < numberOfVertices; ++v )
                {
                    // Calculate the position of the input curve for this rod at the current substep.
                    // This is used to set fixed vertices and also to calculate the hairspray force.
                    m_rodData[ r ]->currVertexPositions[ v ] = ( i_interpolationFactor * m_rodData[ r ]->nextVertexPositions[ v ] +
                                ( 1.0 - i_interpolationFactor ) * m_rodData[ r ]->prevVertexPositions[ v ] );
                }
                
		/*   std::cerr << "rod " << r << "\n";
                
                for ( size_t s=0; s<m_rodData[ r ]->currVertexPositions.size(); ++s )
                {
                    std::cerr << "CURR: " << m_rodData[ r ]->currVertexPositions[ s ] << std::endl;                    
                }
                
                for ( size_t s=0; s<m_rodData[ r ]->prevVertexPositions.size(); ++s )
                {
                    std::cerr << "PREV: " << m_rodData[ r ]->prevVertexPositions[ s ] << std::endl;                    
                }
                
                for ( size_t s=0; s<m_rodData[ r ]->nextVertexPositions.size(); ++s )
                {
                    std::cerr << "NEXT: " << m_rodData[ r ]->nextVertexPositions[ s ] << std::endl;                    
                } 

		*/               
            }
        }
    }

    BASim::Vec3d currentVertexPosition( const int i_rodIndex, const int i_vertexIndex )
    {
        return m_rodData[ i_rodIndex ]->currVertexPositions[ i_vertexIndex ];
    }

    void addKinematicEdge( int i_rodIndex, int i_edgeIndex, MaterialFrame* i_materialframe = NULL )
    {
        m_rodData[ i_rodIndex ]->addKinematicEdge( i_edgeIndex, i_materialframe );
    }

    void removeKinematicEdge( int i_rodIndex, int i_edgeIndex )
    {
        m_rodData[ i_rodIndex ]->removeKinematicEdge( i_edgeIndex );
    }

    void resetKinematicEdge( int i_rodIndex, unsigned int i_edgeNumber, MaterialFrame& i_materialframe )
    {
        m_rodData[ i_rodIndex ]->resetKinematicEdge( i_edgeNumber, i_materialframe );
    }

    void updateKinematicEdge( int i_rodIndex, int i_edgeIndex, MaterialFrame& i_materialframe )
    {
        m_rodData[ i_rodIndex ]->updateKinematicEdge( i_edgeIndex, i_materialframe );
    }

    void setRodParameters( double i_radiusA, double i_radiusB, double i_youngsModulus,
                           double i_shearModulus, double i_viscosity, double i_density )
    {
        for ( int r=0; r< m_rodData.size(); ++r )
        {
            m_rodData[ r ]->setRodParameters( i_radiusA, i_radiusB, i_youngsModulus,
                                              i_shearModulus, i_viscosity, i_density );
        }
    }

    void updateAllBoundaryConditions()
    {
        for ( int r=0; r< m_rodData.size(); ++r )
        {
            if ( shouldSimulateRod( r ) )
            {
                m_rodData[ r ]->updateBoundaryConditions();
            }
        }
    }

    void setNextVertexPosition( int i_rodIndex, int i_vertexIndex, BASim::Vec3d& i_newPosition )
    {
        m_rodData[ i_rodIndex ]->setNextVertexPosition( i_vertexIndex, i_newPosition );
    }

    ElasticRod* elasticRod( int i_rodIndex )
    {
        if(m_rodData.size() == 0)
        {
            return NULL;
        }
        return m_rodData[ i_rodIndex ]->elasticRod();
    }

    /*RodCollisionTimeStepper* collisionStepper( int i_rodIndex )
    {
        return m_rodData[ i_rodIndex ]->collisionStepper();
    }*/

    RodTimeStepper* stepper( int i_rodIndex )
    {
        return m_rodData[ i_rodIndex ]->stepper();
    }


    void setDrawScale( double i_drawScale )
    {
        for ( int r=0; r<m_rodData.size(); ++r )
        {
            m_rodData[ r ]->setDrawScale( i_drawScale );
        }
    }

    void setDrawMode( RodRenderer::DrawMode i_drawMode )
    {
        for ( int r=0; r<m_rodData.size(); ++r )
        {
            m_rodData[ r ]->setDrawMode( i_drawMode );            
        }
    }

    void render()
    {
        for ( int r=0; r<m_rodData.size(); ++r )
        {
            m_rodData[ r ]->render();
        }
    }

    void setColorForSimpleRender( const Color &i_root,
        Color &i_tip )
    {
        for ( int r=0; r<m_rodData.size(); ++r )
        {
            // No rod renderer for fake rods so don't try and change them
            if ( !isPlaceHolderRod( r ) )
            {
                m_rodData[ r ]->rodRenderer()->setColorInSimpleMode( i_root, i_tip);
            }
        }
    }
    
    RodRenderer* rodRenderer( const int i_index )
    {
        return m_rodData[ i_index ]->rodRenderer();
    }
    
    void disableAllRods()
    {
        for ( int r=0; r < m_rodData.size(); ++r )
        {
            m_rodData[ r ]->disableRod();
        }
    }
    
    void addExternalForceToRod( const size_t i_rodIndex, const size_t i_rodVertexIndex, 
                                const BASim::Vec3d i_force )
    {
        m_rodData[ i_rodIndex ]->addExternalForceToVertex( i_rodVertexIndex, i_force );
    }
    
    void resetAllExternalForcesOnRod( const size_t i_rodIndex )
    {
        m_rodData[ i_rodIndex ]->resetExternalForcesOnVertices();        
    }
    
    void resetAllExternalForcesOnRods()
    {
        for ( size_t r = 0; r < m_rodData.size(); ++r )
        {
            m_rodData[ r ]->resetExternalForcesOnVertices();
        }
    }
    
    void setBridsonStepperOnAllRodData( BARodStepper* i_bridsonStepper )
    {
        for ( size_t r = 0; r < m_rodData.size(); ++r )
        {
            m_rodData[ r ]->setBridsonStepper( i_bridsonStepper );
        }
    }

private:
    std::vector< RodData* > m_rodData;

    double m_massDamping;
    bool m_isReadingFromCache;
    bool m_simulationNeedsReset;
};

typedef std::tr1::unordered_map<int, WmFigRodGroup* > RodDataMap;
typedef RodDataMap::iterator RodDataMapIterator;


#endif
