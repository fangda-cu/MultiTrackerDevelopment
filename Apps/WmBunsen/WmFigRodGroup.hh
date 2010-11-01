#ifndef WMFIGRODGROUP_HH_
#define WMFIGRODGROUP_HH_

#include "RodData.hh"

class WmFigRodGroup
{
public:
    WmFigRodGroup();
    ~WmFigRodGroup();

  //  void setRodOptions( int i_rodNumber, RodOptions i_rodOptions );

    void setUndeformedVertexPosition( int i_rodIndex, int i_vertexIndex, Vec3d& i_newPosition )
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
    int addRod( std::vector<Vec3d>& i_rodVertices, RodOptions& i_rodOptions, 
                   double i_massDamping, BASim::Vec3d& i_gravity,  RodTimeStepper::Method i_solverType, bool i_isFromCache = false );

    void addRodsFromCache( vector<vector<BASim::Vec3d> >& i_rodVertices, RodOptions& i_rodOptions, 
                             double i_massDamping )
    {
        for ( int r=0; r < i_rodVertices.size(); ++r )
        {
            RodOptions rodOptions = i_rodOptions;
            rodOptions.numVertices = (int)i_rodVertices[ r ].size();

            Vec3d gravity( 0, 0, 0 );
            addRod( i_rodVertices[ r ], rodOptions, i_massDamping, gravity, RodTimeStepper::NONE, true );
        }
    }

    RodOptions& getRodOptions( int i_rodIndex )
    {
        return m_rodData[ i_rodIndex ]->rodOptions;
    }

    Vec3d& getGravity( int i_rodIndex )
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
        if ( !m_isReadingFromCache && !m_simulationNeedsReset && !isPlaceHolderRod( i_rodIndex ) )
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

    Vec3d nextVertexPosition( int i_rodIndex, int i_vertexIndex )
    {
        return m_rodData[ i_rodIndex ]->nextVertexPosition( i_vertexIndex );
    }
    
    void updateRodNextVertexPositions( int i_rodIndex, vector<Vec3d>& inputStrandVertices )
    {
        m_rodData[ i_rodIndex ]->updateNextRodVertexPositions( inputStrandVertices );
    }

    void setUndeformedMaterialFrame( int i_rodIndex, int i_frameIndex, Vec3d i_v1, Vec3d i_v2, Vec3d i_v3 )
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
            }
        }
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

    void setNextVertexPosition( int i_rodIndex, int i_vertexIndex, Vec3d& i_newPosition )
    {
        m_rodData[ i_rodIndex ]->setNextVertexPosition( i_vertexIndex, i_newPosition );
    }

    ElasticRod* elasticRod( int i_rodIndex )
    {
        return m_rodData[ i_rodIndex ]->elasticRod();
    }

     RodCollisionTimeStepper* collisionStepper( int i_rodIndex )
    {
        return m_rodData[ i_rodIndex ]->collisionStepper();
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

    RodRenderer* rodRenderer( const int i_index )
    {
        return m_rodData[ i_index ]->rodRenderer();
    }

private:
    std::vector< RodData* > m_rodData;

    double m_massDamping;
    bool m_isReadingFromCache;
    bool m_simulationNeedsReset;
};

typedef tr1::unordered_map<int, WmFigRodGroup* > RodDataMap;
typedef RodDataMap::iterator RodDataMapIterator;


#endif
