#include "RodData.hh"

RodData::RodData() : shouldSimulate( true )
{
    rod = NULL; stepper = NULL; rodRenderer = NULL; 
}

RodData::RodData( ElasticRod* i_rod, RodCollisionTimeStepper* i_stepper, RodRenderer* i_rodRenderer ) :
    shouldSimulate( true )
{
    rod = i_rod; stepper = i_stepper; rodRenderer = i_rodRenderer; 
}

RodData::~RodData()
{
    if ( rod != NULL )
        delete rod;
    if ( stepper != NULL )
        delete stepper;
    if ( rodRenderer != NULL )
        delete rodRenderer;
    
    for ( KinematicEdgeDataMap::iterator it = kinematicEdgeDataMap.begin();
              it != kinematicEdgeDataMap.end();
              it++ )
        if ( it->second != NULL )
            delete it->second;
}
         
void RodData::removeKinematicEdge( unsigned int i_edgeNumber )
{ 
    if ( kinematicEdgeDataMap.find( i_edgeNumber ) != kinematicEdgeDataMap.end() )
    {
        delete kinematicEdgeDataMap[ i_edgeNumber ];
        kinematicEdgeDataMap.erase( i_edgeNumber );
    }
}

void RodData::addKinematicEdge( unsigned int i_edgeNumber, ElasticRod* i_rod,
    MaterialFrame* i_materialframe )
{
    removeKinematicEdge( i_edgeNumber );
    KinematicEdgeData* kinematicEdgeData = new KinematicEdgeData( i_edgeNumber, i_rod, i_materialframe );
    kinematicEdgeDataMap[ i_edgeNumber ] = kinematicEdgeData;
}

void RodData::resetKinematicEdge( unsigned int i_edgeNumber, ElasticRod* i_rod, 
    MaterialFrame& i_materialframe )
{
    if ( kinematicEdgeDataMap.find( i_edgeNumber ) != kinematicEdgeDataMap.end() )
    {
        kinematicEdgeDataMap[ i_edgeNumber ]->resetMaterialFrame( i_rod, i_materialframe );
    }
}

void RodData::updateKinematicEdge( unsigned int i_edgeNumber, MaterialFrame& i_materialframe )
{
    if ( kinematicEdgeDataMap.find( i_edgeNumber ) != kinematicEdgeDataMap.end() )
    {
        kinematicEdgeDataMap[ i_edgeNumber ]->updateMaterialFrame( i_materialframe );
    }
}

void RodData::allocateStorage( size_t i_numCVs )
{
    // Make sure we have enough space to store the date for each CV.
    undeformedVertexPositions.resize( i_numCVs );
    initialVertexPositions.resize( i_numCVs );
    prevVertexPositions.resize( i_numCVs );
    currVertexPositions.resize( i_numCVs );
    nextVertexPositions.resize( i_numCVs );

    undeformedMaterialFrame.resize( i_numCVs - 1 );
}

void RodData::resetVertexPositions( vector< Vec3d >& i_vertexPositions )
{
    undeformedVertexPositions = i_vertexPositions;
    initialVertexPositions = i_vertexPositions;
    prevVertexPositions = i_vertexPositions;
    currVertexPositions = i_vertexPositions;
    nextVertexPositions = i_vertexPositions;
}

void RodData::updateNextRodVertexPositions( vector< Vec3d >& i_vertexPositions )
{
    prevVertexPositions = currVertexPositions;

    // Set the current position to be the prev as it will be moved forward in substeps by
    // Beaker as it takes simulation steps.
    //FIXME: Doesn't currVertexPositions already equal nextVertexPositions at the end of a time
    // step?
    currVertexPositions = nextVertexPositions;
    nextVertexPositions = i_vertexPositions;
}

    double KinematicEdgeData::getAngleBetweenReferenceAndStrand()
    { 
          double pi = 3.14159265358979323846264338327950288;
          double theta = 0.0;
          
          if ( rod != NULL )
          {   
                Vec3d m = rod->getMaterial2( edgeNumber );
                m.normalize();
                Vec3d b = materialFrame.m3;
                b.normalize();
    
                Vec3d x = rod->getEdge( edgeNumber );
                x.normalize();
                            
                Vec3d c = m.cross( b );
                double dot = x.dot(c);
                double sign;
                if ( dot == 0 )
                    sign = 0;
                else if ( dot > 0 )
                    sign = 1;
                else {
                    assert( dot < 0 );
                    sign = -1;
                }
                
                theta = atan2( c.norm(), m.dot(b) );
                
                if ( sign < 0 && theta != 0 )
                   theta = 360*(pi/180.0) - theta;
               
          //     theta = acos( rod->getMaterial2( edgeNumber ).dot( materialFrame.m3 ) );
               //theta = acos( rod->getReferenceDirector2( edgeNumber ).dot( materialFrame.m3 ) );
          }
          
          double tt = (theta - 2*pi);
      //    cerr << "tt = " << tt << endl;
          if ( tt < 0.0001 && tt > -0.0001 )
              theta = 0.0;
          
       //   cerr << "Offset is " << theta*(180.0/pi) << endl;
          
          return theta;
    }

    // around that edge.
    double KinematicEdgeData::getAngleFromRodmaterialFrame()
    {
        if ( rod != NULL )
        {
            // Annoyingly the vectors we care about are different axes of the strand frame and
            // the material frame.
            Vec3d rodRefMaterial = rod->getMaterial2( edgeNumber );
            
            double pi = 3.14159265358979323846264338327950288;

            Vec3d x = rod->getEdge( edgeNumber );
            x.normalize();
            Vec3d r = rodRefMaterial;
            r.normalize();
            //Vec3d m = rod->getMaterial2( edgeNumber );
            Vec3d m = materialFrame.m3;
            m.normalize();
            
            Vec3d c = r.cross( m );
            double dot = x.dot(c);
            double sign;
            if ( dot == 0 )
                sign = 0;
            else if ( dot > 0 )
                sign = 1;
            else {
                assert( dot < 0 );
                sign = -1;
            }
            
            double theta = atan2( c.norm(), r.dot(m) );
            
            if ( sign < 0 && theta != 0 ) 
               theta = 360*(pi/180.0) - theta;
          
            double tt = (theta - 2*pi);
         //   cerr << "tt = " << tt << endl;
            if ( tt < 0.0001 && tt > -0.0001 )
              theta = 0.0;
          
            //double theta = acos( rod->getMaterial2( edgeNumber ).dot( materialFrame.m3 ) );
          
            double rTheta = rod->getTheta( edgeNumber );
           // cerr << "rTheta = " << rTheta*(180.0/pi) << endl;
            //cerr << "Theta = " << theta*(180.0/pi) << endl;
            
            //cerr << "total angle = " << (rTheta + ( offsetFromRodRefFrame - theta ))*(180.0/pi) << endl;
            //return rTheta + ( offsetFromRodRefFrame - theta );
            return rTheta +  theta - offsetFromRodRefFrame;
        }
        else
            return 0;
    }
