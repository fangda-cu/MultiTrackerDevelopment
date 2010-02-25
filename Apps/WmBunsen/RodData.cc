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
        kinematicEdgeDataMap.erase( 0 );
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
    

    // around that edge.
    double KinematicEdgeData::getAngleFromRodmaterialFrame()
    {
        if ( rod != NULL )
        {
            // Annoyingly the vectors we care about are different axes of the strand frame and
            // the material frame.
            Vec3d rodRefMaterial = rod->getReferenceDirector2( edgeNumber );
            
            double pi = 3.14159265358979323846264338327950288;
            //double angle = acos( rodRefMaterial.dot( materialFrame.m3 ) );
            
            //double angle = atan2(norm(cross(a,b)),dot(a,b));
            Vec3d x =  rod->getEdge( edgeNumber );
            x.normalize();
            Vec3d y = rodRefMaterial;
            y.normalize();
            Vec3d z = materialFrame.m3;
            
            Vec3d c = y.cross( z );
            double dot = x.dot(c);
            double sign;
            if ( dot == 0 )
                sign = 0;
            else if ( dot > 0 )
                sign = 1;
            else if ( dot < 0 )
                sign = -1;
            
            cerr << "sign = " << sign << endl;
            double angle = atan2( c.norm(), y.dot(z) );
            
            if ( sign < 0 )
                angle = 360*(pi/180.0) - angle;
            
            //angle = angle - offsetFromRodRefFrame;
            
            
            cerr << "angle = " << angle << endl;
            
            //if ( angle < 0 )
              //  angle = 180.0*(pi/180.0) - angle;
            
            
        }
        else
            return 0;
    }
