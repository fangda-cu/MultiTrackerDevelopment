#include "RodData.hh"

RodData::RodData() : m_rod( NULL), m_stepper( NULL ), m_rodRenderer( NULL ), m_gravity( 0.0, 0.0, 0.0 )
{
    m_isPlaceHolderRod = true;
}

RodData::RodData( RodOptions& i_rodOptions, std::vector<BASim::Vec3d>& i_rodVertexPositions,
                  double i_massDamping, BASim::Vec3d& i_gravity, RodTimeStepper::Method i_solverType, 
                  bool i_isReadingFromCache, bool i_doReverseHairdo ) : 
m_rod( NULL), m_stepper( NULL ), m_rodRenderer( NULL ), m_massDamping( i_massDamping )
{
    m_rod = setupRod( i_rodOptions,
                      i_rodVertexPositions,
                      i_rodVertexPositions );

    m_rodRenderer = new RodRenderer( *m_rod );

    // If the rod is coming from the cache file then we don't need the stepper or forces.
    // FIXME: The above setup rod code does not need to be called either really. We need
    // a fake setupRod function that just gives space for the vertices.
    if ( !i_isReadingFromCache )
    {
      // If Adaptive Time Stepping is used,
	    if(1) {
		    RodTimeStepper* stepper = new RodTimeStepper( *m_rod );
		    stepper->setDiffEqSolver( i_solverType );
	    
		    // These get deleted by RodTimeStepper
		    stepper->addExternalForce( new RodMassDamping( m_massDamping ) );
		    
		    if ( i_gravity.norm() > 0)
		    {
		      stepper->addExternalForce( new RodGravity( i_gravity ) );
		      m_gravity = i_gravity;
		    }
    
        
		    AdaptiveBinaryStepper* adpstep = new AdaptiveBinaryStepper( m_rod, stepper );
		    //m_steppers.push_back(adpstep);
	    
		    m_stepper = new RodCollisionTimeStepper( adpstep, m_rod );        
        
            // REVERSE HAIRDO
            if ( i_doReverseHairdo )
            {
                m_rod->doReverseHairdo(stepper);
            }
        } 
        else
	    {
		    RodTimeStepper* stepper = new RodTimeStepper( *m_rod );
		    stepper->setDiffEqSolver( i_solverType );
	    
		    // These get deleted by RodTimeStepper
		    stepper->addExternalForce( new RodMassDamping( m_massDamping ) );
		    
		    if ( i_gravity.norm() > 0)
		    {
		      stepper->addExternalForce( new RodGravity( i_gravity ) );
		      m_gravity = i_gravity;
		    }
	    
		    m_stepper = new RodCollisionTimeStepper( stepper, m_rod );        
	    }
    }
    
    m_isPlaceHolderRod = false;

    // Create space to store the data for each input cv for the previous, current and next positions
    // Used in substepping and for collisions.
    allocateStorage( (int)i_rodVertexPositions.size() );
    resetVertexPositions( i_rodVertexPositions );

}

RodData::~RodData()
{
    if ( m_rod != NULL )
        delete m_rod;
    if ( m_stepper != NULL )
        delete m_stepper;
    if ( m_rodRenderer != NULL )
        delete m_rodRenderer;
    
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

void RodData::addKinematicEdge( unsigned int i_edgeNumber, MaterialFrame* i_materialframe )
{
    removeKinematicEdge( i_edgeNumber );
    KinematicEdgeData* kinematicEdgeData = new KinematicEdgeData( i_edgeNumber, m_rod, i_materialframe );
    kinematicEdgeDataMap[ i_edgeNumber ] = kinematicEdgeData;
}

void RodData::resetKinematicEdge( unsigned int i_edgeNumber, MaterialFrame& i_materialframe )
{
    if ( kinematicEdgeDataMap.find( i_edgeNumber ) != kinematicEdgeDataMap.end() )
    {
        kinematicEdgeDataMap[ i_edgeNumber ]->resetMaterialFrame( m_rod, i_materialframe );
    }
}

void RodData::updateKinematicEdge( unsigned int i_edgeNumber, MaterialFrame& i_materialframe )
{
    if ( kinematicEdgeDataMap.find( i_edgeNumber ) != kinematicEdgeDataMap.end() )
    {
        kinematicEdgeDataMap[ i_edgeNumber ]->updateMaterialFrame( i_materialframe );
    }
}

void RodData::allocateStorage( int i_numCVs )
{
    // Make sure we have enough space to store the date for each CV.
    undeformedVertexPositions.resize( i_numCVs );
    initialVertexPositions.resize( i_numCVs );
    prevVertexPositions.resize( i_numCVs );
    currVertexPositions.resize( i_numCVs );
    nextVertexPositions.resize( i_numCVs );

    undeformedMaterialFrame.resize( i_numCVs - 1 );
}

void RodData::resetVertexPositions( vector< BASim::Vec3d >& i_vertexPositions )
{
    undeformedVertexPositions = i_vertexPositions;
    initialVertexPositions = i_vertexPositions;
    prevVertexPositions = i_vertexPositions;
    currVertexPositions = i_vertexPositions;
    nextVertexPositions = i_vertexPositions;
}

void RodData::updateNextRodVertexPositions( vector< BASim::Vec3d >& i_vertexPositions )
{
    prevVertexPositions = currVertexPositions;

    //FIXME: Doesn't currVertexPositions already equal nextVertexPositions at the end of a time
    // step?
    currVertexPositions = nextVertexPositions;

    nextVertexPositions = i_vertexPositions;
}

void RodData::updateBoundaryConditions()
{
    // Look and see if we have any user controlled vertices or edges, if so we
    // need to apply those as boundary conditions
    
    RodBoundaryCondition* boundary = m_stepper->getBoundaryCondition();
                
    // It looks like we need to erase boundary conditions that were previously set. 
    // When did this happen, it didn't used to be like that?!

    for ( int v=0; v<m_rod->nv(); ++v )
    {
        // Probably most of these were not fixed, this may be slow to do.
        // ??? Benchmark this and check what is going on.
        boundary->releaseVertex( v );
    }
    

    int b = 0;

    for ( KinematicEdgeDataMap::iterator it = kinematicEdgeDataMap.begin();
            it != kinematicEdgeDataMap.end();
            it++ )
    {
        // First make sure these vertices are marked as fixed on the rod
        // or they'll get taken into account on collision calculations.
        unsigned int edgeNum = it->first;
        KinematicEdgeData* kinematicEdgeData = it->second;
        m_rod->fixEdge( edgeNum );
        m_rod->fixVert( edgeNum );
        m_rod->fixVert( edgeNum +1 );

        boundary->setDesiredVertexPosition( edgeNum, currVertexPositions[ edgeNum ] );
        boundary->setDesiredVertexPosition( edgeNum + 1, currVertexPositions[ edgeNum + 1 ]);
        
        // js - TO BE DELETED
        // i don't want to simulate when exporting positions 
        if (0) {
          BASim::Vec3d x = currVertexPositions[ edgeNum + 1 ] - m_rod->getVertex(1);
          for(int i=2; i<m_rod->nv(); i++) {
            boundary->setDesiredVertexPosition( i, m_rod->getVertex(i) + x );
          }
        }

/*
        if ( kinematicEdgeData->rootFrameDefined )
        {
            MaterialFrame m;

            m.m1 = m_rod->getMaterial1( edgeNum );
            m.m2 = m_rod->getMaterial2( edgeNum );
            m.m3 = m_rod->getEdge( edgeNum );
            m.m3.normalize();
            
            MaterialFrame vm = kinematicEdgeData->materialFrame;
            
            vm.m1 = m_rod->getReferenceDirector1( edgeNum );
            vm.m2 = m_rod->getReferenceDirector2( edgeNum );
            vm.m3 = m.m3;
            
            // work out the angle that the material frame is off the reference frame by.
            //cerr << "kinematicEdgeData->getAngleFromRodmaterialFrame() = " << kinematicEdgeData->getAngleFromRodmaterialFrame() << endl;

            boundary->setDesiredEdgeAngle( edgeNum, kinematicEdgeData->getAngleFromRodmaterialFrame() );
        }
        else
        {
            boundary->setDesiredEdgeAngle( edgeNum, m_rod->getTheta( edgeNum ) );
        }
*/
    }
}















    double KinematicEdgeData::getAngleBetweenReferenceAndStrand()
    { 
          double pi = 3.14159265358979323846264338327950288;
          double theta = 0.0;
          
          if ( rod != NULL )
          {   
                BASim::Vec3d m = rod->getMaterial2( edgeNumber );
                m.normalize();
                BASim::Vec3d b = materialFrame.m3;
                b.normalize();
    
                BASim::Vec3d x = rod->getEdge( edgeNumber );
                x.normalize();
                            
                BASim::Vec3d c = m.cross( b );
                double dot = x.dot(c);
                double sign;
                if ( dot == 0 )
                    sign = 0;
                else if ( dot > 0 )
                    sign = 1;
                else {
                    // ? Why was this an assert when it's dealt with below?
                   // assert( dot < 0 );
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
            BASim::Vec3d rodRefMaterial = rod->getMaterial2( edgeNumber );
            
            double pi = 3.14159265358979323846264338327950288;

            BASim::Vec3d x = rod->getEdge( edgeNumber );
            x.normalize();
            BASim::Vec3d r = rodRefMaterial;
            r.normalize();
            //BASim::Vec3d m = rod->getMaterial2( edgeNumber );
            BASim::Vec3d m = materialFrame.m3;
            m.normalize();
            
            BASim::Vec3d c = r.cross( m );
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
