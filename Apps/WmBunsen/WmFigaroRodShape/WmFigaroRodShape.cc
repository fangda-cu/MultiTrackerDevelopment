#include <math.h>           
#include <maya/MIOStream.h>
#include <maya/MFnSingleIndexedComponent.h>
#include <maya/MFnCompoundAttribute.h>
#include <maya/MFnPointArrayData.h>
#include <maya/MDagPath.h>


#include <weta/Wvec3.hh>
#include <weta/Wvec4.hh>
#include <weta/Wnunbsc.hh>

#include "WmFigaroRodShape.hh"         
#include "WmFigaroRodShapeUI.hh"     
#include "WmFigaroRodShapeIterator.hh"


////////////////////////////////////////////////////////////////////////////////
//
// Shape implementation
//
////////////////////////////////////////////////////////////////////////////////

MTypeId WmFigaroRodShape::id( ( 0x00101A, 0x07 ) );
MString WmFigaroRodShape::typeName( "wmFigaroRodShape" );
MObject WmFigaroRodShape::oa_cv;
MObject WmFigaroRodShape::oa_cvX;
MObject WmFigaroRodShape::oa_cvY;
MObject WmFigaroRodShape::oa_cvZ;

MObject WmFigaroRodShape::i_bboxCorner1;
MObject WmFigaroRodShape::i_bboxCorner2;

MObject WmFigaroRodShape::ia_triggerSolve;
MObject WmFigaroRodShape::ia_doPlasticDeformations;
MObject WmFigaroRodShape::ia_constraintStrength;
MObject WmFigaroRodShape::ia_percentLocked;
MObject WmFigaroRodShape::ia_cylinderDrawScale;
MObject WmFigaroRodShape::ia_youngsModulus;
MObject WmFigaroRodShape::ia_shearModulus;
MObject WmFigaroRodShape::ia_viscosity;
MObject WmFigaroRodShape::ia_density;
MObject WmFigaroRodShape::ia_radiusA;
MObject WmFigaroRodShape::ia_radiusB;
MObject WmFigaroRodShape::ia_vertexSpacing;

MObject WmFigaroRodShape::ia_numberOfIterations;

MObject WmFigaroRodShape::ca_sync;

WmFigaroRodShape::WmFigaroRodShape() : m_initialisedRod( false ), m_constraintStrength( 20.0 ), 
    m_drawRodTube( false ), m_numberOfLockedVertices( 0 ), m_cylinderDrawScale( 10.0 ),
    m_youngsModulus( 100.0 * 1e7 /* pascal */ ), m_shearModulus( 34.0 * 1e7  /* pascal */ ),
    m_viscosity( 10.0 /* poise */ ), m_density( 1.3 /* grams per cubic centimeter */ ),
    m_radiusA( 0.05 * 1e-1 /* centimetre */ ), m_radiusB( 0.05 * 1e-1 /* centimetre */ ),
    m_numberOfIterations( 6 ), m_vertexSpacing(1.0)
{        
}

WmFigaroRodShape::~WmFigaroRodShape() {}

void WmFigaroRodShape::resetSimulation( MVectorArray* i_controlPoints )
{
    m_beaker.resetEverything();

    MVectorArray* controlPoints = getControlPoints();

    if ( i_controlPoints == NULL )
    {
        initialiseRod( controlPoints );   
    }
    else
    {
        initialiseRod( i_controlPoints );
    }

    m_rodGroup.setSimulationNeedsReset( false );
    
    m_beaker.addRodsToWorld( 0, &m_rodGroup );
}

void WmFigaroRodShape::initialiseRod( MVectorArray* i_controlPoints )
{
    MStatus stat;

    MVectorArray vertices;

    m_rodGroup.removeAllRods();

    vector<Vec3d> rodVertices;
    rodVertices.resize( i_controlPoints->length() );
    //cout<<" add new rod : "<<i_controlPoints->length()<<" number of CVs "<<endl;
    
    for ( int v=0; v<i_controlPoints->length(); v++ )
    {
        rodVertices[ v ] = Vec3d( (*i_controlPoints)[ v ].x, (*i_controlPoints)[ v ].y, (*i_controlPoints)[ v ].z );
        //cout<<" -- > cv : "<<v<<" : "<<(*i_controlPoints)[ v ]<<endl;
    }

    double massDamping = 10.0;
    //Vec3d gravity( 0.0, -981.0, 0.0 );
    Vec3d gravity( 0.0, 0.0, 0.0 );
    RodTimeStepper::Method solverType = RodTimeStepper::SYM_IMPL_EULER;
    //RodTimeStepper::Method solverType = RodTimeStepper::STATICS;
    RodOptions rodOptions;
    rodOptions.YoungsModulus = m_youngsModulus;
    rodOptions.ShearModulus = m_shearModulus;
    rodOptions.viscosity = m_viscosity;
    rodOptions.density = m_density;
    rodOptions.radiusA = m_radiusA;
    rodOptions.radiusB = m_radiusB;
    rodOptions.refFrame = BASim::ElasticRod::TimeParallel;
    rodOptions.numVertices = (int)rodVertices.size();

    size_t rodIndex = m_rodGroup.addRod( rodVertices, rodOptions, massDamping, gravity, solverType );

//    m_rodGroup.addKinematicEdge( rodIndex, 0 );

    // This will get passed to Beaker later, we are saying that the 0th verex must remain locked
    // where it started.
    m_lockedVertexMap[ 0 ] = rodVertices[ 0 ];

    m_rodGroup.setDrawMode( RodRenderer::SMOOTH );
    m_rodGroup.setDrawScale( 50.0 );

    // Store all the rod edge lengths so we can stop the user moving a vertex and stretching an edge
    m_edgeLengths.resize( rodVertices.size() -1 );
    for ( size_t e=0; e<rodVertices.size() - 1; ++e )
    {
        m_edgeLengths[ e ] = ( rodVertices[ e + 1 ] - rodVertices[ e ] ).norm();        
//        cerr << "m_edgeLength[ e ] = " << m_edgeLengths[ e ] << endl;        
    }

    m_initialisedRod = true;
}

void WmFigaroRodShape::solve( FixedRodVertexMap& i_fixedRodVertexMap )
{
    LockedRodVertexMap lockedRodVertexMap;
    lockedRodVertexMap[ 0 ] = m_lockedVertexMap;
    m_beaker.takeTimeStep( 8, 1.0/24.0, 10, false, false, false, 10, 1.0, &i_fixedRodVertexMap, false, m_constraintStrength, &lockedRodVertexMap );
}

void WmFigaroRodShape::postConstructor()
{
    // The control points are now set by the command when we are generated from a nurbs curve,
    

    // Set the control points of the shape to match the vertices of the rod
    // no need to set them here.
    //MVectorArray controlPoints;
    //getRodVertices( controlPoints );
    //buildControlPoints( controlPoints );
}

void WmFigaroRodShape::getRodVertices( MVectorArray& o_controlPoints )
{
    if ( m_rodGroup.numberOfRods() == 0 )
    {
        return;
    }

    ElasticRod* rod = m_rodGroup.elasticRod( 0 );

    o_controlPoints.setLength( rod->nv() );
    for ( int p=0; p<rod->nv(); ++p )
    {
        Vec3d vertex = rod->getVertex( p );
        o_controlPoints[ p ] = MVector( vertex[ 0 ], vertex[ 1 ], vertex[ 2 ] );
    }
}

MStatus WmFigaroRodShape::buildControlPoints( MVectorArray& i_controlPoints )
{
    MStatus stat;

    MDataBlock datablock = forceCache();

    MArrayDataHandle cpH = datablock.outputArrayValue( mControlPoints, &stat );
    MCHECKERROR( stat, "compute get cpH" )
        
    // Create a builder to aid in the array construction efficiently.
    //
    MArrayDataBuilder bOutArray = cpH.builder( &stat );
    CHECK_MSTATUS( stat );

    for ( int vtx=0; vtx<i_controlPoints.length(); vtx++ )
    {
        double3 &pt = bOutArray.addElement( vtx ).asDouble3();
        pt[ 0 ] = i_controlPoints[ vtx ][ 0 ];
        pt[ 1 ] = i_controlPoints[ vtx ][ 1 ];
        pt[ 2 ] = i_controlPoints[ vtx ][ 2 ];
    }

    cpH.set( bOutArray );
        
    cpH.setAllClean();

    return stat;
}

void WmFigaroRodShape::updateControlPointsFromRod()
{
    MStatus status;

    MVectorArray controlPoints;
    getRodVertices( controlPoints );

    status = setControlPoints( &controlPoints );
    CHECK_MSTATUS( status );

    // Now get the input attribute of the rod and set
    MPlug controlPointsPlugArr( thisMObject(), oa_cv );
    CHECK_MSTATUS( status );

    for ( int p=0; p<controlPoints.length(); ++p )
    {
        MPlug cvPlug = controlPointsPlugArr.elementByLogicalIndex( p, &status );
        CHECK_MSTATUS( status );
        MPlug xPlug = cvPlug.child( 0, &status );
        CHECK_MSTATUS( status );
        MPlug yPlug = cvPlug.child( 1, &status );
        CHECK_MSTATUS( status );
        MPlug zPlug = cvPlug.child( 2, &status );
        CHECK_MSTATUS( status );
        xPlug.setValue( controlPoints[ p ].x );
        yPlug.setValue( controlPoints[ p ].y );
        zPlug.setValue( controlPoints[ p ].z );        
    }
}

void WmFigaroRodShape::drawRod()
{
    if ( m_rodGroup.numberOfRods() == 0 )
    {
        return;
    }

    m_rodGroup.render();

    RodRenderer* rodRenderer = m_rodGroup.rodRenderer( 0 );
    rodRenderer->drawSmoothPartialRod( 0, m_numberOfLockedVertices, Vec3d( 1.0, 0.0, 0.0 ) );
}

void* WmFigaroRodShape::creator()
//
// Description
//
//    Called internally to create a new instance of the users MPx node.
//
{
	return new WmFigaroRodShape();
}

MStatus WmFigaroRodShape::initialize()
//
// Description
//
//    Attribute (static) initialization.
//    See api_macros.h.
//
{ 
    MStatus status;

    // bbox attributes
    //
    MAKE_NUMERIC_ATTR(  i_bboxCorner1, "bboxCorner1", "bb1",
                        MFnNumericData::k3Double, 0,
                        false, false, false );
    MAKE_NUMERIC_ATTR(  i_bboxCorner2, "bboxCorner2", "bb2",
                        MFnNumericData::k3Double, 0,
                        false, false, false );

    {
        MFnNumericAttribute nAttr;
        oa_cvX = nAttr.create( "cvX", "cx", MFnNumericData::kDouble, 0.0, & status );
        nAttr.setConnectable( true );
        status = addAttribute( oa_cvX );
        oa_cvY = nAttr.create( "cvY", "cy", MFnNumericData::kDouble, 0.0, & status );
        nAttr.setConnectable( true );
        status = addAttribute( oa_cvY );
        oa_cvZ = nAttr.create( "cvZ", "cz", MFnNumericData::kDouble, 0.0, & status );
        nAttr.setConnectable( true );
        status = addAttribute( oa_cvZ );
    }

    {
        MFnCompoundAttribute cAttr;
        oa_cv = cAttr.create( "controlVertex", "cv" );
        cAttr.addChild( oa_cvX );
        cAttr.addChild( oa_cvY );
        cAttr.addChild( oa_cvZ );
        cAttr.setArray( true );
        status = addAttribute( oa_cv );
        CHECK_MSTATUS( status );
    }

    attributeAffects( mControlPoints, oa_cv ); 
    attributeAffects( mControlPoints, oa_cvX );
    attributeAffects( mControlPoints, oa_cvY );
    attributeAffects( mControlPoints, oa_cvZ );

    {
        MFnNumericAttribute nAttr;
        ca_sync = nAttr.create( "sync", "syn", MFnNumericData::kDouble, 0.0, & status );
        nAttr.setConnectable( true );
        status = addAttribute( ca_sync );
    }

    /////////////////////////////////////////////////////////////////////
    //
    // Inputs
    //
    ////////////////////////////////////////////////////////////////////

    {
        MFnNumericAttribute nAttr;
        ia_triggerSolve = nAttr.create( "triggerSolve", "ts", MFnNumericData::kDouble, 0.0, & status );
        nAttr.setConnectable( true );
        nAttr.setReadable( false );
        nAttr.setWritable(true );
        status = addAttribute( ia_triggerSolve );        
    }

    {
        MFnNumericAttribute nAttr;
        ia_doPlasticDeformations = nAttr.create( "doPlasticDeformations", "dpd", MFnNumericData::kBoolean, false, & status );
        nAttr.setConnectable( true );
        nAttr.setReadable( false );
        nAttr.setWritable(true );
        status = addAttribute( ia_doPlasticDeformations );        
    }

    {
        MFnNumericAttribute nAttr;
        ia_constraintStrength = nAttr.create( "constraintStrength", "cos", MFnNumericData::kDouble, 10.0, & status );
        nAttr.setConnectable( true );
        nAttr.setReadable( false );
        nAttr.setWritable(true );
        status = addAttribute( ia_constraintStrength );        
    }

    {
        MFnNumericAttribute nAttr;
        ia_percentLocked = nAttr.create( "percentLocked", "plo", MFnNumericData::kInt, 0, & status );
        nAttr.setConnectable( true );
        nAttr.setReadable( false );
        nAttr.setWritable(true );
        nAttr.setMin( 0 );
        nAttr.setMax( 100 );
        status = addAttribute( ia_percentLocked );        
    }

    {
        MFnNumericAttribute nAttr;
        ia_cylinderDrawScale = nAttr.create( "cylinderDrawScale", "cds", MFnNumericData::kDouble, 10.0, & status );
        nAttr.setConnectable( true );
        nAttr.setReadable( false );
        nAttr.setWritable(true );
        nAttr.setMin( 1 );
        nAttr.setMax( 50 );
        status = addAttribute( ia_cylinderDrawScale );        
    }

    {
        MFnNumericAttribute nAttr;
        ia_youngsModulus = nAttr.create( "youngsModulus", "ymo", MFnNumericData::kDouble, 100.0 , & status );
        nAttr.setConnectable( true );
        nAttr.setReadable( false );
        nAttr.setWritable(true );
        status = addAttribute( ia_youngsModulus );        
    }

    {
        MFnNumericAttribute nAttr;
        ia_shearModulus = nAttr.create( "shearModulus", "sho", MFnNumericData::kDouble, 34.0, & status );
        CHECK_MSTATUS( status );
        nAttr.setConnectable( true );
        nAttr.setReadable( false );
        nAttr.setWritable(true );
        status = addAttribute( ia_shearModulus );        
        CHECK_MSTATUS( status );
    }

    {
        MFnNumericAttribute nAttr;
        ia_viscosity = nAttr.create( "viscosity", "vis", MFnNumericData::kDouble, 10.0, & status );
        nAttr.setConnectable( true );
        nAttr.setReadable( false );
        nAttr.setWritable(true );
        status = addAttribute( ia_viscosity );        
    }

    {
        MFnNumericAttribute nAttr;
        ia_density = nAttr.create( "density", "den", MFnNumericData::kDouble, 1.3, & status );
        nAttr.setConnectable( true );
        nAttr.setReadable( false );
        nAttr.setWritable(true );
        status = addAttribute( ia_density );        
    }

    {
        MFnNumericAttribute nAttr;
        ia_radiusA = nAttr.create( "minorRadius", "mir", MFnNumericData::kDouble, 0.05, & status );
        nAttr.setConnectable( true );
        nAttr.setReadable( false );
        nAttr.setWritable(true );
        status = addAttribute( ia_radiusA );        
    }

    {
        MFnNumericAttribute nAttr;
        ia_radiusB = nAttr.create( "majorRadius", "mar", MFnNumericData::kDouble, 0.05, & status );
        nAttr.setConnectable( true );
        nAttr.setReadable( false );
        nAttr.setWritable(true );
        status = addAttribute( ia_radiusB );        
    }

    {
        MFnNumericAttribute nAttr;
        ia_numberOfIterations = nAttr.create( "numberOfIterations", "noi", MFnNumericData::kInt, 6, & status );
        nAttr.setConnectable( true );
        nAttr.setReadable( false );
        nAttr.setWritable(true );
        status = addAttribute( ia_numberOfIterations );        
    }

    {
         MFnNumericAttribute nAttr;
         ia_vertexSpacing = nAttr.create( "vertexSpacing", "vsp", MFnNumericData::kDouble, 1.0, &status );
         nAttr.setConnectable( true );
         nAttr.setReadable( false );
         nAttr.setWritable(true );
         status = addAttribute( ia_vertexSpacing );
    }

    attributeAffects( ia_numberOfIterations, ca_sync );
    attributeAffects( ia_numberOfIterations, oa_cv ); 
    attributeAffects( ia_numberOfIterations, oa_cvX );
    attributeAffects( ia_numberOfIterations, oa_cvY );
    attributeAffects( ia_numberOfIterations, oa_cvZ );
    

    attributeAffects( ia_youngsModulus, ca_sync );
    attributeAffects( ia_youngsModulus, oa_cv ); 
    attributeAffects( ia_youngsModulus, oa_cvX );
    attributeAffects( ia_youngsModulus, oa_cvY );
    attributeAffects( ia_youngsModulus, oa_cvZ );
    attributeAffects( ia_shearModulus, ca_sync );
    attributeAffects( ia_shearModulus, oa_cv ); 
    attributeAffects( ia_shearModulus, oa_cvX );
    attributeAffects( ia_shearModulus, oa_cvY );
    attributeAffects( ia_shearModulus, oa_cvZ );
    attributeAffects( ia_viscosity, ca_sync );
    attributeAffects( ia_viscosity, oa_cv ); 
    attributeAffects( ia_viscosity, oa_cvX );
    attributeAffects( ia_viscosity, oa_cvY );
    attributeAffects( ia_viscosity, oa_cvZ );
    attributeAffects( ia_density, ca_sync );
    attributeAffects( ia_density, oa_cv ); 
    attributeAffects( ia_density, oa_cvX );
    attributeAffects( ia_density, oa_cvY );
    attributeAffects( ia_density, oa_cvZ );
    attributeAffects( ia_radiusA, ca_sync );
    attributeAffects( ia_radiusA, oa_cv ); 
    attributeAffects( ia_radiusA, oa_cvX );
    attributeAffects( ia_radiusA, oa_cvY );
    attributeAffects( ia_radiusA, oa_cvZ );
    attributeAffects( ia_radiusB, ca_sync );
    attributeAffects( ia_radiusB, oa_cv ); 
    attributeAffects( ia_radiusB, oa_cvX );
    attributeAffects( ia_radiusB, oa_cvY );
    attributeAffects( ia_radiusB, oa_cvZ );

    attributeAffects( ia_triggerSolve, ca_sync );
    attributeAffects( ia_triggerSolve, oa_cv ); 
    attributeAffects( ia_triggerSolve, oa_cvX );
    attributeAffects( ia_triggerSolve, oa_cvY );
    attributeAffects( ia_triggerSolve, oa_cvZ );

    attributeAffects( ia_doPlasticDeformations, ca_sync ); 
    attributeAffects( ia_doPlasticDeformations, oa_cv ); 
    attributeAffects( ia_doPlasticDeformations, oa_cvX );
    attributeAffects( ia_doPlasticDeformations, oa_cvY );
    attributeAffects( ia_doPlasticDeformations, oa_cvZ );

    attributeAffects( ia_constraintStrength, ca_sync ); 
    attributeAffects( ia_constraintStrength, oa_cv ); 
    attributeAffects( ia_constraintStrength, oa_cvX );
    attributeAffects( ia_constraintStrength, oa_cvY );
    attributeAffects( ia_constraintStrength, oa_cvZ );
    
    attributeAffects( ia_percentLocked, ca_sync ); 
    attributeAffects( ia_percentLocked, oa_cv ); 
    attributeAffects( ia_percentLocked, oa_cvX );
    attributeAffects( ia_percentLocked, oa_cvY );
    attributeAffects( ia_percentLocked, oa_cvZ );

    attributeAffects( ia_cylinderDrawScale, ca_sync ); 
    attributeAffects( ia_cylinderDrawScale, oa_cv ); 
    attributeAffects( ia_cylinderDrawScale, oa_cvX );
    attributeAffects( ia_cylinderDrawScale, oa_cvY );
    attributeAffects( ia_cylinderDrawScale, oa_cvZ );

    attributeAffects( ia_vertexSpacing, ca_sync );
    attributeAffects( ia_vertexSpacing, oa_cv );
    attributeAffects( ia_vertexSpacing, oa_cvX );
    attributeAffects( ia_vertexSpacing, oa_cvY );
    attributeAffects( ia_vertexSpacing, oa_cvZ );
    

    attributeAffects( ca_sync, oa_cv ); 
    attributeAffects( ca_sync, oa_cvX );
    attributeAffects( ca_sync, oa_cvY );
    attributeAffects( ca_sync, oa_cvZ );

/*    {
        MFnNumericAttribute nAttr;
        ia_inCVX = nAttr.create( "inCVX", "icx", MFnNumericData::kDouble, 0.0, & status );
        nAttr.setConnectable( true );
        status = addAttribute( oa_cvX );
        ia_inCVY = nAttr.create( "inCVY", "icy", MFnNumericData::kDouble, 0.0, & status );
        nAttr.setConnectable( true );
        status = addAttribute( oa_cvY );
        ia_inCVZ = nAttr.create( "inCVZ", "icz", MFnNumericData::kDouble, 0.0, & status );
        nAttr.setConnectable( true );
        status = addAttribute( oa_cvZ );
    }

    {
        MFnCompoundAttribute cAttr;
        ia_inCV = cAttr.create( "inControlVertex", "icv" );
        cAttr.addChild( ia_inCVX );
        cAttr.addChild( ia_inCVY );
        cAttr.addChild( ia_inCVZ );
        cAttr.setArray( true );
        status = addAttribute( ia_inCV );
        CHECK_MSTATUS( status );
    }*/

    return MS::kSuccess;
}

void WmFigaroRodShape::getResampledRodCVs(bool i_needResample, MVectorArray &o_controlCVs )
{

    ElasticRod* elasticRod = m_rodGroup.elasticRod( 0 );
    int numPts  =  elasticRod->nv();
  
    o_controlCVs.setLength( numPts );

    MVectorArray rodCVs;
    rodCVs.setLength( numPts );
    Wvec3d currentCVs[ numPts ];

    for ( int v=0; v<numPts; ++v )
    {
        Vec3d vertex = elasticRod->getVertex( v );
        rodCVs[ v ] = MVector( vertex[ 0 ], vertex[ 1 ], vertex[ 2 ] );
        currentCVs[ v  ] = Wvec3d( vertex[ 0 ], vertex[ 1 ], vertex[ 2 ]  );
        o_controlCVs[ v ] = MVector( vertex[ 0 ], vertex[ 1 ], vertex[ 2 ] );
    }

    if( !i_needResample ||  m_vertexSpacing <= 0 )
    {
        // if do not need resample or the vertex spacing is a invalid value
        // simple copy the old rod CVs and return. 
        return;
    }


    // we need to resample so, clear it.
    o_controlCVs.clear();

    //build a curve for the old rod, cause we need to figure out
    // how many number of cvs for the new rod given different
    // value of vertex spacing. 
    Wnunbsc cubicCurve, sampleCurve;
    cubicCurve.centripetalbesselcubicinterp( currentCVs, numPts );
    cubicCurve.calcpointsperspan( sampleCurve, 3 );

    //new number of cv!
    numPts = int( sampleCurve.chordlength() / m_vertexSpacing );

    // resample!
    resampleRodCVsForControlVertices( rodCVs, numPts, o_controlCVs );
   
}

//
// given a bunch of CVs "i_rodCVs"  either from rods or from curves. and a different number
// "i_numCVs", output the resample CVs "o_outputCVs" based on the input.
//
void WmFigaroRodShape::resampleRodCVsForControlVertices( MVectorArray &i_rodCVs,
          int i_numCVs,
        MVectorArray &o_outputCVs )
{

    MStatus status;
    CHECK_MSTATUS( status );

    int numRodCVs = i_rodCVs.length();
    Wvec3d currentCVs[ numRodCVs ];
    for( int i = 0; i <  numRodCVs; i++ )
    {
        currentCVs[ i ] = Wvec3d( i_rodCVs[i].x,
                i_rodCVs[i].y, i_rodCVs[i].z );

    }

    Wnunbsc cubicCurve, sampleCurve;
    cubicCurve.centripetalbesselcubicinterp( currentCVs, numRodCVs );
    cubicCurve.calcpointsperspan( sampleCurve, 3 );


    double start = sampleCurve.Knot.K[ sampleCurve.Knot.Deg   - 1 ];
    double end   = sampleCurve.Knot.K[ sampleCurve.Knot.NumCV - 1 ];

    double step = ( end - start ) / double( i_numCVs - 1 );
    double t = 0;
    Wvec3d pt;
    for( int i = 0; i < i_numCVs; i++ )
    {
        if ( i == 0 )
        {
            t = start;
        }
        else if( i == i_numCVs - 1  )
        {
             t = end;
        }
        else
        {
            t = start + double(i)*step;
        }

        if ( t > end )
        {
            t = end;
        }

        sampleCurve.calcpoint( pt, t );
        o_outputCVs.append( MVector( pt.X, pt.Y, pt.Z));

    }

}


MStatus WmFigaroRodShape::compute( const MPlug& i_plug, MDataBlock& i_dataBlock )
{
    cerr << "shape called with plug " << i_plug.name() << endl;

    MStatus status;

    if( i_plug == oa_cv || i_plug == oa_cvX ||
       i_plug == oa_cvY || i_plug == oa_cvZ )
    {
        i_dataBlock.inputValue( ca_sync, &status );
        CHECK_MSTATUS( status );

        if ( m_initialisedRod )
        {
            i_dataBlock.inputValue( ia_triggerSolve, &status ).asDouble();
            CHECK_MSTATUS( status );            

            // FIXME: Why not use the datablock....
            MDataHandle pointsH = i_dataBlock.inputValue( mControlPoints, &status );
            CHECK_MSTATUS( status );
            /*MObject pointsDataObj = pointsH.data();
            MFnPointArrayData pointArrayData( pointsDataObj );
            MPointArray pointArray = pointArrayData.array( &status );
            CHECK_MSTATUS( status );*/

            MVectorArray& pointArray = *( getControlPoints() );

            MPlug cvArrPlug( thisMObject(), oa_cv );
            int numOutputCvs = cvArrPlug.numElements();
            MVectorArray outputCVs;
            resampleRodCVsForControlVertices( pointArray, numOutputCvs,  outputCVs );

	    //need to apply the inverse transformation of the curve for the outputCVS
	    MPlugArray plugs;
            i_plug.connectedTo( plugs, false, true, &status);
	    CHECK_MSTATUS( status );
	    if( plugs.length() > 0 )
	    {
		MFnDagNode dagFn( plugs[0].node(), &status);
		CHECK_MSTATUS( status );
	        MDagPath path;
		dagFn.getPath( path );
		CHECK_MSTATUS( status );
		MMatrix xform = path.inclusiveMatrix();
                MMatrix invXform = xform.inverse();
                
		for( unsigned int i = 0; i < outputCVs.length(); i++ )
                {
                        MPoint pt = outputCVs[ i ];
                        pt *= invXform;
                        outputCVs[ i ] = MVector( pt.x, pt.y, pt.z );
		}
		
	    }		

            for( unsigned int cv = 0; cv < outputCVs.length(); cv++ )
            {
                MPlug cvPlug = cvArrPlug.elementByLogicalIndex( cv );
                cvPlug.child( 0, & status ).setDouble( outputCVs[ cv ].x );
                cvPlug.child( 1, & status ).setDouble( outputCVs[ cv ].y );
                cvPlug.child( 2, & status ).setDouble( outputCVs[ cv ].z );

                //cerr << "cv " << cv << " == " << pointArray[ cv ] << endl;
            }
        }

        i_dataBlock.setClean( i_plug );

        return MS::kSuccess;
    }
    else if ( i_plug == ca_sync )
    {
        m_youngsModulus = i_dataBlock.inputValue( ia_youngsModulus ).asDouble() * 1e7; /* pascal */
        m_shearModulus = i_dataBlock.inputValue( ia_shearModulus ).asDouble() * 1e7;  /* pascal */
        m_viscosity = i_dataBlock.inputValue( ia_viscosity ).asDouble(); /* poise */ 
        m_density = i_dataBlock.inputValue( ia_density ).asDouble(); /* grams per cubic centimeter */
        m_radiusA = i_dataBlock.inputValue( ia_radiusA ).asDouble() * 1e-1; /* centimetre */
        m_radiusB = i_dataBlock.inputValue( ia_radiusB ).asDouble() * 1e-1; /* centimetre */

        bool needResample = ( m_vertexSpacing != i_dataBlock.inputValue( ia_vertexSpacing ).asDouble()) ;

        if( needResample )
        {
            m_vertexSpacing = i_dataBlock.inputValue( ia_vertexSpacing ).asDouble();
        }
        //m_numberOfIterations = i_dataBlock.inputValue( ia_numberOfIterations ).asInt();
        
        // Rebuild the rod in this new configuration
        if ( !m_initialisedRod )
        {
            resetSimulation();
        }
        else
        {
            //resample the rod CVs if we need resample ( when vertex spacing is changed)
            MVectorArray rodVertices;
            getResampledRodCVs( needResample, rodVertices );
            resetSimulation( &rodVertices );
        }

        i_dataBlock.inputValue( ia_triggerSolve );

        m_cylinderDrawScale = i_dataBlock.inputValue( ia_cylinderDrawScale ).asDouble();        

        // First we need to flag the vertex that was passed to compute as a fixed point
        // FIXME: These fixed vertices should be in the rodData class so they are attached
        // to a specific rod. Have a set of fixed vertices and a set of constrained Vertices

        m_constraintStrength = i_dataBlock.inputValue( ia_constraintStrength ).asDouble();        

        int percentLocked = i_dataBlock.inputValue( ia_percentLocked ).asInt();
        
        // We fix the 2nd vertex with one of these constraints so that it tries to maintain
        // the flow direction
        ElasticRod* elasticRod = m_rodGroup.elasticRod( 0 );
        m_fixedVertexMap[ 1 ] = elasticRod->getVertex( 1 );

        // Set the display scale to match what the user specified
        elasticRod->setRadiusScale( m_cylinderDrawScale );

        // Fix all vertices in the range 0 -> m_percentLocked 
        double length = 0.0;
        for ( int v=2; v<elasticRod->nv(); ++v )
        {
            length += ( elasticRod->getVertex( v ) - elasticRod->getVertex( v - 1 ) ).norm();
        }

        double lockedLength = ( percentLocked / 100.0 ) * length;

        m_numberOfLockedVertices = 1;

        length = 0.0;
        for ( int v=2; v<elasticRod->nv(); ++v )
        {
            length += ( elasticRod->getVertex( v ) - elasticRod->getVertex( v - 1 ) ).norm();
            if ( length < lockedLength )
            {
                cerr << "locking vertex " << v << endl;
                m_fixedVertexMap[ v ] = elasticRod->getVertex( v );

                m_numberOfLockedVertices = v;
            }
        }

        FixedRodVertexMap fixedRodVertexMap;
        fixedRodVertexMap[ 0 ] = m_fixedVertexMap;

        cerr << "there are " << m_fixedVertexMap.size() << " vertices fixed\n";
        
        cerr << "DOING SOLVE\n";
        doSolverIterations( fixedRodVertexMap );
        //solve( fixedRodVertexMap );
        cerr << "END SOLVE\n";

        // Only plastic deformations make sense and as we want the user to be able to control
        // all the rod attributes we shall recreate the rod every frame.
        /*if ( i_dataBlock.inputValue( ia_doPlasticDeformations, &status ).asBool() )
        {
            cerr << "doing plastic Deformation\n";
            // Make this a plastic deformation
            m_rodGroup.elasticRod( 0 )->updateReferenceProperties();
        }*/        

        // Now copy the simulated positions of all vertices onto the controlPoints attribute
            
        MVectorArray rodVertices;
        getRodVertices( rodVertices );

        MArrayDataHandle pointsH = i_dataBlock.outputValue( mControlPoints, &status );
        CHECK_MSTATUS( status );

        MArrayDataBuilder newBuilder( &i_dataBlock, mControlPoints, 0, &status );
        CHECK_MSTATUS( status );

        for(unsigned int e=0; e<rodVertices.length(); ++e )
        {
            newBuilder.addElement( e, &status).set( rodVertices[ e ] );
            CHECK_MSTATUS( status );

        }

        pointsH.set( newBuilder );
        pointsH.setAllClean();

        //resampleRodCVsForControlVertices( rodVertices, i_dataBlock, pointsH );


//        unsigned int elementCount = pointsH.elementCount( &status );
//        CHECK_MSTATUS( status );
//
//        if ( rodVertices.length() == elementCount )
//        {
//            for ( unsigned int e=0; e<elementCount; ++e )
//            {
//                status = pointsH.jumpToElement( e );
//                CHECK_MSTATUS( status );
//
//                MDataHandle pointHandle = pointsH.outputValue( &status );
//                CHECK_MSTATUS( status );
//
//                pointHandle.set( rodVertices[ e ].x, rodVertices[ e ].y, rodVertices[ e ]. z );
//            }
//        }

        // Only plastic deformations make sense and as we want the user to be able to control
        // all the rod attributes we shall recreate the rod every frame.
        /*if ( i_dataBlock.inputValue( ia_doPlasticDeformations, &status ).asBool() )
        {
            cerr << "doing plastic Deformation\n";
            // Make this a plastic deformation
            m_rodGroup.elasticRod( 0 )->updateReferenceProperties();
        }*/
        
        m_fixedVertexMap.clear();

        i_dataBlock.setClean( i_plug );

        computeBoundingBox( i_dataBlock );

        return MS::kSuccess;
    }
    
    return MS::kUnknownParameter;
}


MPxGeometryIterator* WmFigaroRodShape::geometryIteratorSetup(MObjectArray& componentList,
													MObject& components,
													bool forReadOnly )
//
// Description
//
//    Creates a geometry iterator compatible with his shape.
//
// Arguments
//
//    componentList - list of components to be iterated
//    components    - component to be iterator
//    forReadOnly   -
//
// Returns
//
//    An iterator for the components
//
{
	WmFigaroRodShapeIterator * result = NULL;
	if ( components.isNull() ) 
	{
		result = new WmFigaroRodShapeIterator( getControlPoints(), componentList );
	}
	else 
	{
		result = new WmFigaroRodShapeIterator( getControlPoints(), components );
	}
	return result;
}

bool WmFigaroRodShape::setInternalValueInContext( const MPlug& i_plug, const MDataHandle& handle, MDGContext& i_context )
{
    MStatus status;

    int index = -1;

    // We get called to set one component of the attribute at a time so we use 999999 to
    // signify Beaker to leave that component where it is rather than moving it to some
    // unknown location when it does the solve
    MPoint newPoint( 999999, 999999, 999999 );

    if (  (i_plug == mControlPoints) ||
          (i_plug == mControlValueX) ||
          (i_plug == mControlValueY) ||
          (i_plug == mControlValueZ) )
    {
        if ( i_plug == mControlPoints && !i_plug.isArray() ) 
        {
            index = i_plug.logicalIndex();
            double3& ptData = handle.asDouble3();
            newPoint.x = ptData[0];
            newPoint.y = ptData[1];
            newPoint.z = ptData[2];

            m_fixedVertexMap[ index ] = Vec3d( newPoint.x, newPoint.y, newPoint.z );
        }
        else if ( i_plug == mControlValueX ) 
        {
            MPlug parentPlug = i_plug.parent();
            index = parentPlug.logicalIndex();
            newPoint.x = handle.asDouble();
        
            m_fixedVertexMap[ index ][0] =newPoint.x;
        }
        else if ( i_plug == mControlValueY ) 
        {
            MPlug parentPlug = i_plug.parent();
            index = parentPlug.logicalIndex();
            newPoint.y = handle.asDouble();

            m_fixedVertexMap[ index ][1] =newPoint.y;
        }
        else if ( i_plug == mControlValueZ ) 
        {
            MPlug parentPlug = i_plug.parent();
            index = parentPlug.logicalIndex();
            newPoint.z = handle.asDouble();

            m_fixedVertexMap[ index ][2] =newPoint.z;
        }
    }

    // Signal that Maya still has to store the data
    return false;    
}

bool WmFigaroRodShape::doSolverIterations( FixedRodVertexMap& i_fixedRodVertexMap )
{
    for ( int i=0; i<m_numberOfIterations; ++i )
    {
        solve( i_fixedRodVertexMap );
    }

    return true;
}

bool WmFigaroRodShape::acceptsGeometryIterator( bool writeable )
//
// Description
//
//    Specifies that this shape can provide an iterator for getting/setting
//    control point values.
//
// Arguments
//
//    writable - maya asks for an iterator that can set points if this is true
//
{
	return true;
}

bool WmFigaroRodShape::acceptsGeometryIterator( MObject&, bool writeable, bool forReadOnly )
//
// Description
//
//    Specifies that this shape can provide an iterator for getting/setting
//    control point values.
//
// Arguments
//
//    writable   - maya asks for an iterator that can set points if this is true
//    forReadOnly - maya asking for an iterator for querying only
//
{
	return true;
}

void WmFigaroRodShape::transformUsing( const MMatrix& matrix, const MObjectArray &componentList )
{
    // let the other version of transformUsing do the work
    transformUsing( matrix, componentList, MPxSurfaceShape::kNoPointCaching, NULL );
}

void WmFigaroRodShape::updatePointIfNotStretching( MPointArray& io_controlPoints, 
        const int i_index, const MMatrix& i_matrix )
{
    int edgeIndex = i_index - 1;
    double currLength = ( io_controlPoints[ i_index ] - io_controlPoints[ i_index - 1 ] ).length();

    if ( currLength > m_edgeLengths[ edgeIndex ] )
    {
  //      cerr << "this edge is > rest length, correcting!\n";
        MVector edge = io_controlPoints[ i_index ] - io_controlPoints[ i_index - 1 ];
        edge.normalize();
        edge *= m_edgeLengths[ edgeIndex ];
        io_controlPoints[ i_index ] = io_controlPoints[ i_index - 1 ] + edge;
    }

    MPoint oldPoint = io_controlPoints[ i_index ];
    MPoint newPoint = io_controlPoints[ i_index ] * i_matrix;

    // Check if this movement stretch the rod
    
    if ( edgeIndex > 0 )
    {    
        double newLength = ( io_controlPoints[ i_index ] - io_controlPoints[ i_index - 1 ] ).length();
        
        double lengthPlus10Percent = m_edgeLengths[ edgeIndex ] + 0.0001 * m_edgeLengths[ edgeIndex ];
    //    cerr << "end Length = " << newLength << ", old length was " << m_edgeLengths[ edgeIndex ] << ", .001% = " << lengthPlus10Percent << endl;
    
        if ( newLength < lengthPlus10Percent )
        {
            io_controlPoints[ i_index ] = newPoint;
        }
    }
}


void WmFigaroRodShape::transformUsing(const MMatrix& mat, const MObjectArray& componentList,
                     MPxSurfaceShape::MVertexCachingMode cachingMode, MPointArray* pointCache )
{
    MStatus stat;
    MVectorArray& controlPointsVec = *( getControlPoints() );
    MVectorArray unmovedControlPointsVec = controlPointsVec;
    MPointArray controlPoints;
    controlPoints.setLength( controlPointsVec.length() );
    
    // Turn the vectors into Points or the matrix multiplication doesn't work correctly below
    for ( int p=0; p<controlPointsVec.length(); ++p )
    {
        controlPoints[ p ] = controlPointsVec[ p ];
    }

    bool savePoints    = (cachingMode == MPxSurfaceShape::kSavePoints);
    unsigned int i=0,j=0;
    unsigned int len = componentList.length();
    
    if ( cachingMode == MPxSurfaceShape::kRestorePoints ) 
    {
        // restore the points based on the data provided in the pointCache attribute
        //
        unsigned int cacheLen = pointCache->length();
        if (len > 0) 
        {
            // traverse the component list
            //
            for ( i = 0; i < len && j < cacheLen; i++ )
            {
                MObject comp = componentList[i];
                MFnSingleIndexedComponent fnComp( comp );
                int elemCount = fnComp.elementCount();
                for ( int idx=0; idx<elemCount && j < cacheLen; idx++, ++j ) 
                {
                    int elemIndex = fnComp.element( idx );
                    controlPoints[ elemIndex ] = (*pointCache)[ j ];
                }
            }
        } 
        else
        {
            // if the component list is of zero-length, it indicates that we
            // should transform the entire surface
            //
            len = controlPoints.length();
            for ( unsigned int idx = 0; idx < len && j < cacheLen; ++idx, ++j ) 
            {
                controlPoints[ idx ] = (*pointCache)[ j ];
            }
        }
    }
    else 
    {    
        // Transform the surface vertices with the matrix.
        // If savePoints is true, save the points to the pointCache.
        //
        if (len > 0) 
        {
            // Traverse the componentList 
            //
            for ( i=0; i<len; i++ )
            {
                MObject comp = componentList[i];
                MFnSingleIndexedComponent fnComp( comp );
                int elemCount = fnComp.elementCount();

                if (savePoints && 0 == i) 
                {
                    pointCache->setSizeIncrement( elemCount );
                }
                for ( int idx=0; idx<elemCount; idx++ )
                {
                    int elemIndex = fnComp.element( idx );
                    if ( savePoints ) {
                        pointCache->append( controlPoints[ elemIndex ] );
                    }
                    
                    updatePointIfNotStretching( controlPoints, elemIndex, mat );
                  //  controlPoints[ elemIndex ] *= mat;                                        
                }
            }
        }
        else
        {
            // If the component list is of zero-length, it indicates that we
            // should transform the entire surface
            //
            len = controlPoints.length();
            if ( savePoints ) 
            {
                pointCache->setSizeIncrement( len );
            }
            for ( unsigned int idx = 0; idx < len; ++idx ) 
            {
                if ( savePoints ) 
                {
                    pointCache->append( controlPoints[idx]) ;
                }
                updatePointIfNotStretching( controlPoints, idx, mat );
               // controlPoints[ idx ] *= mat;
            }
        }
    }

    FixedRodVertexMap fixedRodVertexMap;
    FixedVertexMap fixedVertexMap;

    bool haveFixedVertices = false;        
    for ( int p=0; p<controlPointsVec.length(); ++p )
    {
        controlPointsVec[ p ] = controlPoints[ p ];

        if ( controlPointsVec[ p ] != unmovedControlPointsVec[ p ] )
        {
            fixedVertexMap[ p ] = Vec3d( controlPointsVec[ p ][ 0 ],
                                         controlPointsVec[ p ][ 1 ],
                                         controlPointsVec[ p ][ 2 ] );            
            
            haveFixedVertices = true;
        }
    }    

  /*  for ( FixedVertexMap::iterator it=fixedVertexMap.begin(); it!=fixedVertexMap.end(); ++it )
    {
        // The rod seems to get really twisted if we only lock single vertices,
        // lock edges next to each fixed vertex to stop it twisting up nastily
        int fixedVertex = it->first;
        int nextFixed = fixedVertex + 1;
        if ( nextFixed >= controlPointsVec.length() )
        {
            nextFixed = fixedVertex - 1; 
        }
        if ( nextFixed < 0 )
        {
            continue;
        }

        // Check if nextFixed is already fixed (ie moved by the user )
        if ( controlPointsVec[ nextFixed ] == unmovedControlPointsVec[ nextFixed ] )
        {
            // We need to fix this one as the user didn't move it.
            // Since vertex nextFixed is possible not getting moved by the user then we have to
            // move it so it keeps up
            MVector unmovedDelta = unmovedControlPointsVec[ nextFixed ] - unmovedControlPointsVec[ fixedVertex ];
            controlPointsVec[ nextFixed ] = controlPointsVec[ fixedVertex ] + unmovedDelta;
        
            fixedVertexMap[ nextFixed ] = Vec3d( controlPointsVec[ nextFixed ][ 0 ],
                                                 controlPointsVec[ nextFixed ][ 1 ],
                                                 controlPointsVec[ nextFixed ][ 2 ] );
        }
    }
*/
    fixedRodVertexMap[ 0 ] = fixedVertexMap;

    if ( haveFixedVertices )
    {
        doSolverIterations( fixedRodVertexMap );
    }
    else
    {
        // If the user hasn't fixed any vertices then assume they're done moving and make it
        // a plastic deformation
        //m_rodGroup.elasticRod( 0 )->updateReferenceProperties();
      //  cerr << "creating new rod!\n";
      //  resetSimulation();
    }

    updateControlPointsFromRod();

    //stat = setControlPoints( &controlPointsVec );
    //CHECK_MSTATUS( stat );

    // Retrieve the value of the cached surface attribute.
    // We will set the new geometry data into the cached surface attribute
    //
    // Access the datablock directly. This code has to be efficient
    // and so we bypass the compute mechanism completely.
    // NOTE: In general we should always go though compute for getting
    // and setting attributes.
    //
    MDataBlock datablock = forceCache();

    /* We do not currently support history on the shape */
    /*
    
    MDataHandle cachedHandle = datablock.outputValue( cachedSurface, &stat );
    MCHECKERRORNORET( stat, "computeInputSurface error getting cachedSurface")
    apiMeshData* cached = (apiMeshData*) cachedHandle.asPluginData();

    MDataHandle dHandle = datablock.outputValue( mControlPoints, &stat );
    MCHECKERRORNORET( stat, "transformUsing get dHandle" )      
    
    // If there is history then calculate the tweaks necessary for
    // setting the final positions of the vertices.
    // 
    if ( hasHistory() && (NULL != cached) ) {
        // Since the shape has history, we need to store the tweaks (deltas)
        // between the input shape and the tweaked shape in the control points
        // attribute.
        //
        stat = buildControlPoints( datablock, geomPtr->vertices.length() );
        MCHECKERRORNORET( stat, "transformUsing buildControlPoints" )

        MArrayDataHandle cpHandle( dHandle, &stat );
        MCHECKERRORNORET( stat, "transformUsing get cpHandle" )

        // Loop through the component list and transform each vertex.
        //
        for ( i=0; i<len; i++ )
        {
            MObject comp = componentList[i];
            MFnSingleIndexedComponent fnComp( comp );
            int elemCount = fnComp.elementCount();
            for ( int idx=0; idx<elemCount; idx++ )
            {
                int elemIndex = fnComp.element( idx );
                cpHandle.jumpToElement( elemIndex );
                MDataHandle pntHandle = cpHandle.outputValue(); 
                double3& pnt = pntHandle.asDouble3();       

                MPoint oldPnt = cached->fGeometry->vertices[elemIndex];
                MPoint newPnt = geomPtr->vertices[elemIndex];
                MPoint offset = newPnt - oldPnt;

                pnt[0] += offset[0];
                pnt[1] += offset[1];
                pnt[2] += offset[2];                
            }
        }
    }

    // Copy outputSurface to cachedSurface
    //
    if ( NULL == cached ) {
        cerr << "NULL cachedSurface data found\n";
    }
    else {
        *(cached->fGeometry) = *geomPtr;
    }

    MPlug pCPs(thisMObject(),mControlPoints);
    pCPs.setValue(dHandle); */

    // Moving vertices will likely change the bounding box.
    //
    computeBoundingBox( datablock );
}

MStatus WmFigaroRodShape::computeBoundingBox( MDataBlock& datablock )
//
// Description
//
//    Use the larges/smallest vertex positions to set the corners 
//    of the bounding box.
//
{
    MStatus stat = MS::kSuccess;

    // Update bounding box
    //
    MDataHandle lowerHandle = datablock.outputValue( i_bboxCorner1 );
    MDataHandle upperHandle = datablock.outputValue( i_bboxCorner2 );
    double3 &lower = lowerHandle.asDouble3();
    double3 &upper = upperHandle.asDouble3();

    MVectorArray& controlPoints = *( getControlPoints() );
    int cnt = controlPoints.length();
    if ( cnt == 0 ) return stat;

    // This clears any old bbox values
    //
    MPoint tmppnt = controlPoints[ 0 ];
    lower[0] = tmppnt[0]; lower[1] = tmppnt[1]; lower[2] = tmppnt[2];
    upper[0] = tmppnt[0]; upper[1] = tmppnt[1]; upper[2] = tmppnt[2];


    for ( int i=0; i<cnt; i++ )
    {
        MPoint pnt = controlPoints[ i ];

        if ( pnt[0] < lower[0] ) lower[0] = pnt[0];
        if ( pnt[1] < lower[1] ) lower[1] = pnt[1];
        if ( pnt[2] < lower[2] ) lower[2] = pnt[2];
        if ( pnt[0] > upper[0] ) upper[0] = pnt[0];
        if ( pnt[1] > upper[1] ) upper[1] = pnt[1];
        if ( pnt[2] > upper[2] ) upper[2] = pnt[2];
    }
    
    lowerHandle.setClean();
    upperHandle.setClean();

    // Signal that the bounding box has changed.
    //
    childChanged( MPxSurfaceShape::kBoundingBoxChanged );

    return stat;
}




























