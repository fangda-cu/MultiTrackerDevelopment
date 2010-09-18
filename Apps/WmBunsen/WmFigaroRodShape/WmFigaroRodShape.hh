#ifndef _WmFigaroRodShape
#define _WmFigaroRodShape

///////////////////////////////////////////////////////////////////////////////
//
// WmFigaroRodShape.h
//
// Implements a new type of shape node in maya called WmFigaroRodShape.
//
// To use it
//
// loadPlugin WmFigaroRodShape
// string $node = `createNode WmFigaroRodShape`; // You'll see nothing.
//
//
// // Now add some CVs, one
// string $attr = $node + ".controlPoints[0]";
// setAttr $attr 2 2 2; 					// Now you'll have something on screen, in wireframe mode
//
//
// // or a bunch
// int $idx = 0;
// for ( $i=0; $i<100; $i++)
// {
//    for ( $j=0; $j<100; $j++)
//    {
//        string $attr = $node + ".controlPoints[ " + $idx + "]";
//        setAttr $attr $i $j 3;
//        $idx++;
//    }
// }
//
//
// INPUTS
//     NONE
//
// OUTPUTS
// 	   NONE
//
////////////////////////////////////////////////////////////////////////////////

#include <maya/MTypeId.h> 
#include <maya/MPxComponentShape.h> 
#include "../Beaker.hh"
#include "../WmFigRodNode.hh"

class MPxGeometryIterator;

class WmFigaroRodShape : public MPxComponentShape
{
public:
	WmFigaroRodShape();
	virtual ~WmFigaroRodShape(); 

    virtual MStatus compute( const MPlug& i_plug, MDataBlock& i_dataBlock );

	// Associates a user defined iterator with the shape (components)
	//
    virtual void postConstructor();
	virtual MPxGeometryIterator* geometryIteratorSetup( MObjectArray&, MObject&, bool forReadOnly = false );
	virtual bool acceptsGeometryIterator( bool  writeable=true );
	virtual bool acceptsGeometryIterator( MObject&, bool writeable=true, bool forReadOnly = false );
    virtual void transformUsing( const MMatrix& matrix, const MObjectArray &componentList ); 
    virtual void transformUsing( const MMatrix& mat, const MObjectArray& componentList,
                                 MPxSurfaceShape::MVertexCachingMode cachingMode, MPointArray* pointCache );
    bool virtual setInternalValueInContext( const MPlug& i_plug, const MDataHandle& handle, MDGContext& i_context );

	static  void *creator();
	static  MStatus initialize();

    static MObject i_bboxCorner1;
    static MObject i_bboxCorner2;

    static MObject oa_cv;
    static MObject oa_cvX;
    static MObject oa_cvY;
    static MObject oa_cvZ;
    
	static	MTypeId	id;
    static  MString typeName;

    static MObject ia_triggerSolve;

    static MObject ia_doPlasticDeformations;
    static MObject ia_constraintStrength;

    static MObject ia_percentLocked;

    static MObject ia_cylinderDrawScale;

    static MObject ia_youngsModulus;
    static MObject ia_shearModulus;
    static MObject ia_viscosity;
    static MObject ia_density;
    static MObject ia_radiusA;
    static MObject ia_radiusB;

    static MObject ca_sync;

   /* static MObject ia_inCV;
    static MObject ia_inCVX;
    static MObject ia_inCVY;
    static MObject ia_inCVZ;*/

    void drawRod();
    void resetSimulation( MVectorArray* i_controlPoints = NULL );

    int numberOfLockedvertices()
    {
        return m_numberOfLockedVertices;
    }
    
private:
    MStatus computeBoundingBox( MDataBlock& datablock );
    MStatus buildControlPoints( MVectorArray& i_controlPoints );
    void getRodVertices( MVectorArray& o_controlPoints );
    void updateControlPointsFromRod();
    void initialiseRod( MVectorArray* i_controlPoints );
    void updatePointIfNotStretching( MPointArray& io_controlPoints, const size_t i_index, const MMatrix& i_matrix );
    bool doSolverIterations( FixedRodVertexMap& i_fixedRodVertexMap );

    void solve( FixedRodVertexMap& i_fixedRodVertexMap );
    
    WmFigRodGroup m_rodGroup;
    Beaker m_beaker;
    
    std::vector< double > m_edgeLengths;
    bool m_initialisedRod;
    double m_constraintStrength;

    FixedVertexMap m_fixedVertexMap;    
    LockedVertexMap m_lockedVertexMap;    

    bool m_drawRodTube;

    int m_numberOfLockedVertices;
    double m_cylinderDrawScale;
    double m_youngsModulus;
    double m_shearModulus;
    double m_viscosity;
    double m_density;
    double m_radiusA;
    double m_radiusB;
};


//////////////////////////////////////////////////////////////////////
//
// Error checking
//
//    MCHECKERROR       - check the status and print the given error message
//    MCHECKERRORNORET  - same as above but does not return
//
//////////////////////////////////////////////////////////////////////

#define MCHECKERROR(STAT,MSG)       \
    if ( MS::kSuccess != STAT ) {   \
        cerr << MSG << endl;        \
            return MS::kFailure;    \
    }

#define MCHECKERRORNORET(STAT,MSG)  \
    if ( MS::kSuccess != STAT ) {   \
        cerr << MSG << endl;        \
    }

//////////////////////////////////////////////////////////////////////
//
// Attribute creation
//
//       MAKE_TYPED_ATTR   - creates and adds a typed attribute
//       MAKE_NUMERIC_ATTR - creates and adds a numeric attribute
//       ADD_ATTRIBUTE     - adds the given attribute
//       ATTRIBUTE_AFFECTS - calls attributeAffects
//
//////////////////////////////////////////////////////////////////////

#define MAKE_TYPED_ATTR( NAME, LONGNAME, SHORTNAME, TYPE, DEFAULT )         \
                                                                            \
    MStatus NAME##_stat;                                                    \
    MFnTypedAttribute NAME##_fn;                                            \
    NAME = NAME##_fn.create( LONGNAME, SHORTNAME, TYPE, DEFAULT );          \
    NAME##_fn.setHidden( true );                                            \
    NAME##_stat = addAttribute( NAME );                                     \
    MCHECKERROR(NAME##_stat, "addAttribute error");

#define MAKE_NUMERIC_ATTR( NAME, LONGNAME, SHORTNAME, TYPE, DEFAULT,        \
                            ARRAY, BUILDER, KEYABLE )                       \
                                                                            \
    MStatus NAME##_stat;                                                    \
    MFnNumericAttribute NAME##_fn;                                          \
    NAME = NAME##_fn.create( LONGNAME, SHORTNAME, TYPE, DEFAULT );          \
    MCHECKERROR(NAME##_stat, "numeric attr create error");                  \
    NAME##_fn.setArray( ARRAY );                                            \
    NAME##_fn.setUsesArrayDataBuilder( BUILDER );                           \
    NAME##_fn.setHidden( ARRAY );                                           \
    NAME##_fn.setKeyable( KEYABLE );                                        \
    NAME##_stat = addAttribute( NAME );                                     \
    MCHECKERROR(NAME##_stat, "addAttribute error");

#define ADD_ATTRIBUTE( ATTR )                                               \
    MStatus ATTR##_stat;                                                    \
    ATTR##_stat = addAttribute( ATTR );                                     \
    MCHECKERROR( ATTR##_stat, "addAttribute: ATTR" )

#define ATTRIBUTE_AFFECTS( IN, OUT )                                        \
    MStatus IN##OUT##_stat;                                                 \
    IN##OUT##_stat = attributeAffects( IN, OUT );                           \
    MCHECKERROR(IN##OUT##_stat,"attributeAffects:" #IN "->" #OUT);


#endif /* _WmFigaroRodShape */
