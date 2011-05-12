#ifndef WMSWEENEYSHAPE_HH_
#define WMSWEENEYSHAPE_HH_

///////////////////////////////////////////////////////////////////////////////
//
// WmSweeneyShape.h
//
// Implements a new type of shape node in maya called WmSweeneyShape.
//
// To use it
//
// loadPlugin WmSweeneyShape
// string $node = `createNode WmSweeneyShape`; // You'll see nothing.
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

#include "WmSweeneyRodManager.hh"

class MPxGeometryIterator;

class WmSweeneyShape : public MPxComponentShape
{
public:
	WmSweeneyShape();
	virtual ~WmSweeneyShape(); 

	// Associates a user defined iterator with the shape (components)
	//
	virtual MPxGeometryIterator*	geometryIteratorSetup( MObjectArray&, MObject&, bool forReadOnly = false );
	virtual bool                    acceptsGeometryIterator( bool  writeable=true );
	virtual bool                    acceptsGeometryIterator( MObject&, bool writeable=true, bool forReadOnly = false);

	static  void *          creator();
	static  MStatus         initialize();

	static MTypeId id;
    static MString typeName;
    
private:
    void initialiseRodFromBarberShopInput();

    WmSweeneyRodManager rodManager;
};

#endif /* _WmSweeneyShape */
