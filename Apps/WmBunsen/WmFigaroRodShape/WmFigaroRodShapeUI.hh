
///////////////////////////////////////////////////////////////////////////////
//
// WmFigaroRodShapeUI.h
//
// Encapsulates the UI portion of a user defined shape. All of the
// drawing and selection code goes here.
//
////////////////////////////////////////////////////////////////////////////////

#include <maya/MPxSurfaceShapeUI.h> 
#include "WmFigaroRodShape.hh"

class WmFigaroRodShapeUI : public MPxSurfaceShapeUI
{
public:
	WmFigaroRodShapeUI();
	virtual ~WmFigaroRodShapeUI(); 

	/////////////////////////////////////////////////////////////////////
	//
	// Overrides
	//
	/////////////////////////////////////////////////////////////////////

	// Puts draw request on the draw queue
	//
	virtual void	getDrawRequests( const MDrawInfo & info,
									 bool objectAndActiveOnly,
									 MDrawRequestQueue & requests );

	// Main draw routine. Gets called by maya with draw requests.
	//
	virtual void	draw( const MDrawRequest & request,
						  M3dView & view ) const;

	// Main selection routine
	//
	virtual bool	select( MSelectInfo &selectInfo,
							MSelectionList &selectionList,
							MPointArray &worldSpaceSelectPts ) const;

	/////////////////////////////////////////////////////////////////////
	//
	// Helper routines
	//
	/////////////////////////////////////////////////////////////////////

	void	drawVertices( const MDrawRequest & request, M3dView & view ) const;

	bool 	selectVertices( MSelectInfo &selectInfo,
		   		MSelectionList &selectionList,
				MPointArray &worldSpaceSelectPts ) const;

	static  void *      creator();

private:
	// Draw Tokens
	//
	enum {
		kDrawVertices, // component token
		kDrawWireframe,
		kDrawWireframeOnShaded,
		kDrawSmoothShaded,
		kDrawFlatShaded,
		kLastToken
	};
};
