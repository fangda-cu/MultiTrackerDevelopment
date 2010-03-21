/////////////////////////////////////////////
// Command to create contexts
//////////////////////////////////////////////

#ifndef _WMFIGSELECTIONCONTEXTCOMMAND_HH_
#define _WMFIGSELECTIONCONTEXTCOMMAND_HH_

#include "WmFigSelectionContext.hh"

#include <maya/MPxContextCommand.h>
#include <maya/MPxContext.h>
#include <maya/MEvent.h>


class WmFigSelectionContextCommand : public MPxContextCommand
{
public:

    WmFigSelectionContextCommand();
    virtual MPxContext*	makeObj();
    static	void*		creator();
	
    virtual MStatus appendSyntax();
    virtual MStatus doEditFlags();
    virtual MStatus doQueryFlags();
};

#endif 
