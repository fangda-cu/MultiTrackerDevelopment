/////////////////////////////////////////////
// Command to create contexts
//////////////////////////////////////////////

#ifndef _WMFIGCONTROLLERCONTEXTCOMMAND_HH_
#define _WMFIGCONTROLLERCONTEXTCOMMAND_HH_

#include "WmFigControllerContext.hh"

#include <maya/MPxContextCommand.h>
#include <maya/MPxContext.h>
#include <maya/MEvent.h>


class WmFigControllerContextCommand : public MPxContextCommand
{
public:

    WmFigControllerContextCommand();
    virtual MPxContext*	makeObj();
    static	void*		creator();
	
    virtual MStatus appendSyntax();
    virtual MStatus doEditFlags();
    virtual MStatus doQueryFlags();
};

#endif 
