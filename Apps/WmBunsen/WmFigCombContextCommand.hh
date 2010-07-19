/////////////////////////////////////////////
// Command to create contexts
//////////////////////////////////////////////

#ifndef _WMFIGCOMBCONTEXTCOMMAND_HH_
#define _WMFIGCOMBCONTEXTCOMMAND_HH_

#include "WmFigCombContext.hh"

#include <maya/MPxContextCommand.h>
#include <maya/MPxContext.h>
#include <maya/MEvent.h>


class WmFigCombContextCommand : public MPxContextCommand
{
public:

    WmFigCombContextCommand();
    virtual MPxContext*	makeObj();
    static	void*		creator();
	
    virtual MStatus appendSyntax();
    virtual MStatus doEditFlags();
    virtual MStatus doQueryFlags();
};

#endif 
