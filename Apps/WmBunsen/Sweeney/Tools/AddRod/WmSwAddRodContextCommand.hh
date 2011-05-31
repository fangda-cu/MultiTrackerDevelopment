/////////////////////////////////////////////
// Command to create contexts
//////////////////////////////////////////////

#ifndef WMSWADDRODCONTEXTCOMMAND_HH_
#define WMSWADDRODCONTEXTCOMMAND_HH_

#include "WmSwAddRodContext.hh"

#include <maya/MPxContextCommand.h>
#include <maya/MPxContext.h>
#include <maya/MEvent.h>

namespace sweeney {

class WmSwAddRodContextCommand : public MPxContextCommand
{
public:

    WmSwAddRodContextCommand();
    virtual MPxContext*	makeObj();
    static	void*		creator();
	
    virtual MStatus appendSyntax();
    virtual MStatus doEditFlags();
    virtual MStatus doQueryFlags();
};

}

#endif 
