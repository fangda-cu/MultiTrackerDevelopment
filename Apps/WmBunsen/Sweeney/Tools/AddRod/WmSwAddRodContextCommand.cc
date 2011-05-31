#include "WmSwAddRodContextCommand.hh"

namespace sweeney {

WmSwAddRodContextCommand::WmSwAddRodContextCommand() 
{
}

void* WmSwAddRodContextCommand::creator() 
{
	return new WmSwAddRodContextCommand;
}

MPxContext* WmSwAddRodContextCommand::makeObj() 
{ 
    WmSwAddRodContext* contextObject = new WmSwAddRodContext();
    WmSwAddRodContext::activeContext = contextObject;
    return contextObject;
}

MStatus WmSwAddRodContextCommand::appendSyntax() 
{
    MSyntax stx( syntax() );

    // No options so far...
	
    return MS::kSuccess;
}

MStatus WmSwAddRodContextCommand::doEditFlags()
{
    MArgParser prsr( parser() );

    // no flags to set yet
    
    return MS::kSuccess;
}

MStatus WmSwAddRodContextCommand::doQueryFlags() 
{
    MArgParser prsr( parser() );

    // no flags to query yet
    
    return MS::kSuccess;
}

}
