#include "WmFigControllerContextCommand.hh"

WmFigControllerContextCommand::WmFigControllerContextCommand() 
{
}

void* WmFigControllerContextCommand::creator() 
{
	return new WmFigControllerContextCommand;
}

MPxContext* WmFigControllerContextCommand::makeObj() 
{ 
    WmFigControllerContext* contextObject = new WmFigControllerContext();
    WmFigControllerContext::activeContext = contextObject;
    return contextObject;
}

MStatus WmFigControllerContextCommand::appendSyntax() 
{
    MSyntax stx( syntax() );

    // No options so far...
	
    return MS::kSuccess;
}

MStatus WmFigControllerContextCommand::doEditFlags()
{
    MArgParser prsr( parser() );

    // no flags to set yet
    
    return MS::kSuccess;
}

MStatus WmFigControllerContextCommand::doQueryFlags() 
{
    MArgParser prsr( parser() );

    // no flags to query yet
    
    return MS::kSuccess;
}
