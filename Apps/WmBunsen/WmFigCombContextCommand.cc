#include "WmFigCombContextCommand.hh"

WmFigCombContextCommand::WmFigCombContextCommand() 
{
}

void* WmFigCombContextCommand::creator() 
{
	return new WmFigCombContextCommand;
}

MPxContext* WmFigCombContextCommand::makeObj() 
{ 
    WmFigCombContext* contextObject = new WmFigCombContext();
    WmFigCombContext::activeContext = contextObject;
    return contextObject;
}

MStatus WmFigCombContextCommand::appendSyntax() 
{
    MSyntax stx( syntax() );

    // No options so far...
	
    return MS::kSuccess;
}

MStatus WmFigCombContextCommand::doEditFlags()
{
    MArgParser prsr( parser() );

    // no flags to set yet
    
    return MS::kSuccess;
}

MStatus WmFigCombContextCommand::doQueryFlags() 
{
    MArgParser prsr( parser() );

    // no flags to query yet
    
    return MS::kSuccess;
}
