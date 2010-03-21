#include "WmFigSelectionContextCommand.hh"

WmFigSelectionContextCommand::WmFigSelectionContextCommand() 
{
}

void* WmFigSelectionContextCommand::creator() 
{
	return new WmFigSelectionContextCommand;
}

MPxContext* WmFigSelectionContextCommand::makeObj() 
{ 
    WmFigSelectionContext* contextObject = new WmFigSelectionContext();
    WmFigSelectionContext::activeContext = contextObject;
    return contextObject;
}

MStatus WmFigSelectionContextCommand::appendSyntax() 
{
    MSyntax stx( syntax() );

    // No options so far...
	
    return MS::kSuccess;
}

MStatus WmFigSelectionContextCommand::doEditFlags()
{
    MArgParser prsr( parser() );

    // no flags to set yet
    
    return MS::kSuccess;
}

MStatus WmFigSelectionContextCommand::doQueryFlags() 
{
    MArgParser prsr( parser() );

    // no flags to query yet
    
    return MS::kSuccess;
}
