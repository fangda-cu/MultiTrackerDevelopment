#ifndef WMSWEENEYUTILS_HH_
#define WMSWEENEYUTILS_HH_

#include "WmSweeneyNode.hh"

// This is a generic place for utility functions that are needed across files and do not really
// live in any one class.

namespace sweeney {
namespace utils {

MStatus findSelectedSweeneyNodeAndRodManager( WmSweeneyNode* o_wmSweeneyNode, WmSweeneyRodManager* o_wmSweenyRodManager );
MStatus findASelectedNodeByTypeName( MString& i_typeName, MObject* o_selectedNodeObject );

}
}

#endif
