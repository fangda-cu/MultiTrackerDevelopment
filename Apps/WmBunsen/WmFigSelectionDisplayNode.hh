#ifndef WMFIGSELECTIONDISPLAYNODE_H
#define WMFIGSELECTIONDISPLAYNODE_H

#include "WmFigRodComponentList.hh"
#include <maya/MPxLocatorNode.h>
#include <maya/MColorArray.h>
#include <maya/MIntArray.h>
#include <map>

// Notes: David Gould
//
// Why have this class?
// ---------------------
// The rationale for this class is that the rods's own draw function won't allow for the display of the selected
// vertices, etc. That draw function is inside of BASim which would therefore require changing that function. BASim
// is used by Columbia University so I'm hesitant to make changes in there.
//
//
class WmFigSelectionDisplayNode : public MPxLocatorNode
{
public:
	WmFigSelectionDisplayNode();

    virtual void draw( M3dView &view, const MDagPath &path, M3dView::DisplayStyle style, M3dView::DisplayStatus status );
    virtual bool isBounded() const;

    static void *creator();
    static MStatus initialize();

    static const MTypeId TypeId;
    static const MString TypeName;

    WmFigRodComponentList selection;
};

#endif
