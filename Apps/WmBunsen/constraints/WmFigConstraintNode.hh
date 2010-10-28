/*
 * WmFigConstraintNode.hh
 *
 *  Created on: Oct 12, 2010
 *      Author: dgould
 */

#ifndef WMFIGCONSTRAINTNODE_H
#define WMFIGCONSTRAINTNODE_HH

#include <maya/MTypeId.h>
#include <maya/MPxLocatorNode.h>
#include "../../..//BASim/src/Physics/ElasticRods/RodPenaltyForce.hh"
#include <list>

class WmFigConstraintNode : public MPxLocatorNode
{
public:
	WmFigConstraintNode();
	~WmFigConstraintNode();

	virtual MStatus compute( const MPlug &plug, MDataBlock &data );
	static  void *creator();
	static  MStatus initialize();

    virtual void draw( M3dView &view, const MDagPath &path, M3dView::DisplayStyle style, M3dView::DisplayStatus status );
    virtual bool isBounded() const;

	static const MTypeId TypeId;
    static const MString TypeName;

    // Basic attributes
    //
    static MObject enable;
    static MObject constraintType;
    static MObject stiffness;
    static MObject targetWorldPosition;
    static MObject rodVertices;
    static MObject figRodNodeMsg;

    std::list<BASim::RodVertexConstraint *> rodVertexConstraints;
};

#endif
