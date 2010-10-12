/*
 * WmFigConstraintNode.hh
 *
 *  Created on: Oct 12, 2010
 *      Author: dgould
 */

#ifndef WMFIGCONSTRAINTNODE_H
#define WMFIGCONSTRAINTNODE_HH

#include <maya/MTypeId.h>
#include <maya/MPxNode.h>
#include "../../..//BASim/src/Physics/ElasticRods/RodPenaltyForce.hh"
#include <list>

class WmFigConstraintNode : public MPxNode
{
public:
	WmFigConstraintNode();
	~WmFigConstraintNode();

	virtual MStatus compute( const MPlug &plug, MDataBlock &data );
	static  void *creator();
	static  MStatus initialize();

	static const MTypeId TypeId;
    static const MString TypeName;

    // Basic attributes
    //
    static MObject enable;
    static MObject constraintType;
    static MObject stiffness;
    static MObject worldPosition;
    static MObject rodId;
    static MObject vertexId;
    static MObject figRodNodeMsg;

    std::list<BASim::RodVertexConstraint *> rodVertexConstraints;

#if 0
    // Basic attributes
    //
	static MObject applyDeformer;
	static MObject envelope;
	static MObject deformUsingParameterWeights;
	static MObject parameterValues;
	static MObject outWPSDDataRef;

	// During training attributes
	//
	static MObject minDistanceThreshold;

	// Display attributes
	//
    //static MObject displayBlending;
    static MObject displayParameters;
    static MObject displayParametersNonConstant;
    static MObject displayParametersDetails;
    static MObject displayFlagVerticesForExactMatch;

	static MObject displayTraining;
	static MObject displayWhichTrainingSample;
	static MObject displayFontSizeChoice;

	// Internal/compute attributes
	//
	static MObject _numTrainingSamples;
	static MObject _parameterWeightsPerVertex;

	// Private attributes
	//
	static MObject __wpsdData;

	static void DeleteCallback(MObject &node, void *clientData);
	static void AttributeChangedCallback(MNodeMessage::AttributeMessage msg, MPlug  &plug, MPlug  &otherPlug, void *clientData);
#endif
};

#endif
