/*
 * WmFigConstraintNode.cc
 *
 *  Created on: Oct 12, 2010
 *      Author: dgould
 */
#include "WmFigConstraintNode.hh"
#include <maya/MFnEnumAttribute.h>
#include <maya/MFnMatrixAttribute.h>
#include <maya/MFnMessageAttribute.h>
#include <maya/MFnNumericAttribute.h>

const MTypeId WmFigConstraintNode::TypeId( 0x001156, 0x60 ); // Registered in /vol/weta/src/Wcommon/MTypeId.txt
const MString WmFigConstraintNode::TypeName( "wmFigConstraintNode" );

// Basic attributes
//
MObject WmFigConstraintNode::enable;
MObject WmFigConstraintNode::constraintType;
MObject WmFigConstraintNode::stiffness;
MObject WmFigConstraintNode::worldPosition;
MObject WmFigConstraintNode::rodId;
MObject WmFigConstraintNode::vertexId;

MObject WmFigConstraintNode::figRodNodeMsg;
//MObject WmFigConstraintNode::inWorldMatrix;

WmFigConstraintNode::WmFigConstraintNode()
{

}

WmFigConstraintNode::~WmFigConstraintNode()
{

}

MStatus WmFigConstraintNode::compute( const MPlug &plug, MDataBlock &data )
{
	MStatus stat;
	return MS::kUnknownParameter; // Haven't computed the plug, pass it back up the class hierarchy

#if 0
	MDataHandle stateHnd = data.inputValue( state );
	int state = stateHnd.asInt();
	bool hasNoEffect = (state == 1); // No Effect/Pass through

	if( plug == outWPSDDataRef ) {
	    // Even though it is the outWPSSDataRef that is being requested, ensure that the outWPSDData MPxData that it will reference
		// has been created. It will be empty by default.
	    //
		MDataHandle outWPSDDataHnd = data.inputValue( __wpsdData );
		WmWeightedPSDData *wpsdData = (WmWeightedPSDData *)outWPSDDataHnd.asPluginData();
		if( !wpsdData ) {
			MFnPluginData pluginDataFn;
			MObject wpsdDataObj = pluginDataFn.create( WmWeightedPSDData::TypeId );
			outWPSDDataHnd.set( wpsdDataObj );
			outWPSDDataHnd.setClean();

			outWPSDDataHnd = data.inputValue( __wpsdData );
			wpsdData = (WmWeightedPSDData *)outWPSDDataHnd.asPluginData();
		}

		//cout << "------- computing outWPSDDataRef" << std::endl;
		MDataHandle outWPSDDataRefHnd = data.outputValue( outWPSDDataRef );

		//cout << "------- creating instance of WmWeightedPSDDataRef " << std::endl;
		MFnPluginData pluginDataFn;
		MObject wpsdDataRefObj = pluginDataFn.create( WmWeightedPSDDataRef::TypeId, &stat );
		if( !stat ) {
			//cout << "------- unabled to create instance of WmWeightedPSDDataRef" << std::endl;
			MGlobal::displayWarning( MString( "Unable to create instance of ") + WmWeightedPSDDataRef::TypeName );
			return MS::kFailure;
		}

		WmWeightedPSDDataRef *wpsdDataRefPtr = (WmWeightedPSDDataRef *) pluginDataFn.data();
		wpsdDataRefPtr->wpsdDataPtr_weak = wpsdData; // Set the pointer to the the WmWeightedPSDData within this node

		outWPSDDataRefHnd.set( wpsdDataRefObj );
		outWPSDDataRefHnd.setClean();
	} else {
		// Compute the outMesh
		//
		if( plug == outMesh ) {
			//cout << "computing outMesh" << std::endl;

			MDataHandle inMeshHnd = data.inputValue( inMesh );
			MDataHandle outMeshHnd = data.outputValue( outMesh );

			outMeshHnd.copy( inMeshHnd );
			MObject inMesh = inMeshHnd.asMesh();
			MObject newMesh = outMeshHnd.asMesh();

			if( !hasNoEffect ) {
				MDataHandle applyHnd = data.inputValue( applyDeformer );
			    bool applyValue = applyHnd.asBool();

			    MDataHandle envelopeHnd = data.inputValue( envelope );
			    double envelopeValue = envelopeHnd.asDouble();

			    MDataHandle displayTrainingHnd = data.inputValue( displayTraining );
			    bool displayTraining = displayTrainingHnd.asBool();

			    if( displayTraining ) {
				    MDataHandle displayWhichTrainingSampleHnd = data.inputValue( displayWhichTrainingSample );
				    int tsIndex = displayWhichTrainingSampleHnd.asInt();

				    MObject thisObj = thisMObject();
				    WFnWeightedPSD wpsdFn( thisObj );
				    int numTrainingSamples = wpsdFn.numTrainingSamples();
				    if( numTrainingSamples > 0 ) {
						tsIndex--; // Convert from 1-based index to 0-based index

				    	if( tsIndex < 1 )
							tsIndex = 0;
						else {
							if (tsIndex >= numTrainingSamples )
								tsIndex = numTrainingSamples-1;
						}

						//cout << "updating mesh with training sample mesh " << tsIndex << std::endl;

						// Get the training shape data
						WmWeightedPSDData::TrainingSample &trainingSample = *wpsdFn.getTrainingSample( tsIndex );

						// Set the mesh vertices to the training shape positions
						//
						MFnMesh meshFn( newMesh );
//						MPointArray trainingPoints( trainingSample.positionValues );
//						if( wpsdFn.getPositionsAreFaceRelative() ) {
//							ConvertFaceRelativePositionsToLocalPositions( inMesh, trainingPoints );
//						}

						MPointArray trainingPoints;
						if( wpsdFn.getPositionsAreFaceRelative() ) {
							ConvertFaceRelativePositionsToLocalPositions( inMesh,
																		  trainingSample.sparseVertexIndices,
																		  trainingSample.sparseVertexPositions,
																		  trainingSample._mapVertexIndexToSparseIndex,
																		  trainingPoints );
						} else { // Absolute positions (and not sparse)
							trainingPoints = trainingSample.sparseVertexPositions;
						}

						meshFn.setPoints( trainingPoints, MSpace::kObject );
						meshFn.updateSurface();
					}
			    } else {
			    	if( applyValue ) {
			    		// Get the wpsdData
			    		//
			    		MDataHandle outWPSDDataRefHnd = data.inputValue( outWPSDDataRef );
			    		WmWeightedPSDDataRef &wpsdDataRef = *(WmWeightedPSDDataRef *)outWPSDDataRefHnd.asPluginData();
			    		WmWeightedPSDData &wpsdData = *wpsdDataRef.wpsdDataPtr_weak;

			    		// There has been some training
			    		//
			    		if( wpsdData.trainingSamples.size() ) {
							// Get the input parameters (they are not normalized since they could come from attribute connections
			    			// to say, a rotation which is an angle in degrees)
							//
							MArrayDataHandle parameterValuesHnd = data.inputArrayValue( parameterValues );
							unsigned int nParameterValues = parameterValuesHnd.elementCount();

							MDataHandle parameterValueHnd;
							MDoubleArray parameterValues;
							parameterValues.setLength( nParameterValues );
							unsigned int i;
							for( i=0; i < nParameterValues; i++ ) {
								parameterValueHnd = parameterValuesHnd.inputValue();
								parameterValues[i] = parameterValueHnd.asDouble();
								parameterValuesHnd.next();
							}

//							MString txt( "Current parameters: " );
//							for( i=0; i < parameterValues.length(); i++ )
//								txt = txt + parameterValues[i] + " ";
//							MGlobal::displayInfo( txt );

							// Check that the parameters are within the parameter bounds
							//
							if( wpsdData.areParametersWithinBounds( parameterValues, true ) ) {

								// Normalize parameter values so that they are each in the range [0,1]
								//
								wpsdData.normalizeParameters( parameterValues );

								MDataHandle deformUsingParameterWeightsHnd = data.inputValue( deformUsingParameterWeights );
							    bool useParameterWeights = deformUsingParameterWeightsHnd.asBool();

								MFnMesh meshFn( newMesh );
								unsigned int nVerts = meshFn.numVertices();
								assert( wpsdData.numVertices() != nVerts );

								//cout << "ENVELOPE: " << envelopeValue << std::endl;

								MPointArray newLocalPositions( nVerts );
								//wpsdData.calcWPSDPositions( parameterValues, useParameterWeights, blendedPositions );
								wpsdData.setUseParameterWeights( useParameterWeights );
								//wpsdData.calcWPSDPositions( newMesh, parameterValues, useParameterWeights, newLocalPositions );
								wpsdData.calcWPSDPositions( newMesh, parameterValues, envelopeValue, newLocalPositions );

	//							cout << "--> nVerts " << meshFn.numVertices() << std::endl;
	//							cout << "--> nParameterWeightsPerVertex " << wpsdData.parameterWeightsPerVertex.size() << std::endl;

/*
								// The weighted psd blended the face-relative positions. Convert these to local-space positions.
								//
								if( wpsdData.positionsAreFaceRelative ) {
									ConvertFaceRelativePositionsToLocalPositions( inMesh, blendedPositions );
								}
*/
								meshFn.setPoints( newLocalPositions, MSpace::kObject );
								meshFn.updateSurface();
							}
			    		}
			    	}
			    }
			}

//			cout << "computed mesh " << (computedMesh ? "yes" : "no") << std::endl;

//			// No new mesh created, just pass the inMesh through to the outMesh
//			if( !computedMesh ) {
//				MFnMesh meshFn( inMeshCopy );
//				meshFn.updateSurface();
//			}

			outMeshHnd.set( newMesh );
			outMeshHnd.setClean();
		} else {
			if( plug == _parameterWeightsPerVertex ) {
				cout << "compute _parameterWeightsPerVertex" << std::endl;

				MArrayDataHandle paramWeightsPerVertHnd = data.outputArrayValue( _parameterWeightsPerVertex );

				MObject thisObj = thisMObject();
				WFnWeightedPSD wpsdFn( thisObj );
				WmWeightedPSDData::ParameterWeightsPerVertex *paramWeightsPerVertex = wpsdFn.getParameterWeightsPerVertex();

				const unsigned int nVertices = paramWeightsPerVertex->size();

				MArrayDataBuilder builder = paramWeightsPerVertHnd.builder();

				// Remove the existing elements
				unsigned int i;
				for( i=0; i < nVertices; i++ )
					builder.removeElement(0);

				builder.growArray( nVertices );

				//cout << "num verts " << paramWeightsPerVertex->size() << std::endl;

				unsigned int vi;
				for( vi=0; vi < paramWeightsPerVertex->size(); vi++ ) {
					MDataHandle elemDataHandle = builder.addLast();

					MFnDoubleArrayData weightsArrayFn;
					MObject paramWeightsObj = weightsArrayFn.create( (*paramWeightsPerVertex)[vi] );
					elemDataHandle.set( paramWeightsObj );
				}

				paramWeightsPerVertHnd.set( builder );
				paramWeightsPerVertHnd.setAllClean();
			} else
				return MS::kUnknownParameter; // Haven't computed the plug, pass it back up the class hierarchy
		}
	}
#endif

	data.setClean( plug );
	return stat;
}

void *WmFigConstraintNode::creator()
{
	return new WmFigConstraintNode();
}

MStatus WmFigConstraintNode::initialize()
{
	MStatus stat;

	// Create and initialize custom attributes
	//
	MFnEnumAttribute eAttr;
    //MFnMatrixAttribute mAttr;
    MFnMessageAttribute msgAttr;
    MFnNumericAttribute nAttr;

	enable = nAttr.create( "enable", "e", MFnNumericData::kBoolean, true );

	constraintType = eAttr.create( "constraintType", "ct", 1 );
    eAttr.addField( "Fixed", 0 );
    eAttr.addField( "Rest", 1 );
    //eAttr.addField( "Distance", 2 );

	stiffness = nAttr.create( "stiffness", "stf", MFnNumericData::kDouble, 50.0 );
	nAttr.setKeyable( true );

	worldPosition = nAttr.create( "worldPosition", "wp", MFnNumericData::k3Double );
	nAttr.setDefault( 1.0, 0.0, 0.0 );
	nAttr.setKeyable( true );

    rodId = nAttr.create( "rodId", "ri", MFnNumericData::kInt, -1 );
    vertexId = nAttr.create( "vertexId", "vi", MFnNumericData::kInt, -1 );

    figRodNodeMsg = msgAttr.create( "figRodNodeMsg", "frm" );
    msgAttr.setHidden( true );

    //inWorldMatrix = mAttr.create( "inWorldMatrix", "iwm" );

    addAttribute( enable );
    addAttribute( constraintType );
    addAttribute( stiffness );
    addAttribute( worldPosition );
    addAttribute( rodId );
    addAttribute( vertexId );
    addAttribute( figRodNodeMsg );
    //addAttribute( inWorldMatrix );


#if 0
	MFnNumericAttribute nAttr;
	MFnTypedAttribute tAttr;


	// Basic attributes
	applyDeformer = nAttr.create( "applyDeformer", "a", MFnNumericData::kBoolean, false );

	envelope = nAttr.create( "envelope", "env", MFnNumericData::kDouble, 1.0 );
	nAttr.setMin( 0.0 );
	nAttr.setMax( 1.0 );

	deformUsingParameterWeights = nAttr.create( "deformUsingParameterWeights", "dupw", MFnNumericData::kBoolean, true );
	parameterValues = nAttr.create( "parameterValues", "pvs", MFnNumericData::kDouble, 0.0 );
	nAttr.setArray( true );
	nAttr.setUsesArrayDataBuilder( true );

	outWPSDDataRef = tAttr.create( "outWPSDDataRef", "owpsdref", WmWeightedPSDDataRef::TypeId, MObject::kNullObj, &stat );
	if (!stat) {
		MGlobal::displayWarning( "Unable to create \"outWPSDDataRef\" attribute" );
		return stat;
	}
	tAttr.setStorable(false); // This is just a runtime reference (pointer) so no need to store

	minDistanceThreshold = nAttr.create( "minDistanceThreshold", "mdt", MFnNumericData::kDouble, 0.01 );
	nAttr.setMin( 0.0 );
	nAttr.setMax( 5.0 );

	// Display attributes
	displayParameters = nAttr.create( "displayParameters", "dbp", MFnNumericData::kBoolean, false );
	displayParametersNonConstant = nAttr.create( "displayParametersNonConstant", "dbn", MFnNumericData::kBoolean, true );
	displayFlagVerticesForExactMatch = nAttr.create( "displayFlagVerticesForExactMatch", "dbe", MFnNumericData::kBoolean, true );
	displayParametersDetails = nAttr.create( "displayParametersDetails", "dpr", MFnNumericData::kBoolean, false );
	displayTraining = nAttr.create( "displayTraining", "dtrn", MFnNumericData::kBoolean, false );
	displayWhichTrainingSample = nAttr.create( "displayWhichTrainingSample", "dwts", MFnNumericData::kInt, 0 );
	nAttr.setMin( 1 );
    displayFontSizeChoice = eAttr.create( "displayFontSizeChoice", "dfsc", 0 );
    eAttr.addField( "10 point", 0 );
    eAttr.addField( "12 point", 1 );
    eAttr.addField( "18 point", 2 );

	// Internal attributes
	//
	_numTrainingSamples = nAttr.create( "_numTrainingSamples", "_nts", MFnNumericData::kInt, 0 );
	nAttr.setStorable(false);
	nAttr.setWritable(false);
	nAttr.setInternal(true);
	_parameterWeightsPerVertex = tAttr.create( "_parameterWeightsPerVertex", "_pwpv", MFnData::kDoubleArray );
	tAttr.setArray(true);
	tAttr.setStorable(false);
	tAttr.setWritable(false);
	//tAttr.setInternal(true);
	tAttr.setUsesArrayDataBuilder(true);

	// Private attributes
	//
	__wpsdData = tAttr.create( "__wpsdData", "__wpsd", WmWeightedPSDData::TypeId, MObject::kNullObj, &stat );
	if (!stat) {
		MGlobal::displayWarning( "Unable to create \"__wpsdData\" attribute" );
		return stat;
	}
	tAttr.setHidden(true);
	// Use the outWPSDataRef to make connections to other nodes rather than directly to this
	// attribute otherwise getting the attribute will more likely result in copying of the data
	// which can be heavy.
	//
	tAttr.setConnectable(false);

	// Basic attributes
	//
	addAttribute( applyDeformer );
	addAttribute( envelope );
	addAttribute( deformUsingParameterWeights );
	addAttribute( parameterValues );
	addAttribute( outWPSDDataRef );

	// During training Attributes
	//
	addAttribute( minDistanceThreshold );

	// Display Attributes
	//
	//addAttribute( displayBlending );
	addAttribute( displayParameters );
	addAttribute( displayParametersNonConstant );
	addAttribute( displayFlagVerticesForExactMatch );
	addAttribute( displayParametersDetails );
	addAttribute( displayTraining );
	addAttribute( displayWhichTrainingSample );
	addAttribute( displayFontSizeChoice );

	// Internal/compute attributes
	//
	addAttribute( _numTrainingSamples );
	addAttribute( _parameterWeightsPerVertex );

	// Private attributes
	//
	addAttribute( __wpsdData );

	// N.B. Maya doesn't understand "implicit" affects very well. e.g. A affects B affect C. It runs great the first time but then if A is changed it doesn't update C.
	// To overcome this specify all affects explcitly, i.e. A affect B, B affects C, A affects C
	//
	MObject affected;
	unsigned int i;
	MObjectArray isAffectedList;
	isAffectedList.append( outWPSDDataRef );
	isAffectedList.append( outMesh );
	for( i=0; i < isAffectedList.length(); i++ ) {
		affected = isAffectedList[i];

		attributeAffects( inMesh, affected );
		attributeAffects( applyDeformer, affected );
		attributeAffects( envelope, affected );
		attributeAffects( deformUsingParameterWeights, affected );
		attributeAffects( parameterValues, affected );
	}

	attributeAffects( __wpsdData, outWPSDDataRef );
	attributeAffects( __wpsdData, outMesh );
	attributeAffects( outWPSDDataRef, outMesh );

	attributeAffects( inMesh, _parameterWeightsPerVertex );

	attributeAffects( displayTraining, outMesh );
	attributeAffects( displayWhichTrainingSample, outMesh );
#endif

	return MS::kSuccess;
}
