#ifndef WMFIGRODNURBSINPUT_H_
#define WMFIGRODNURBSINPUT_H_

/**
  * @file Apps/WmFigaro/WmFigRodNurbsInput.h
  * @author Alasdair Coull (acoull@wetafx.co.nz)
  * @date 20-04-2010
  */

#include <maya/MDataBlock.h>
#include <maya/MDataHandle.h>
#include <maya/MFnVectorArrayData.h>
#include <maya/MVectorArray.h>
#include <maya/MFnNurbsCurve.h>
#include <maya/MPoint.h>

#include <set>

#include "WmFigRodInputType.hh"
#include "WmFigRodGroup.hh"

/**
  * @brief Class to convert Maya NURBS curves into a format suitable for creating rods
  *
  *  
  *
  *
  */ 

class WmFigRodNurbsInput : public WmFigRodInputType
{ 
public:
     /**
      * @brief Default constructor
      */
     WmFigRodNurbsInput( MObject& i_nurbsAttribute, bool i_lockFirstEdgeToInput, WmFigRodGroup& i_rodGroup,
                         double i_vertexSpacing, double i_minimumRodLength, RodOptions& i_rodOptions,
                         double i_massDamping, BASim::Vec3d& i_gravity,  
                         RodTimeStepper::Method i_solverType, std::set< int >& i_simulationSet,
                         const bool i_doReverseHairdo );

    /**
      * @brief Default destructor
      */
     ~WmFigRodNurbsInput();

    /**
      * @brief Initialise rod data from input types.
      */
    virtual void initialiseRodDataFromInput( MDataBlock& i_dataBlock );
    virtual void updateRodDataFromInput( MDataBlock& i_dataBlock );
    virtual int numberOfInputs( MDataBlock& i_dataBlock );

    void getAndResampleInputCurves(  MDataBlock& i_dataBlock, vector< vector<BASim::Vec3d > >& o_inputCurveVertices );
private:

    MObject& m_inputNurbsAttribute;
    bool m_lockFirstEdgeToInput;
    WmFigRodGroup& m_rodGroup;
    double m_vertexSpacing;
    double m_minimumRodLength;
    RodOptions m_rodOptions;
    double m_massDamping;
    BASim::Vec3d m_gravity;
    RodTimeStepper::Method m_solverType;
    std::set< int >& m_simulationSet;
    bool m_doReverseHairdo;
};

#endif // WMFIGRODNURBSINPUT_H_

 
