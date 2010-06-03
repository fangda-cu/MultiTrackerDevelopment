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

#include "WmFigRodInputType.hh"

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
     WmFigRodNurbsInput( MObject& i_nurbsAttribute, bool i_lockFirstEdgeToInput );

    /**
      * @brief Default destructor
      */
     ~WmFigRodNurbsInput();

    /**
      * @brief Initialise rod data from input types.
      */
    virtual void initialiseRodDataFromInput( MDataBlock& i_dataBlock );
    virtual void updateRodDataFromInput( MDataBlock& i_dataBlock );
    virtual size_t numberOfInputs( MDataBlock& i_dataBlock );

private:

    MObject& m_inputNurbsAttribute;
    bool m_lockFirstEdgeToInput;
};

#endif // WMFIGRODNURBSINPUT_H_

 
