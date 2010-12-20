#ifndef WMFIGRODFILEINPUT_H_
#define WMFIGRODFILEINPUT_H_

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
#include "WmFigRodFileIO.hh"
#include "WmFigRodGroup.hh"

/**
  * @brief Class to load rod data from a cached file and make it look like live data for the rod node
  *
  *  
  *
  *
  */ 

class WmFigRodFileInput : public WmFigRodInputType, WmFigRodFileIO
{ 
public:
     /**
      * @brief Default constructor
      */
     WmFigRodFileInput( MString& i_cacheFilename, WmFigRodGroup& i_rodGroup, 
                        RodOptions& i_rodOptions );
   
    /**
      * @brief Default destructor
      */
     ~WmFigRodFileInput ();

    /**
      * @brief Initialise rod data from input types.
      */
    virtual void initialiseRodDataFromInput( MDataBlock& i_dataBlock );    

private:
    MString m_cacheFilename;
    WmFigRodGroup& m_rodGroup;
    RodOptions m_rodOptions;
};

#endif // WMFIGRODFILEINPUT_H_
