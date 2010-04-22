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
     WmFigRodFileInput( MString& i_cacheFilename );
   
    /**
      * @brief Default destructor
      */
     ~WmFigRodFileInput ();

    /**
      * @brief Initialise rod data from input types.
      */
    virtual void initialiseRodDataFromInput( MDataBlock& i_dataBlock, std::vector<RodData*>* i_pRodData  );    
    
private:
    MString m_cacheFilename;
};

#endif // WMFIGRODFILEINPUT_H_
