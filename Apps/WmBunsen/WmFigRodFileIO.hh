#ifndef WMFIGRODFILEIO_H_
#define WMFIGRODFILEIO_H_

/**
  * @file Apps/WmFigaro/WmFigRodFileIO.h
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

// (W)eta (F)igaro (R)od (C)ache
#define ROD_CACHE_FILE_TAG "WFRC"

/**
  * @brief Base class providing file IO operations for Figaro rod cache files
  *
  *  
  *
  *
  */ 

class WmFigRodFileIO
{ 
public:
     /**
      * @brief Default constructor
      */
     WmFigRodFileIO();
   
    /**
      * @brief Default destructor
      */
     ~WmFigRodFileIO();

    static  bool checkTagAndVersion( const MString i_cacheFilename, int &versionNumber);

    static FILE* readNumberOfRodsFromFile( const MString i_cacheFilename, int& o_numRodsInFile, 
            bool closeFileAfterReading = true );
    
    static bool readDataFromRodCacheFile( const MString i_cacheFilename, int& o_numRodsInFile,
          vector<vector<BASim::Vec3d> >& o_rodVertices, vector<vector<BASim::Vec3d> >& o_unsimulatedRodVertices );
    
    static void updateRodDataFromCacheFile( MString i_cacheFileName, WmFigRodGroup& i_rodGroup );

    static void writeRodDataToCacheFile( MString& i_cacheFileame, WmFigRodGroup& i_rodGroup );

private:

    static const int m_versionNumber;
};

#endif // WMFIGRODFILEIO_H_

    
