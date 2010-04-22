#include "WmFigRodFileInput.hh"

#include <maya/MGlobal.h>
#include <maya/MTime.h>

WmFigRodFileInput::WmFigRodFileInput( MString& i_cacheFilename ) :
    m_cacheFilename( i_cacheFilename )
{
}

WmFigRodFileInput::~WmFigRodFileInput()
{
}

void WmFigRodFileInput::initialiseRodDataFromInput( MDataBlock& i_dataBlock, 
    std::vector<RodData*>* i_pRodData  )
{
    if ( i_pRodData == NULL )
        return;     // No rod data, nowhere to store rod data

    size_t numRodsInFile;
    vector<vector<Vec3d> > rodVertices;
    vector<vector<Vec3d> > unsimulatedRodVertices;

    if ( !readDataFromRodCacheFile( m_cacheFilename, numRodsInFile, rodVertices, 
                                    unsimulatedRodVertices ) )
        return;

    if ( numRodsInFile != i_pRodData->size() )
    {
        MGlobal::displayError( "Wrong number of rods in file, have not implemented changing it yet" );
        MGlobal::displayError( MString( "Expected " ) +  ( double) i_pRodData->size() + ", got " 
                               + ( double )numRodsInFile );
        return;
    }

    for ( size_t r=0; r<numRodsInFile; r++ )
    {
        size_t numVertices = rodVertices[ r ].size();
        
        // Resize all the data vectors to make sure they're large enough for the number of cvs
        // we're getting from the file
        (*i_pRodData)[ r ]->allocateStorage( numVertices );

        (*i_pRodData)[ r ]->resetVertexPositions( rodVertices[ r ] );
    }
}