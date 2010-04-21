#include "WmFigRodFileInput.hh"

#include <maya/MGlobal.h>
#include <maya/MTime.h>

WmFigRodFileInput::WmFigRodFileInput( MObject& i_filePathAttribute, MObject& i_timeAttribute ) :
    m_filePathAttribute( i_filePathAttribute ), m_timeAttribute( i_timeAttribute )
{
}

WmFigRodFileInput::~WmFigRodFileInput()
{
}

void WmFigRodFileInput::initialiseRodDataFromInput( MDataBlock& i_dataBlock, 
    std::vector<RodData*>* i_pRodData  )
{
    MString cacheFilename = getCacheFilename( i_dataBlock );

    if ( i_pRodData == NULL )
        return;     // No rod data, nowhere to store rod data

    size_t numRodsInFile;
    vector<vector<Vec3d> > rodVertices;
    vector<vector<Vec3d> > unsimulatedRodVertices;

    if ( !readDataFromRodCacheFile( cacheFilename, numRodsInFile, rodVertices, 
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

MString WmFigRodFileInput::getCacheFilename( MDataBlock& i_dataBlock )
{
    MStatus stat;

    // FIXME: This doesn't deal nicely with fractional frames in the way Barbershop does.

    MString filePath = i_dataBlock.inputValue( m_filePathAttribute, &stat ).asString();
    CHECK_MSTATUS( stat );

    if ( filePath == "" )
        return "";

    double currentTime = i_dataBlock.inputValue( m_timeAttribute, &stat ).asTime().value();
    CHECK_MSTATUS( stat );

    MString fileName = ( filePath + "." + currentTime + ".fig" );

    return fileName;
}

FILE* WmFigRodFileInput::readNumberOfRodsFromFile( const MString i_cacheFilename, 
    size_t& o_numRodsInFile, bool closeFileAfterReading )
{
    FILE *fp = NULL;

    if ( i_cacheFilename == "" )
    {
        MGlobal::displayError( "Empty cache path, can't read in rod data!" );
        return NULL;
    }

    fp = fopen( i_cacheFilename.asChar(), "r" );
    if ( fp == NULL )
    {
        MGlobal::displayError( MString( "Problem opening file " + i_cacheFilename + " for reading." ) );
        return NULL;
    }

/*    int fileFormat;
    fread( &fileFormat, sizeof( int ), 1, fp );
    if ( fileFormat != FILE_FORMAT_VERSION )
    {
        Mglobal::displayError( "Unsupported cache file format version\n" );
        fclose( fp );
        fp = NULL
        o_numRodsInFile = 0;

    }*/

    size_t ret = fread( &o_numRodsInFile, sizeof( size_t ), 1, fp );

    if ( closeFileAfterReading )
    {
        fclose( fp );
        fp = NULL;
    }

    return fp;
}

bool WmFigRodFileInput::readDataFromRodCacheFile( const MString i_cacheFilename, size_t& o_numRodsInFile,
    vector<vector<Vec3d> >& o_rodVertices, vector<vector<Vec3d> >& o_unsimulatedRodVertices )
{

    FILE* fp = readNumberOfRodsFromFile( i_cacheFilename, o_numRodsInFile, false );
    if ( fp == NULL )
        return false;

    o_rodVertices.resize( o_numRodsInFile );
    o_unsimulatedRodVertices.resize( o_numRodsInFile );

    for ( size_t r=0; r<o_numRodsInFile; r++ )
    {
        size_t numVertices;
        size_t ret = fread( &numVertices, sizeof( size_t ), 1, fp );

        o_rodVertices[ r ].resize( numVertices );
        o_unsimulatedRodVertices[ r ].resize( numVertices );

        for ( size_t v=0; v<numVertices; v++ )
        {
            double pos[3];

            ret = fread( &pos[0], sizeof( double ), 3, fp );
            o_rodVertices[ r ][ v ] = Vec3d( pos[0], pos[1], pos[2] );

            ret = fread( &pos[0], sizeof( double ), 3, fp );
            o_unsimulatedRodVertices[ r ][ v ] = Vec3d( pos[0], pos[1], pos[2] );
        }
    }

    fclose( fp );
    return true;
}
