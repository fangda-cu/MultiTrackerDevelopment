#include "WmFigRodFileIO.hh"

#include <maya/MGlobal.h>


const int WmFigRodFileIO::m_versionNumber = 1;

/*
FILE FORMAT NOTES
=================

For now it's just a four character tag ("WFRC"), followed
by a binary integer giving the file format version, followed
by a binary size_t field giving the number of rods in the file.

Then for each rod in the file, number of vertices as a size_t,
then the vertices are written as six floats: 3 for the
xyz of the simulated position, followed by 3 for the xyz
of the unsimulated position.

TAG
VERSION NUMBER
NUM RODS

NUM VERTICES
6 COORDS
6 COORDS
6 COORDS
...
NUM VERTICES
6 COORDS
6 COORDS
6 COORDS
...

etc.

*/


WmFigRodFileIO::WmFigRodFileIO()
{
}

WmFigRodFileIO::~WmFigRodFileIO()
{
}

/* static */ FILE* WmFigRodFileIO::readNumberOfRodsFromFile( const MString i_cacheFilename, 
    int& o_numRodsInFile, bool closeFileAfterReading )
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

    // skip tag and version number
    fseek(fp, 4 + sizeof(int), SEEK_SET);

    size_t numRodsInFile;
    size_t ret = fread( &numRodsInFile, sizeof( size_t ), 1, fp );
    o_numRodsInFile = (int)numRodsInFile;
    assert ( numRodsInFile == o_numRodsInFile );

    if ( closeFileAfterReading )
    {
        fclose( fp );
        fp = NULL;
    }

    return fp;
}

/* static */ bool WmFigRodFileIO::checkTagAndVersion( const MString i_cacheFilename, int &versionNumber )
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

     char tagBuf[5];
     size_t ret = fread( tagBuf, 1, 4, fp );
     tagBuf[4] = '\0';

     if(0 != strcmp(tagBuf, ROD_CACHE_FILE_TAG))
     {
         fclose( fp );
         MGlobal::displayError( MString( i_cacheFilename + " doesn't appear to be a figaro cache file (tag missing)." ) );
         return false;
     }


     ret = fread( &versionNumber, sizeof( int ), 1, fp );

     if(ret != 1)
     {
         fclose( fp );
         MGlobal::displayError( MString( "Error reading version number from figaro cache file " + i_cacheFilename + ".") );
         return false;
     }

     return true;
 }


/* static */ bool WmFigRodFileIO::readDataFromRodCacheFile( const MString i_cacheFilename,
    int& o_numRodsInFile, vector<vector<BASim::Vec3d> >& o_rodVertices, 
    vector<vector<BASim::Vec3d> >& o_unsimulatedRodVertices )
{
    int versionNumber = -1;
    bool isValid = checkTagAndVersion( i_cacheFilename, versionNumber );

    if(!isValid)
    {
        // this may an un-versioned rod cache file; let's warn, but still try to read it
        MGlobal::displayError( i_cacheFilename + " doesn't appear to be a valid rod cache file." );

        return false;
    }

    if(versionNumber != m_versionNumber)
    {
        char versionFromFile[10];
        sprintf(versionFromFile, "%d", versionNumber);
        char versionWeNeed[10];
        sprintf(versionWeNeed, "%d", m_versionNumber);
        MGlobal::displayError( MString( "Error: " + i_cacheFilename + " is file version " + MString(versionFromFile) + "; we need version " + MString(versionWeNeed) ) );

        return false;
    }

    FILE* fp = readNumberOfRodsFromFile( i_cacheFilename, o_numRodsInFile, false );
    if ( fp == NULL )
    {
        return false;
    }

    o_rodVertices.resize( o_numRodsInFile );
    o_unsimulatedRodVertices.resize( o_numRodsInFile );

    for ( int r=0; r<o_numRodsInFile; r++ )
    {
        size_t numVertices;
        size_t ret = fread( &numVertices, sizeof( size_t ), 1, fp );

        o_rodVertices[ r ].resize( numVertices );
        o_unsimulatedRodVertices[ r ].resize( numVertices );

        for ( size_t v=0; v<numVertices; v++ )
        {
            double pos[3];

            ret = fread( &pos[0], sizeof( double ), 3, fp );
            o_rodVertices[ r ][ v ] = BASim::Vec3d( pos[0], pos[1], pos[2] );

            ret = fread( &pos[0], sizeof( double ), 3, fp );
            o_unsimulatedRodVertices[ r ][ v ] = BASim::Vec3d( pos[0], pos[1], pos[2] );
        }
    }

    fclose( fp );
    return true;
}


/* static */ bool WmFigRodFileIO::updateRodDataFromCacheFile( MString i_cacheFileName, WmFigRodGroup& i_rodGroup )
{
    int numRodsInFile;
    vector<vector<BASim::Vec3d> > rodVertices;
    vector<vector<BASim::Vec3d> > unsimulatedRodVertices;

    if ( !readDataFromRodCacheFile( i_cacheFileName, numRodsInFile, rodVertices,
                                    unsimulatedRodVertices ) )
    {
        cerr << "Failed to read cache file " << i_cacheFileName << ".\n";
        return false;
    }

    //cerr << "numRodsInFile << " << numRodsInFile << endl;
    //cerr << "i_rodGroup.numberOfRods()  << " << i_rodGroup.numberOfRods() << endl;

    if ( numRodsInFile != i_rodGroup.numberOfRods() )
    {
        MGlobal::displayError( "Rewind simulation to reset before cache file can be read" );
        return false;
    }

    for ( int r=0; r<numRodsInFile; r++ )
    {        
        BASim::ElasticRod* rod = i_rodGroup.elasticRod( r );
        if ( rod == NULL )
        {
            MGlobal::displayError( "Rewind simulation to reset before cache file can be read" );
            continue;
        }

        int numVertices = (int)rodVertices[ r ].size();

        if ( (int)numVertices != rod->nv() )
        {
            MGlobal::displayError( "Can't read in cached rod data as number of vertices do not match rods in scene\n" );
            continue;
        }

        for ( int v=0; v<numVertices; v++ )
        {
            rod->setVertex( v, rodVertices[ r ][ v ] );
        }
    }
    
    return true;
}

/* static */ void WmFigRodFileIO::writeRodDataToCacheFile( MString& i_cacheFileame, 
     WmFigRodGroup& i_rodGroup  )
{
    if ( i_cacheFileame == "" )
    {
        MGlobal::displayError( "Empty cache path, not writing anything to disk." );
        return;
    }

    cerr << "Writing sim data to disk in file: '" << i_cacheFileame << "'\n";

    FILE *fp;
    fp = fopen( i_cacheFileame.asChar(), "w" );
    if ( fp == NULL )
    {
        MGlobal::displayError( MString( "Problem opening file " + i_cacheFileame + " for writing." ) );
        return;
    }


    fwrite( ROD_CACHE_FILE_TAG, 1, 4, fp );
    fwrite( &m_versionNumber, sizeof(int), 1, fp );
    
    size_t numberOfRealRods = i_rodGroup.numberOfRealRods();
    int totalNumberOfRods = i_rodGroup.numberOfRods();
    fwrite( &numberOfRealRods, sizeof( size_t ), 1, fp );

    for ( int r=0; r<totalNumberOfRods; r++ )
    {
        if ( !i_rodGroup.isPlaceHolderRod( r ) )
        {
            // Write number of vertices in this rod
            BASim::ElasticRod* rod = i_rodGroup.elasticRod( r );
    
            if ( rod == NULL )
            {
                cerr << "WTF, rod is NULL \n";
                continue;
            }

            size_t numVertices = rod->nv();
            fwrite( &numVertices, sizeof( size_t ), 1, fp );
            
            for ( int v=0; v<numVertices; v++ )
            {
                double pos[3];

                // Wonder if its safe to write BASim::Vec3ds. Need to check what's in them.
                // Really should package all this and write it as one.
                BASim::Vec3d vertex = rod->getVertex( (unsigned int)v );
                pos[ 0 ] = vertex.x();
                pos[ 1 ] = vertex.y();
                pos[ 2 ] = vertex.z();

                fwrite( &pos[0], sizeof( double ), 3, fp );

                // Now store the unsimulated positions so that we can send them to Barbershop
                // when we want to simulate with it.
                // Boundary Condition
                vertex = i_rodGroup.nextVertexPosition( r, v );
                pos[ 0 ] = vertex.x();
                pos[ 1 ] = vertex.y();
                pos[ 2 ] = vertex.z();

                fwrite( &pos[0], sizeof( double ), 3, fp );

                // We should write velocities and anything else we need into the file too so that
                // we can restart a simulation from a cached file.
            }
        }
    }
    fclose ( fp );
}
