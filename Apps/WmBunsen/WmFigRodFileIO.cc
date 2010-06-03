#include "WmFigRodFileIO.hh"

#include <maya/MGlobal.h>

WmFigRodFileIO::WmFigRodFileIO()
{
}

WmFigRodFileIO::~WmFigRodFileIO()
{
}

/* static */ FILE* WmFigRodFileIO::readNumberOfRodsFromFile( const MString i_cacheFilename, 
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

    cerr << "Reading file " << i_cacheFilename << endl;

    size_t ret = fread( &o_numRodsInFile, sizeof( size_t ), 1, fp );

    if ( closeFileAfterReading )
    {
        fclose( fp );
        fp = NULL;
    }

    return fp;
}

/* static */ bool WmFigRodFileIO::readDataFromRodCacheFile( const MString i_cacheFilename, 
    size_t& o_numRodsInFile, vector<vector<Vec3d> >& o_rodVertices, 
    vector<vector<Vec3d> >& o_unsimulatedRodVertices )
{
    FILE* fp = readNumberOfRodsFromFile( i_cacheFilename, o_numRodsInFile, false );
    if ( fp == NULL )
    {
        return false;
    }

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


/* static */ void WmFigRodFileIO::updateRodDataFromCacheFile( MString i_cacheFileName, WmFigRodGroup& i_rodGroup )
{
    size_t numRodsInFile;
    vector<vector<Vec3d> > rodVertices;
    vector<vector<Vec3d> > unsimulatedRodVertices;

    if ( !readDataFromRodCacheFile( i_cacheFileName, numRodsInFile, rodVertices,
                                    unsimulatedRodVertices ) )
    {
        cerr << "Failed to read\n";
        return;
    }

    cerr << "numRodsInFile << " << numRodsInFile << endl;
    cerr << "i_rodGroup.numberOfRods()  << " << i_rodGroup.numberOfRods() << endl;

    if ( numRodsInFile != i_rodGroup.numberOfRods() )
    {
        MGlobal::displayError( "Rewind simulation to reset before cache file can be read" );
        return;
    }

    for ( size_t r=0; r<numRodsInFile; r++ )
    {
        
        BASim::ElasticRod* rod = i_rodGroup.elasticRod( r );
        if ( rod == NULL )
        {
            MGlobal::displayError( "Rewind simulation to reset before cache file can be read" );
            continue;
        }

        size_t numVertices = rodVertices[ r ].size();

        if ( (int)numVertices != rod->nv() )
        {
            MGlobal::displayError( "Can't read in cached rod data as number of vertices do not match rods in scene\n" );
            continue;
        }

        for ( size_t v=0; v<numVertices; v++ )
        {
            rod->setVertex( (int)v, rodVertices[ r ][ v ] );
        }
    }
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

    // FIXME: create a proper cache file format. This one is bad, it at least needs a proper header.
    
    size_t numberOfRealRods = i_rodGroup.numberOfRealRods();
    size_t totalNumberOfRods = i_rodGroup.numberOfRods();
    fwrite( &numberOfRealRods, sizeof( size_t ), 1, fp );

    for ( size_t r=0; r<totalNumberOfRods; r++ )
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
            
            for ( size_t v=0; v<numVertices; v++ )
            {
            double pos[3];

            // Wonder if its safe to write Vec3ds. Need to check what's in them.
            // Really should package all this and write it as one.
            Vec3d vertex = rod->getVertex( (unsigned int)v );
            pos[ 0 ] = vertex.x();
            pos[ 1 ] = vertex.y();
            pos[ 2 ] = vertex.z();

            fwrite( &pos[0], sizeof( double ), 3, fp );

            // Now store the unsimulated positions so that we can send them to Barbershop
            // when we want to simulate with it.
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
