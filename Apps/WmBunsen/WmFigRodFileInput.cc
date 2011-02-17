#include "WmFigRodFileInput.hh"

#include <maya/MGlobal.h>
#include <maya/MTime.h>

using namespace std;

WmFigRodFileInput::WmFigRodFileInput( MString& i_cacheFilename, WmFigRodGroup& i_rodGroup,
                                      RodOptions& i_rodOptions ) :
    m_cacheFilename( i_cacheFilename ), m_rodGroup( i_rodGroup ), m_rodOptions( i_rodOptions )
{
    m_simulating = false;
}

WmFigRodFileInput::~WmFigRodFileInput()
{
}

void WmFigRodFileInput::initialiseRodDataFromInput( MDataBlock& i_dataBlock )
{
    int numRodsInFile;
    vector<vector<BASim::Vec3d> > rodVertices;
    vector<vector<BASim::Vec3d> > unsimulatedRodVertices;

    if ( !readDataFromRodCacheFile( m_cacheFilename, numRodsInFile, rodVertices, 
                                    unsimulatedRodVertices ) )
    {
        cerr << "Failed to initialise rod data from cache file\n";
        return;
    }
 
    cerr << "initialising rods, " << numRodsInFile << " rods in file\n";
    m_rodGroup.setIsReadingFromCache( true );

    // FIXME: We need a new addRod that lets us create these cached rods without all the other stuff
    m_rodGroup.addRodsFromCache( rodVertices, m_rodOptions, 10.0 );    
}
