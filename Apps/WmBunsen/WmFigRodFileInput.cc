#include "WmFigRodFileInput.hh"

#include <maya/MGlobal.h>
#include <maya/MTime.h>

WmFigRodFileInput::WmFigRodFileInput( MString& i_cacheFilename, WmFigRodGroup& i_rodGroup,
                                      RodOptions& i_rodOptions ) :
    m_cacheFilename( i_cacheFilename ), m_rodGroup( i_rodGroup ), m_rodOptions( i_rodOptions )
{
}

WmFigRodFileInput::~WmFigRodFileInput()
{
}

void WmFigRodFileInput::initialiseRodDataFromInput( MDataBlock& i_dataBlock )
{
    size_t numRodsInFile;
    vector<vector<Vec3d> > rodVertices;
    vector<vector<Vec3d> > unsimulatedRodVertices;

    if ( !readDataFromRodCacheFile( m_cacheFilename, numRodsInFile, rodVertices, 
                                    unsimulatedRodVertices ) )
    {
        return;
    }
 
    m_rodGroup.setIsReadingFromCache( true );

    // FIXME: We need a new addRod that lets us create these cached rods without all the other stuff
    m_rodGroup.addRodsFromCache( rodVertices, m_rodOptions, 10.0 );    
}