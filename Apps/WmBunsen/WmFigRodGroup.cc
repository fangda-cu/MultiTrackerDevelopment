#include "WmFigRodGroup.hh"

WmFigRodGroup::WmFigRodGroup() : m_isReadingFromCache( false ), m_simulationNeedsReset( false )
{
}

WmFigRodGroup::~WmFigRodGroup()
{
    removeAllRods();
}

size_t WmFigRodGroup::addRod()
{
    // This is a place holder rod really has no data but exists so that the input curves/strand
    // indices line up with which rod is which.
    RodData* rodData = new RodData();
    m_rodData.push_back( rodData );

    return m_rodData.size() - 1;
}

size_t WmFigRodGroup::addRod( std::vector<Vec3d>& i_rodVertices, RodOptions& i_rodOptions, 
                              double i_massDamping, bool i_isFromCache )
{
    // We have rod vertices so lets build a real rod.
    
    Vec3d gravity( 0.0, -980.0, 0.0 );

    RodData* rodData = new RodData( i_rodOptions, i_rodVertices, i_massDamping, gravity, i_isFromCache );

    m_rodData.push_back( rodData );

    return m_rodData.size() - 1;
}