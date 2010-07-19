#ifndef SCENEXML_HH
#define SCENEXML_HH

/**
 * \file SceneXML.hh
 *
 * \author acoull@wetafx.co.nz
 * \date 06/23/2010
 */

#ifdef WETA
#include "../Physics/ElasticRods/ElasticRod.hh"
#include "../Physics/ElasticRods/AnisotropicRod.hh"
#include "../Physics/ElasticRods/RodUtils.hh"
#include "tinyxml/tinyxml.h"
#else
// Put Columbia includes here
#endif

#include <vector>
#include <tr1/unordered_map>

namespace BASim {

typedef std::tr1::unordered_map< int, Vec3d > FixedVertexData;

struct FrameData
{
    // Data that changes per frame. Currently only fixed vertices
    FixedVertexData fixedVertices;
};

struct InitialRodConfiguration
{
    RodOptions rodOptions;
    Vec3d gravity;
    double massDamping;    
    std::vector< Vec3d > initialRodVertexPositions;    
};

/** Class that handles loading and saving of BASim scene data. It's called scene because it comes
    from Maya. Perhaps WorldXML might be better? */
class SceneXML
{
public:

    /// Default constructor
    SceneXML();

    /// Default Destructor
    ~SceneXML();

    void setInitialSceneState( std::string i_xmlFilepath, 
                               std::vector< InitialRodConfiguration >& i_initialRodConfigurations,
                               std::string i_mayaSceneFilename, double i_stepSize );
    void createHeader();
    void createInitialConfiguration();
    void addFrameData( std::vector< FrameData >& i_frameData );
    void writeFile();

private:
    TiXmlDocument* m_xmlDocument;
    TiXmlElement* m_xmlSceneElement;
    
    std::string m_mayaSceneFilename;
    std::string m_xmlFilePath;
    double m_stepSize;

    std::vector< InitialRodConfiguration > m_initialRodConfigurations;    
  
//    std::vector< FrameData > m_frameData;
};

}

#endif
