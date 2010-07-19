#include "SceneXML.hh"

namespace BASim {

SceneXML::SceneXML() : m_xmlDocument( NULL ), m_xmlSceneElement( NULL ), m_mayaSceneFilename( "" ),
    m_xmlFilePath( "" )
{
    m_initialRodConfigurations.clear();
    //m_fixedRodVertices.clear();
}

SceneXML::~SceneXML()
{
    if ( m_xmlDocument != NULL )
    {   
        delete m_xmlDocument;
    }
}

void SceneXML::setInitialSceneState( std::string i_xmlFilePath, std::vector< InitialRodConfiguration >& i_initialRodConfigurations,
    std::string i_mayaSceneFilename, double i_stepSize )
{
    if ( m_xmlDocument != NULL )
    {   
        delete m_xmlDocument;
    }
    m_xmlDocument = new TiXmlDocument;

    m_xmlFilePath = i_xmlFilePath;
    m_initialRodConfigurations = i_initialRodConfigurations;
    m_mayaSceneFilename = i_mayaSceneFilename;
    m_stepSize = i_stepSize;

    createHeader();
    createInitialConfiguration();
}

void SceneXML::addFrameData( std::vector< FrameData >& i_frameData )
{
    TiXmlElement* frameElement = new TiXmlElement( "Frame" );
    m_xmlSceneElement->LinkEndChild( frameElement );

    for ( size_t r=0; r<i_frameData.size(); ++r )
    {
        TiXmlElement* rodElement = new TiXmlElement( "Rod" );
        frameElement->LinkEndChild( rodElement );
       
        for ( FixedVertexData::iterator it=i_frameData[ r ].fixedVertices.begin();
              it != i_frameData[ r ].fixedVertices.end(); ++it )
        {
            int vertexIndex = it->first;
            Vec3d vertexPosition = it->second;

            TiXmlElement* vertexElement = new TiXmlElement( "Vertex" );
            rodElement->LinkEndChild( vertexElement );

            vertexElement->SetAttribute( "index", vertexIndex );
            vertexElement->SetDoubleAttribute( "x", vertexPosition[ 0 ] );
            vertexElement->SetDoubleAttribute( "y", vertexPosition[ 1 ] );
            vertexElement->SetDoubleAttribute( "z", vertexPosition[ 2 ] );
        }
    }
}

void SceneXML::createHeader()
{
    TiXmlDeclaration* decl = new TiXmlDeclaration( "1.0", "", "" );

    TiXmlElement* fileFormatVersionElement = new TiXmlElement( "File Format Version" );  
    fileFormatVersionElement->LinkEndChild( new TiXmlText( "1.0" ) );  
    m_xmlDocument->LinkEndChild( fileFormatVersionElement );

    m_xmlSceneElement = new TiXmlElement( "Scene" );
    m_xmlDocument->LinkEndChild( m_xmlSceneElement );
}

void SceneXML::createInitialConfiguration()
{
    TiXmlElement* mayaSceneElement = new TiXmlElement( "Maya File" );
    mayaSceneElement->LinkEndChild( new TiXmlText( m_mayaSceneFilename.c_str() ) );
    m_xmlSceneElement->LinkEndChild( mayaSceneElement );  

    TiXmlElement* initializationElement = new TiXmlElement( "Initialization" );
    m_xmlSceneElement->LinkEndChild( initializationElement );

    //TiXmlElement* optionsElement = new TiXmlElement( "Options" );
    //initializationElement->LinkEndChild( optionsElement );
    
    for ( size_t r=0; r<m_initialRodConfigurations.size(); ++r )
    {
        TiXmlElement* rodElement = new TiXmlElement( "Rod" );
        initializationElement->LinkEndChild( rodElement );
        
        TiXmlElement* rodOptionsElement = new TiXmlElement( "Options" );
        rodElement->LinkEndChild( rodOptionsElement );
        
        RodOptions& rodOptions = m_initialRodConfigurations[ r ].rodOptions;

        rodOptionsElement->SetAttribute( "numVertices", rodOptions.numVertices );
        rodOptionsElement->SetDoubleAttribute( "Density", rodOptions.density );
        rodOptionsElement->SetDoubleAttribute( "radiusA", rodOptions.radiusA );
        rodOptionsElement->SetDoubleAttribute( "radiusB", rodOptions.radiusB );
        rodOptionsElement->SetDoubleAttribute( "radiusScale", rodOptions.radiusScale );    
        rodOptionsElement->SetDoubleAttribute( "YoungsModulus", rodOptions.YoungsModulus );
        rodOptionsElement->SetDoubleAttribute( "ShearModulus", rodOptions.ShearModulus );
        rodOptionsElement->SetDoubleAttribute( "Viscosity", rodOptions.viscosity );
        rodOptionsElement->SetAttribute( "anisotropic", rodOptions.anisotropic );
        rodOptionsElement->SetAttribute( "inextensible", rodOptions.inextensible );
        rodOptionsElement->SetAttribute( "refFrame", rodOptions.refFrame );
        
        rodOptionsElement->SetDoubleAttribute( "MassDamping", m_initialRodConfigurations[ r ].massDamping );

        Vec3d gravity = m_initialRodConfigurations[ r ].gravity;

        TiXmlElement* gravityElement = new TiXmlElement( "Gravity" );
        rodOptionsElement->LinkEndChild( gravityElement );
        gravityElement->SetDoubleAttribute( "x", gravity[ 0 ] );
        gravityElement->SetDoubleAttribute( "y", gravity[ 1 ] );
        gravityElement->SetDoubleAttribute( "z", gravity[ 2 ] );
    
        for ( size_t v=0; v<m_initialRodConfigurations[ r ].initialRodVertexPositions.size(); ++v )
        {
            TiXmlElement* vertexElement = new TiXmlElement( "Vertex" );
            rodElement->LinkEndChild( vertexElement );

            Vec3d vertex = m_initialRodConfigurations[ r ].initialRodVertexPositions[ v ];

            vertexElement->SetAttribute( "index", (int)v );
            vertexElement->SetDoubleAttribute( "x", vertex[ 0 ] );
            vertexElement->SetDoubleAttribute( "y", vertex[ 1 ] );
            vertexElement->SetDoubleAttribute( "z", vertex[ 2 ] );
        }
    }
}

void SceneXML::writeFile()
{
    m_xmlDocument->SaveFile( m_xmlFilePath.c_str() );    
}


} // namespace BASim

