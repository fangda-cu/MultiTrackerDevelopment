/**
 * \file TopologicalObjectSerializer.cc
 *
 * \author smith@cs.columbia.edu
 * \date 07/16/2010
 */

#include "TopologicalObjectSerializer.hh"

namespace BASim 
{

void TopologicalObjectSerializer::saveTopologicalObject( const TopologicalObject& topobj, const std::string& filename )
{
  // Attempt to open the file for binary output
  std::ofstream of(filename.c_str(),std::ios::binary);
  if( !of.is_open() ) 
  {
    std::cerr << "\033[31;1mWARNING IN TOPOLOGICALOBJECTSERIALIZER:\033[m Failed to open " << filename << " for output. No data will be saved." << std::endl;
    return;
  }

  appendTopologicalObjectToFile( topobj, of );
  
  of.close();

  //std::cout << "\033[35;1mTOPOLOGICALOBJECTSERIALIZED MESSAGE:\033[m Saved topological object to: " << filename << std::endl;
}

void TopologicalObjectSerializer::appendTopologicalObjectToFile( const TopologicalObject& topobj, std::ofstream& of )
{
  assert( of.is_open() );
  
  // Write a header containing the object id (identifies topobj as a rod, cloth, etc)
  serializeHeader( of, topobj );
  // Write the number of object, vertex, face, and edge properties
  serializeNumProperties( of, topobj );
  // Write out the new graph structure (which was stored as properties in the old structure)
  TopologicalObject::serializeStructure(of, topobj);
  // Write the size of vertex, edge, and face properties (object properties are always size 1)
  //serializePropertySizes( of, topobj );
  
  // Write the vertex properties
  const PropertyContainer::Properties& vprops = getVertexPropertyContainer(topobj).getProperties();
  for( int i = 0; i < (int) vprops.size(); ++i ) BASim::serializeProperty( of, vprops[i] );
  // Write the edge properties
  const PropertyContainer::Properties& eprops = getEdgePropertyContainer(topobj).getProperties();
  for( int i = 0; i < (int) eprops.size(); ++i ) BASim::serializeProperty( of, eprops[i] );
  // Write the face properties
  const PropertyContainer::Properties& fprops = getFacePropertyContainer(topobj).getProperties();
  for( int i = 0; i < (int) fprops.size(); ++i ) BASim::serializeProperty( of, fprops[i] );
  // Write the object properties
  const PropertyContainer::Properties& oprops = getObjectPropertyContainer(topobj).getProperties();
  for( int i = 0; i < (int) oprops.size(); ++i ) BASim::serializeProperty( of, oprops[i] );
}
  
void TopologicalObjectSerializer::loadTopologicalObject( ElasticRod** tpobj, const std::string& filename )
{
  assert( *tpobj == NULL );

  // Attempt to open the file for binary output
  std::ifstream ifs(filename.c_str(),std::ios::binary);
  if( !ifs.is_open() ) 
  {
    std::cerr << "\033[31;1mWARNING IN TOPOLOGICALOBJECTSERIALIZER:\033[m Failed to open " << filename << " for reading. No data will be loaded." << std::endl;
    return;
  }

  loadTopologicalObjectFromFile( tpobj, ifs );

  ifs.close();
}

void TopologicalObjectSerializer::loadTopologicalObjectFromFile( ElasticRod** tpobj, std::ifstream& ifs )
{
  TopologicalObject* tempptr = *tpobj;
  loadGenericTopologicalObject( &tempptr, ifs );
  assert( tempptr != NULL );
  *tpobj = dynamic_cast<ElasticRod*>(const_cast<TopologicalObject*>(tempptr));
  assert( *tpobj != NULL );
  
  // 'Reattach' properties within rods
  ElasticRod::RodForces& forces = (*tpobj)->getForces();
  ElasticRod::RodForces::iterator fIt;
  for( fIt = forces.begin(); fIt != forces.end(); ++fIt ) (*fIt)->reattatchProperties();
}

void TopologicalObjectSerializer::loadGenericTopologicalObject( TopologicalObject** tpobj, std::ifstream& ifs )
{
  assert( (*tpobj) == NULL );
  assert( ifs.is_open() );
  
  //std::cout << "Inside loadTopologicalObjectFromFile" << std::endl;
  
  // Load the id that identifies the object as a rod, cloth, etc
  ObjectID objtypid = UNRECOGNIZEDOBJECTID;
  loadHeader( ifs, tpobj, objtypid );
  assert( (*tpobj) != NULL );
  assert( objtypid == RODID || objtypid == TRIMESHID );
  
  //std::cout << "Created a new rod" << std::endl;
  
  // Load the number of object, vertex, edge, and face properties
  int noprps = -1; int nvprps = -1; int neprps = -1; int nfprps = -1;
  loadNumProperties( ifs, noprps, nvprps, neprps, nfprps );
  assert( noprps >= 0 ); assert( nvprps >= 0 ); assert( neprps >= 0 ); assert( nfprps >= 0 );
  
  //std::cout << "<oprops,vprops,eprops,fprops>: " << noprps << " " << nvprps << " " << neprps << " " << nfprps << std::endl;
  
  //load the graph structure data
  TopologicalObject::loadStructure(ifs, **tpobj);
  
  // Load the number of vertices, edges, and faces (which correspond to the property sizes) 
  /*int nv = -1; int ne = -1; int nf = -1;
  loadPropertySizes( ifs, nv, ne, nf );
  assert( nv >= 0 ); assert( ne >= 0 ); assert( nf >= 0 );*/
  
  //std::cout << "<nv,ne,nf>: " << nv << " " << ne << " " << nf << std::endl;
  
  // Note: PropertyContainer::size() gives the size of the container, NOT the size of the properties.
  PropertyContainer& vctr = getVertexPropertyContainer(**tpobj);
  vctr.clear();
  PropertyContainer& ectr = getEdgePropertyContainer(**tpobj);
  ectr.clear();
  PropertyContainer& fctr = getFacePropertyContainer(**tpobj);
  fctr.clear();
  PropertyContainer& octr = getObjectPropertyContainer(**tpobj);
  octr.clear();
  
  // Load vertex properties
  assert( vctr.size() == 0 );    
  for( int i = 0; i < (int) nvprps; ++i )
  {
    PropertyBase* prop = NULL;
    loadProperty( ifs, &prop, *tpobj, (**tpobj).nv() );
    vctr.getProperties().push_back(prop);
  }
  assert( (int) vctr.size() == nvprps );
  
  //std::cout << "Loaded vertex properties." << std::endl;
  
  // Load edge properties
  assert( ectr.size() == 0 );    
  for( int i = 0; i < (int) neprps; ++i )
  {
    PropertyBase* prop = NULL;
    loadProperty( ifs, &prop, *tpobj, (**tpobj).ne() );
    ectr.getProperties().push_back(prop);
  }
  assert( (int) ectr.size() == neprps );
  
  //std::cout << "Loaded edge properties." << std::endl;
  
  // Load face properties
  assert( fctr.size() == 0 );    
  for( int i = 0; i < (int) nfprps; ++i )
  {
    PropertyBase* prop = NULL;
    loadProperty( ifs, &prop, *tpobj, (**tpobj).nf() );
    fctr.getProperties().push_back(prop);
  }
  assert( (int) fctr.size() == nfprps );
  
  //std::cout << "Loaded face properties." << std::endl;
  
  // Load object properties
  assert( octr.size() == 0 );    
  for( int i = 0; i < (int) noprps; ++i )
  {
    PropertyBase* prop = NULL;
    loadProperty( ifs, &prop, *tpobj, 1 );
    octr.getProperties().push_back(prop);
  }
  assert( (int) octr.size() == noprps );
}

void TopologicalObjectSerializer::reloadRodFromFile( ElasticRod* rod, const std::string& filename )
{
  std::cout << "Not implemented" << std::endl;
  exit(1);

//  assert( rod != NULL );
//
//  // Attempt to open the file for binary output
//  std::ifstream ifs(filename.c_str(),std::ios::binary);
//  if( !ifs.is_open() ) 
//  {
//    std::cerr << "\033[31;1mWARNING IN TOPOLOGICALOBJECTSERIALIZER:\033[m Failed to open " << filename << " for reading. No data will be loaded." << std::endl;
//    return;
//  }
//
//  ObjectID objtypid = UNRECOGNIZEDOBJECTID;
//  ifs.read((char*)&objtypeid,sizeof(ObjectID));
//  assert( objtypid == RODID );
//
//  // ... pick up here
//
//  ifs.close();  
}






  
  
  
  
  
  
  
BASim::ObjectID TopologicalObjectSerializer::computeObjectId( const TopologicalObject& topobj ) const
{
  ObjectID objid;
  if(dynamic_cast<ElasticRod*>(const_cast<TopologicalObject*>(&topobj))) 
  {
    objid = RODID;
  }
  else if(dynamic_cast<TriangleMesh*>(const_cast<TopologicalObject*>(&topobj)))
  {
    objid = TRIMESHID;
  }
  else
  {
    std::cerr << "\033[31;1mWARNING IN TOPOLOGICALOBJECTSERIALIZER:\033[m Attempt to serialize invalid object type. Output will be corrupt." << std::endl;
    objid = UNRECOGNIZEDOBJECTID;
  }
  return objid;
}

void TopologicalObjectSerializer::serializeHeader( std::ofstream& of, const TopologicalObject& topobj ) const
{
  assert( of.is_open() );
  ObjectID objid = computeObjectId(topobj);    
  of.write((char*)&objid,sizeof(ObjectID));
}

void TopologicalObjectSerializer::loadHeader( std::ifstream& ifs, TopologicalObject** tpobj, ObjectID& objtypeid )
{
  assert( ifs.is_open() );
  assert( *tpobj == NULL );

  ifs.read((char*)&objtypeid,sizeof(ObjectID));

  switch (objtypeid) 
  {
    case RODID:
    {
      *tpobj = new AnisotropicRod;
      break;
    }
    case TRIMESHID:
    {
      *tpobj = new TriangleMesh;
      break;
    }
    default:
    {
      std::cerr << "\033[31;1mERROR IN TOPOLOGICALOBJECTSERIALIZER:\033[m Invalid object id encountered." << std::endl;
      exit(1);
    }
  }
}
  
void TopologicalObjectSerializer::serializeNumProperties( std::ofstream& of, const TopologicalObject& topobj ) const
{
  assert( of.is_open() );
  
  // Dump number of object properties
  int numobjprops = getObjectPropertyContainer(topobj).size();
  assert( numobjprops >= 0 );
  of.write((char*)&numobjprops,sizeof(int));
  
  // Dump number of vertex properties
  int numvertexprops = getVertexPropertyContainer(topobj).size();
  assert( numvertexprops >= 0 );
  of.write((char*)&numvertexprops,sizeof(int));
  
  // Dump number of edge properties
  int numedgeprops = getEdgePropertyContainer(topobj).size();
  assert( numedgeprops >= 0 );
  of.write((char*)&numedgeprops,sizeof(int));
  
  // Dump number of face properties
  int numfaceprops = getFacePropertyContainer(topobj).size();
  assert( numfaceprops >= 0 );
  of.write((char*)&numfaceprops,sizeof(int));
}

void TopologicalObjectSerializer::loadNumProperties( std::ifstream& ifs, int& noprps, int& nvprps, int& neprps, int& nfprps ) const
{
  assert( ifs.is_open() );

  // Load number of object properties
  ifs.read((char*)&noprps,sizeof(int));
  assert( noprps >= 0 );
  
  // Load number of vertex properties
  ifs.read((char*)&nvprps,sizeof(int));
  assert( nvprps >= 0 );
  
  // Load number of edge properties
  ifs.read((char*)&neprps,sizeof(int));
  assert( neprps >= 0 );
  
  // Load number of face properties
  ifs.read((char*)&nfprps,sizeof(int));
  assert( nfprps >= 0 );
}

void TopologicalObjectSerializer::serializePropertySizes( std::ofstream& of, const TopologicalObject& topobj ) const
{
  assert( of.is_open() );

  // Don't write the object property size, its always 1
  int sze = getObjectPropertyContainer(topobj).getProperties()[0]->size();
  assert( sze == 1 );

  // Write the vertex property size
  sze = getVertexPropertyContainer(topobj).getProperties()[0]->size();
  of.write((char*)&sze,sizeof(int));

  // Write the edge property size
  sze = getEdgePropertyContainer(topobj).getProperties()[0]->size();
  of.write((char*)&sze,sizeof(int));

  // Write the face property size
  sze = getFacePropertyContainer(topobj).getProperties()[0]->size();
  of.write((char*)&sze,sizeof(int));
}

void TopologicalObjectSerializer::loadPropertySizes( std::ifstream& ifs, int& nv, int& ne, int& nf ) const
{
  assert( ifs.is_open() );
  
  ifs.read((char*)&nv,sizeof(int));
  assert( nv >= 0 );
  
  ifs.read((char*)&ne,sizeof(int));
  assert( ne >= 0 );
  
  ifs.read((char*)&nf,sizeof(int));
  assert( nf >= 0 );
}

const PropertyContainer& TopologicalObjectSerializer::getObjectPropertyContainer( const TopologicalObject& rod ) const
{
  ObjPropHandle<int> oh;
  return rod.container(oh);
}

PropertyContainer& TopologicalObjectSerializer::getObjectPropertyContainer( TopologicalObject& rod )
{
  ObjPropHandle<int> oh;
  return rod.container(oh);
}

const PropertyContainer& TopologicalObjectSerializer::getVertexPropertyContainer( const TopologicalObject& rod ) const
{
  VPropHandle<int> vh;
  return rod.container(vh);
}

PropertyContainer& TopologicalObjectSerializer::getVertexPropertyContainer( TopologicalObject& rod )
{
  VPropHandle<int> vh;
  return rod.container(vh);
}  

const PropertyContainer& TopologicalObjectSerializer::getEdgePropertyContainer( const TopologicalObject& rod ) const
{
  EPropHandle<int> eh;
  return rod.container(eh);
}

PropertyContainer& TopologicalObjectSerializer::getEdgePropertyContainer( TopologicalObject& rod )
{
  EPropHandle<int> eh;
  return rod.container(eh);
}  

const PropertyContainer& TopologicalObjectSerializer::getFacePropertyContainer( const TopologicalObject& rod ) const
{
  FPropHandle<int> fh;
  return rod.container(fh);
}

PropertyContainer& TopologicalObjectSerializer::getFacePropertyContainer( TopologicalObject& rod )
{
  FPropHandle<int> fh;
  return rod.container(fh);
}

} // namespace BASim


