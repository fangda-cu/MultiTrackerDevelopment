#include "ObjectSerializer.hh"

namespace BASim 
{

void ObjectSerializer::saveObject( const ObjectBase& obj, const std::string& filename )
{
  // I've never run this method!
  
  // Attempt to open the file for binary output
  std::ofstream of("example.bin",std::ios::binary);
  if( !of.is_open() ) 
  {
    std::cerr << "\033[31;1mWARNING IN OBJECTSERIALIZER:\033[m Failed to open " << filename << " for output. No data will be saved." << std::endl;
    return;
  }
  appendObjectToFile( obj, of );
  
  of.close();
}

void ObjectSerializer::appendObjectToFile( const ObjectBase& obj, std::ofstream& of )
{
  //std::cout << "ObjectSerializer::appendObjectToFile() (should be world object)" << std::endl;

  assert( of.is_open() );

  const PropertyContainer::Properties& oprops = obj.getPropertyContainer().getProperties();
  
  // Dump number of object properties
  int numobjprops = (int)( oprops.size() );
  serializeVal(of,numobjprops);

  // Write the object properties
  for( PropertyContainer::Properties::size_type i = 0; i < oprops.size(); ++i ) BASim::serializeProperty( of, oprops[i] );
}
  
void ObjectSerializer::loadObject( ObjectBase** obj, const std::string& filename )
{
  // I've never run this method!
  assert( *obj == NULL );

  // Attempt to open the file for binary output
  std::ifstream ifs("example.bin",std::ios::binary);
  if( !ifs.is_open() ) 
  {
    std::cerr << "\033[31;1mWARNING IN OBJECTSERIALIZER:\033[m Failed to open " << filename << " for reading. No data will be loaded." << std::endl;
    return;
  }
  loadObjectFromFile( obj, ifs );

  ifs.close();
}

void ObjectSerializer::loadObjectFromFile( ObjectBase** obj, std::ifstream& ifs )
{
  assert( ifs.is_open() );
  
  // Load the number of object properties
  int numprops = -1;
  loadVal(ifs,numprops);
  assert( numprops >= 0 );

  std::cout << numprops << std::endl;
  
  // Note: PropertyContainer::size() gives the size of the container, NOT the size of the properties.
  // Load object properties
  PropertyContainer::Properties& oprops = (*obj)->getPropertyContainer().getProperties();
  oprops.clear();
  assert( oprops.size() == 0 );    
  for( int i = 0; i < numprops; ++i )
  {
    PropertyBase* prop = NULL;
    loadProperty( ifs, &prop, *obj, 1 );
    std::cout << prop->name() << std::endl;
    oprops.push_back(prop);
  }
  assert( (int) oprops.size() == numprops );
}

} // namespace BASim


