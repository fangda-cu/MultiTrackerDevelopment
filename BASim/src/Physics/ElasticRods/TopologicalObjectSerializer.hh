/**
 * \file TopologicalObjectSerializer.hh
 *
 * \author smith@cs.columbia.edu
 * \date 07/16/2010
 */

// DONT USE THIS CLASS YET STUFF IS STILL CHANGING A BUNCH. IT SHOULD BE READY
// VERY SOON MY FRIENDS. B.S. 07/16/2010

#ifndef TOPOLOGICALOBJECTSERIALIZER_HH
#define TOPOLOGICALOBJECTSERIALIZER_HH

#include <fstream>

#ifdef WETA
#include "../..//Core/Property.hh"
#include "ElasticRod.hh"
#include "RodStretchingForce.hh"
#include "RodTwistingForceSym.hh"
#include "RodBendingForceSym.hh"
#include "../../IO/SerializationUtils.hh"
#include "../../Core/TriangleMesh.hh"
#include "AnisotropicRod.hh"
#else
#include "BASim/src/Core/Property.hh"
#include "BASim/src/Physics/ElasticRods/ElasticRod.hh"
#include "BASim/src/Physics/ElasticRods/RodStretchingForce.hh"
#include "BASim/src/Physics/ElasticRods/RodTwistingForceSym.hh"
#include "BASim/src/Physics/ElasticRods/RodBendingForceSym.hh"
#include "BASim/src/IO/SerializationUtils.hh"
#endif

namespace BASim 
{

// TODO: When loading properties, I probably have an unneeded set of copies. Clean this up.
// TODO: Clean up handling of rods vs. cloth vs. etc.
// TODO: Handle errors in a consistent manner.
// TODO: Could add extra error checks to input/output.
// TODO: Rename this all, a lot is actually specific to rods

class TopologicalObjectSerializer
{
public:

  // These are fairly generic
  void saveTopologicalObject( const TopologicalObject& topobj, const std::string& filename );
  void appendTopologicalObjectToFile( const TopologicalObject& topobj, std::ofstream& of );

  // These actually have to do some rod-specific stuff, rename at some point
  void loadTopologicalObject( ElasticRod** tpobj, const std::string& filename );
  void loadTopologicalObjectFromFile( ElasticRod** tpobj, std::ifstream& ifs );

  void reloadRodFromFile( ElasticRod* rod, const std::string& filename );

  void loadGenericTopologicalObject( TopologicalObject** tpobj, std::ifstream& ifs );

private:

  // Determines the type of this topological object. 
  ObjectID computeObjectId( const TopologicalObject& topobj ) const;
  
  // Writes an object ID that identifies this specific topoological object as a rod, cloth, etc.
  void serializeHeader( std::ofstream& of, const TopologicalObject& topobj ) const;
  
  // Load the object id, creates the object.
  void loadHeader( std::ifstream& ifs, TopologicalObject** tpobj, ObjectID& objtypeid );
  
  // Writes the number of object, vertex, edge, and face properties present in this container.
  void serializeNumProperties( std::ofstream& of, const TopologicalObject& topobj ) const;

  void loadNumProperties( std::ifstream& ifs, int& noprps, int& nvprps, int& neprps, int& nfprps ) const;

  // Dumps number of vertices, edges, face.
  void serializePropertySizes( std::ofstream& of, const TopologicalObject& topobj ) const;
  
  void loadPropertySizes( std::ifstream& ifs, int& nv, int& ne, int& nf ) const;


  const PropertyContainer& getObjectPropertyContainer( const TopologicalObject& rod ) const;
  PropertyContainer& getObjectPropertyContainer( TopologicalObject& rod );
  
  const PropertyContainer& getVertexPropertyContainer( const TopologicalObject& rod ) const;
  PropertyContainer& getVertexPropertyContainer( TopologicalObject& rod );
  
  const PropertyContainer& getEdgePropertyContainer( const TopologicalObject& rod ) const;
  PropertyContainer& getEdgePropertyContainer( TopologicalObject& rod );
  
  const PropertyContainer& getFacePropertyContainer( const TopologicalObject& rod ) const;
  PropertyContainer& getFacePropertyContainer( TopologicalObject& rod );
};

} // namespace BASim

#endif // TOPOLOGICALOBJECTSERIALIZER_HH

