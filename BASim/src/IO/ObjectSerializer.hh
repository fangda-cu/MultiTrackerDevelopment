/**
 * \file ObjectSerializer.hh
 *
 * \author smith@cs.columbia.edu
 * \date 07/16/2010
 */

// TODO: Test test test! 

#ifndef OBJECTSERIALIZER_HH
#define OBJECTSERIALIZER_HH

#include <fstream>

#ifdef WETA
#include "../Core/Property.hh"
#include "../Core/ObjectBase.hh"
#include "../Physics/ElasticRods/RodStretchingForce.hh"
#include "../Physics/ElasticRods/RodTwistingForceSym.hh"
#include "../Physics/ElasticRods/RodBendingForceSym.hh"
#include "../IO/SerializationUtils.hh"
#else
#include "BASim/src/Core/Property.hh"
#include "BASim/src/Core/ObjectBase.hh"
//#include "BASim/src/Physics/ElasticRods/RodStretchingForce.hh"
//#include "BASim/src/Physics/ElasticRods/RodTwistingForceSym.hh"
//#include "BASim/src/Physics/ElasticRods/RodBendingForceSym.hh"
#include "BASim/src/IO/SerializationUtils.hh"
#endif

namespace BASim 
{

// TODO: When loading properties, I probably have an unneeded set of copies. Clean this up.
// TODO: Clean up handling of rods vs. cloth vs. etc.
// TODO: Handle errors in a consistent manner.
// TODO: Consolidate a number of the serializing methods, they are basically the same.
// TODO: Could add extra error checks to input/output.

class ObjectSerializer
{
public:

  void saveObject( const ObjectBase& obj, const std::string& filename );

  void appendObjectToFile( const ObjectBase& obj, std::ofstream& of );

  void loadObject( ObjectBase** obj, const std::string& filename );

  void loadObjectFromFile( ObjectBase** obj, std::ifstream& ifs );
};

}

#endif

