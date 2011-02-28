/**
 * \file SerializationTests.hh
 *
 * \author smith@cs.columbia.edu
 * \date 08/30/2010
 */

#ifndef SERIALIZATIONTESTS_HH
#define SERIALIZATIONTESTS_HH

#include <iostream>
#include <fstream>

#include "ProblemBase.hh"
#include "BASim/src/Physics/ElasticRods/TopologicalObjectSerializer.hh"

class SerializationTests : public Problem
{
public:
  
  SerializationTests();
  virtual ~SerializationTests();
  
  void serialize( std::ofstream& of );
  void resumeFromfile( std::ifstream& ifs );
  
protected:

  void setupTwoVertexNoMovement();
  void setupTwoVertexConstantVelocityNoForcesExcited();
  
  void Setup();

private:
  ElasticRod* m_rod;
  RodTimeStepper* m_stepper;
};

#endif 

