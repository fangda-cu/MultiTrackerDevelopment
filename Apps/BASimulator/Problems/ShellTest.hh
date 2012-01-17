/**
 * \file BentTwisting.hh
 *
 * \author batty@cs.columbia.edu
 * \date April 20, 2011
 */

#ifndef SHELLTEST_HH
#define SHELLTEST_HH

#include "ProblemBase.hh"

#include "BASim/src/Physics/DeformableObjects/DeformableObject.hh"
#include "BASim/src/Physics/DeformableObjects/Shells/ElasticShell.hh"
#include "BASim/src/Physics/DeformableObjects/DefoObjTimeStepper.hh"
/**
 * Test problem for the shells...
 */

class ShellTest : public Problem
{
public:

  ShellTest();
  virtual ~ShellTest();

  std::ofstream sphere_data;

  //virtual void serialize( std::ofstream& of );
  //virtual void resumeFromfile( std::ifstream& ifs );

protected:

  void Setup();
  void AtEachTimestep();

  DeformableObject* shellObj;
  ElasticShell* shell;
  DefoObjTimeStepper* stepper;

  Scalar m_timestep;
  Scalar m_initial_thickness;
  int m_active_scene;

public:
  void setupScene1(); //vertical flat sheet
  void setupScene2(); //vertical cylindrical sheet
  void setupScene3(); //spherical sheet
  void setupScene4(); //two-triangle bending test
  void setupScene5(); //catenary 
  void setupScene6(); //viscous hemispherical bubble
  void setupScene7(); //sheet between two circles
  void setupScene8(); //torus
  void setupScene9(); //non-manifold edge example
  void setupScene10(); //flat stretching sheet with inflow/deletion
  void setupScene11(); //a cube
  void setupScene12(); //popping low viscosity hemispherical bubble with Bernoulli air pressure model
  void setupScene13(); //a constant inflow hitting a solid boundary
  void setupScene14(); //a rectangular sheet falling onto a sphere
  void setupScene15(); //Two annular rings spreading apart
  void setupScene16(); //a pouring sheet falling onto conveyor belt with objects(sphere)
  void setupScene17(); //a sheet falling onto a static object defined by a signed distance field
  void setupScene18(); //a constant inflow hitting a solid boundary
};

#endif // SHELLTEST_HH

