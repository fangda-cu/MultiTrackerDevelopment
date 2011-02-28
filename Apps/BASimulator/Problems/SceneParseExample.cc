/**
 * \file SceneParseExample.cc
 *
 * \author smith@cs.columbia.edu
 * \date 04/18/2010
 */

#include "SceneParseExample.hh"
#include <vector>

#include <stdio.h>
#include <iostream>
#include <cassert>
#include <string>
#include <fstream>


SceneParseExample::SceneParseExample()
: Problem("Scene Parse Example", "An example of how to load a scene from an xml file.")
{
  addDynamicsProps();
  addRodOptions();
  addRodTimeStepperOptions();
  
  GetIntOpt("nv") = 13;
  GetScalarOpt("dt") = 0.001;
}

SceneParseExample::~SceneParseExample()
{
}

void SceneParseExample::Setup()
{
  loadDynamicsProps();
  
  RodOptions opts;
  getRodOptions(opts);
  
  std::vector<Vec3d> undeformed;
  for( int i = 0; i < opts.numVertices; ++i ) undeformed.push_back(Vec3d((double) i,0.0,0.0));
  
  Vec3d up(1.0/sqrt(2.0),1.0/sqrt(2.0),0.0);
  Vec3d down(1.0/sqrt(2.0),-1.0/sqrt(2.0),0.0);
  std::vector<Vec3d> vertices;
  vertices.push_back(Vec3d(0.0,0.0,0.0));
  for( int i = 1; i < opts.numVertices; ++i ) 
  {
    if( i%2 == 1 ) vertices.push_back( vertices[i-1]+up );
    else vertices.push_back( vertices[i-1]+down );
  }
  
  ElasticRod* newrod = setupRod(opts, vertices, undeformed);
  m_world->addObject(newrod);
  RodTimeStepper* newstepper = getRodTimeStepper(*newrod);
  m_world->addController(newstepper);
  
  for( int i = 0; i < (int) vertices.size(); ++i ) vertices[i] += Vec3d(0.0,1.0,0.0);
  newrod = setupRod(opts, vertices, undeformed);
  m_world->addObject(newrod);
  newstepper = getRodTimeStepper(*newrod);
  m_world->addController(newstepper);
  
  ObjParser objparser;
  TriangleMesh* tri_mesh = new TriangleMesh();
  objparser.loadTriangularMesh( "assets/TriangulatedTorus.obj", *tri_mesh );
  m_world->addObject(tri_mesh);
  
  
  // Play with the scene loader
//  XMLSceneParser sceneparser( "assets/examplexmlscene.xml" );
//
//  std::vector<std::string> rodfilenames;
//  sceneparser.loadRodFileNames( rodfilenames );
//  
//  for( int i = 0; i < (int) rodfilenames.size(); ++i ) std::cout << rodfilenames[i] << std::endl;

  // Play with the xml parser
//  RodTextFileParser testparser("assets/rod0.txt");
//
//  std::cout << "nv: " << testparser.getNumVerts() << std::endl;
//  std::vector<Vec3d> deformed;
//  testparser.getDeformedVertices(deformed);
//  for( int i = 0; i < (int) deformed.size(); ++i ) std::cout << deformed[i] << std::endl;

  // Play with the xml scene outputter
}

void SceneParseExample::AtEachTimestep()
{
  XMLSceneOutputter sceneoutputter;
  sceneoutputter.outputScene( "output", "scene.xml", *this );
}
