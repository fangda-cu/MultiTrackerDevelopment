/**
 * \file SceneParseExample.hh
 *
 * \author smith@cs.columbia.edu
 * \date 04/18/2010
 */

#ifndef SCENEPARSEEXAMPLE_HH
#define SCENEPARSEEXAMPLE_HH

#include "ProblemBase.hh"
#include "BASim/src/IO/XMLSceneParser.hh"
#include "BASim/src/IO/RodTextFileParser.hh"
#include "BASim/src/IO/XMLSceneOutputter.hh"

class SceneParseExample: public Problem
{
public:

  SceneParseExample();
  virtual ~SceneParseExample();

protected:

  void Setup();
  void AtEachTimestep();

};

#endif
