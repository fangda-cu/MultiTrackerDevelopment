/**
 * \file ImplicitEulerTestsOne.hh
 *
 * \author smith@cs.columbia.edu
 * \date 05/03/2010
 */

#ifndef IMPLICITEULERTESTS_HH
#define IMPLICITEULERTESTS_HH

#include "ProblemBase.hh"

class ImplicitEulerTestsOne : public Problem
{
public:
  
  ImplicitEulerTestsOne();
  virtual ~ImplicitEulerTestsOne();
  
protected:

  void Setup();
  void AtEachTimestep();

private:
  std::vector<ElasticRod*> m_rods;
};

#endif 

