#include "Problems/BentTwisting.hh"
#include "Problems/SHO.hh"

void CreateProblemVector()
{
  problems.push_back(NULL);
  problems.push_back(new BentTwisting());
  problems.push_back(new SHO());
}
