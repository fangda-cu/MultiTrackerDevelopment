#include "Problems/BentTwisting.hh"
#include "Problems/SHO.hh"
#include "Problems/CollisionTest.hh"
#include "Problems/HairyBall.hh"
#include "Problems/MayaSceneTest.hh"


void CreateProblemVector()
{
  problems.push_back(NULL);
  problems.push_back(new BentTwisting());
  problems.push_back(new SHO());
  problems.push_back(new CollisionTest());
  problems.push_back(new HairyBall());
  problems.push_back(new MayaSceneTest());
}
