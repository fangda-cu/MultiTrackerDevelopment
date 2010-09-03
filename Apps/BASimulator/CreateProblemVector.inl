#include "Problems/BentTwisting.hh"
#include "Problems/CollisionTest.hh"
#include "Problems/HairyBall.hh"
#include "Problems/MayaSceneTest.hh"

//#include "Problems/CollisionTestTwo.hh"

#ifdef HAVE_PARDISO
#include "Problems/MicrotuboliDNA.hh"
#endif

void CreateProblemVector()
{
  problems.push_back(NULL);
  problems.push_back(new BentTwisting());
//  problems.push_back(new CollisionTest());
  problems.push_back(new HairyBall());
  problems.push_back(new MayaSceneTest());

  #ifdef HAVE_PARDISO
  problems.push_back(new MicrotuboliDNA());
  #endif
}
