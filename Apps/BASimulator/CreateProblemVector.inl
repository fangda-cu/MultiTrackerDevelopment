#include "Problems/BentTwisting.hh"
#include "Problems/SHO.hh"
#include "Problems/ViscosityTests.hh"
#include "Problems/RegressionTests.hh"
#include "Problems/HairyBall.hh"
#include "Problems/TestCase6.hh"
//#include "Problems/PlantRootGrowth.hh"
//#include "Problems/HandTest.hh"

//#include "Problems/SerializationTests.hh"

#include "Problems/CollisionTestTwo.hh"
#include "Problems/CollisionTestWeta.hh"
#include "Problems/HairyBallWeta.hh"

#ifdef HAVE_PARDISO
//#include "Problems/MicrotuboliDNA.hh"
#endif

void CreateProblemVector()
{
  problems.push_back(NULL);
  problems.push_back(new BentTwisting());

  problems.push_back(new SHO());
  problems.push_back(new ViscosityTests());
  problems.push_back(new RegressionTests());
  problems.push_back(new HairyBall());
  problems.push_back(new TestCase6());

  //problems.push_back(new PlantRootGrowth());
  //problems.push_back(new HandTest());

 // problems.push_back(new SerializationTests());
  
//  problems.push_back(new CollisionTestTwo());
//  problems.push_back(new CollisionTestWeta());
 // problems.push_back(new HairyBallWeta());
 
  #ifdef HAVE_PARDISO
  //problems.push_back(new MicrotuboliDNA());
  #endif
}
