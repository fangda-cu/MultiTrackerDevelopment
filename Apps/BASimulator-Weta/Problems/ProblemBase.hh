/**
 * \file ProblemBase.hh
 *
 * \author miklos@cs.columbia.edu (based on problem-setup.h from sar2120@columbia.edu)
 * \date 09/09/2009
 */

#ifndef PROBLEMBASE_HH
#define PROBLEMBASE_HH

#include "Option.hh"

#include <map>
#include <string>
#include <queue>
#include <BASim/BASim>

#include "BASim/src/Physics/ElasticRods/MultipleRodTimeStepper.hh"

#include "BASim/src/IO/ObjectSerializer.hh"

namespace BASim
{

class Problem
{
public:

  //typedef std::map<std::string, Option> OptionMap;

  Problem(const std::string& name = "", const std::string& desc = "");
  virtual ~Problem();

  void BaseSetup(int argc, char** argv);
  void BaseFinalize();
  void BaseAtEachTimestep();

  const std::string& ProblemName() const { return m_problemName; }
  const std::string& ProblemDescription() const { return m_problemDesc; }

  template <typename T>
  int AddOption(const std::string& name, const std::string& desc, const T& def);

  Option* GetOption(const std::string& name);
  bool& GetBoolOpt(const std::string& name);
  int& GetIntOpt(const std::string& name);
  Scalar& GetScalarOpt(const std::string& name);
  Vec3d& GetVecOpt(const std::string& name);
  std::string& GetStringOpt(const std::string& name);

  int LoadOptions(const char* filename);
  int LoadOptions(const std::string& filename)
  {
    return LoadOptions(filename.c_str());
  }
  int LoadOptions(int argc, char** argv);

  void PrintOptions(std::ostream& os);

  const World& getWorld() const { return *m_world; }
  World& getWorld() { return *m_world; }

  /**
   * Adds properties to the World that are generally used in dynamics
   * problems.
   */
  void addDynamicsProps();

  /**
   * Loads in dynamics properties from the options. setupDynamics must
   * be called before this can be used.
   */
  void loadDynamicsProps();

  bool dynamicsPropsLoaded() const { return m_dynamicsProps; }

  const Scalar& getTime() const;
  void setTime(const Scalar& time);

  const Scalar& getDt() const;
  void setDt(const Scalar& dt);

  const Vec3d& getGravity() const;
  void setGravity(const Vec3d& gravity);

  void addRodOptions();
  void getRodOptions(RodOptions& opts);
  void addRodTimeStepperOptions();
  RodTimeStepper* getRodTimeStepper(ElasticRod& rod);
  MultipleRodTimeStepper* getMultipleRodTimeStepper();

  std::queue<double>& getBreakpoints() { return m_sim_breakpoints; }
  void insertBreakpoint( double t ) { m_sim_breakpoints.push(t); }

  // Children that support serialization must override this method
  virtual void serialize( std::ofstream& of )
  {
    std::cerr << "\033[31;1mERROR IN PROBLEMBASE:\033[m Serialization not implemented for problem: " << m_problemName << ". No output saved. Exiting." << std::endl;
    exit(1);
  }
  
  // Children that support serialization must override this method
  virtual void resumeFromfile( std::ifstream& ifs )
  {
    std::cerr << "\033[31;1mERROR IN PROBLEMBASE:\033[m Serialization not implemented for problem: " << m_problemName << ". Not resuming problem. Exiting." << std::endl;
    exit(1);
  }

  void serializeProblem( std::ofstream& of )
  {
    assert( of.is_open() );
    assert( m_world != NULL );

    // Serialize the world class
    int WORLDMAGIC = 0xE0E0E0E0;
    of.write((char*)&WORLDMAGIC,sizeof(int));
    of.write((char*)&WORLDMAGIC,sizeof(int));

    ObjectSerializer objserializer;
    objserializer.appendObjectToFile(*m_world,of);

    of.write((char*)&WORLDMAGIC,sizeof(int));
    of.write((char*)&WORLDMAGIC,sizeof(int));

    int PROBLEMBASEMAGIC = 0xD0D0D0D0;
    of.write((char*)&PROBLEMBASEMAGIC,sizeof(int));
    of.write((char*)&PROBLEMBASEMAGIC,sizeof(int));

    // Serialize this class
    //std::cout << "Problem::serializeProblem()" << std::endl;

    // Serialize the remaining state from this problem
    // std::string m_problemName < taken care of by constructor
    // std::string m_problemDesc < taken care of by constructor
    serializeMapStringOption(of,m_options);
    // World* m_world < taken care of previously
    serializeBool(of,m_dynamicsProps);
    // ObjPropHandle<Scalar> m_time < taken care of by World
    // ObjPropHandle<Scalar> m_dt < taken care of by World
    // ObjPropHandle<Vec3d> m_gravity < taken care of by World
    serializeDoubleQueue(of,m_sim_breakpoints);
    
    of.write((char*)&PROBLEMBASEMAGIC,sizeof(int));
    of.write((char*)&PROBLEMBASEMAGIC,sizeof(int));
  }
  
  void resumeProblem( std::ifstream& ifs )
  {
//    std::cout << "Loading from resumeProblem" << std::endl;  
    assert( ifs != NULL );
    assert( m_world != NULL );
    
    int byteeater;
    ifs.read((char*)&byteeater,sizeof(int));
    ifs.read((char*)&byteeater,sizeof(int));
    
    // Re-load the world object
    ObjectSerializer objserializer;
    ObjectBase* toload = m_world;
    objserializer.loadObjectFromFile( &toload, ifs );

    ifs.read((char*)&byteeater,sizeof(int));
    ifs.read((char*)&byteeater,sizeof(int));

    ifs.read((char*)&byteeater,sizeof(int));
    ifs.read((char*)&byteeater,sizeof(int));

    PropertyContainer::Properties& oprops = m_world->getPropertyContainer().getProperties();
    int numprops = oprops.size();
    
    //std::cout << numprops << std::endl;
    //for( int i = 0; i < (int) oprops.size(); ++i ) std::cout << oprops[i]->name() << std::endl;

    loadMapStringOption(ifs,m_options);
    loadBool(ifs,m_dynamicsProps);
    loadDoubleQueue(ifs,m_sim_breakpoints);

//    this->PrintOptions(std::cout);
//    std::cout << "m_dynamicsProps: " << m_dynamicsProps << std::endl;
//    std::cout << "num breakpoints: " << m_sim_breakpoints.size() << std::endl;
//    std::cout << "    " << m_sim_breakpoints.front() << std::endl;

    if( m_dynamicsProps )
    {
      assert( m_world->property_exists(m_time,"time") );
      assert( m_world->property_exists(m_dt,"time-step") );
      assert( m_world->property_exists(m_gravity,"gravity") );

      m_world->property_handle(m_time,"time");
      m_world->property_handle(m_dt,"time-step");
      m_world->property_handle(m_gravity,"gravity");
    }

//    std::cout << "m_time: " << m_world->property(m_time) << std::endl;
//    std::cout << "m_dt: " << m_world->property(m_dt) << std::endl;
//    std::cout << "m_gravity: " << m_world->property(m_gravity).transpose() << std::endl;

    ifs.read((char*)&byteeater,sizeof(int));
    ifs.read((char*)&byteeater,sizeof(int));
  }


protected:

  virtual void Setup() {}
  virtual void AtEachTimestep() {}
  virtual void AfterLoad() {}
  virtual void AfterStep() {}

  std::string m_problemName;
  std::string m_problemDesc;
  std::map<std::string, Option> m_options;

  World* m_world;

  bool m_dynamicsProps;
  ObjPropHandle<Scalar> m_time;
  ObjPropHandle<Scalar> m_dt;
  ObjPropHandle<Vec3d> m_gravity;
  
  // Times to pause the simulation for debugging purpouses
  std::queue<double> m_sim_breakpoints;
};

}

#include "ProblemBase.tcc"

#endif // PROBLEMBASE_HH
