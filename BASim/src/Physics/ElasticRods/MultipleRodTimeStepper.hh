/**
 * \file MultipleRodTimeStepper.hh
 *
 * \author smith@cs.columbia.edu
 * \date 05/29/2010
 */

// TODO: Make computing forces faster by not allocating extra vectors and just passing a base index
// TODO: Update the documentation!
// TODO: Don't recompute all of the degrees of freedom everytime ndof is called...
// TODO: Get this working with statics solver

#ifndef MULTIPLERODTIMESTEPPER_HH
#define MULTIPLERODTIMESTEPPER_HH

#include "../../Core/ObjectControllerBase.hh"
#include "../../Math/TimeSteppingBase.hh"
#include "../../Math/SolverUtils.hh"
#include "RodBoundaryCondition.hh"
#include "RodExternalForce.hh"
#include "RodRodExternalForce.hh"
#include "MinimalRodStateBackup.hh"

namespace BASim
{

/** 
 * Steps multiple rods, with inter-rod forces (e.g. springs connecting rods), forward in time.
 */
class MultipleRodTimeStepper : public ObjectControllerBase
{
public:

  //enum Method { SYMPL_EULER, IMPL_EULER, STATICS, NONE };
  enum Method { SYM_IMPL_EULER, SYMPL_EULER, IMPL_EULER, NONE };

  MultipleRodTimeStepper();
  ~MultipleRodTimeStepper();

  bool execute()
  {
    assert( m_diffEqSolver != NULL );
    return m_diffEqSolver->execute();
  }

  void setTime( Scalar time )
  {
    assert( m_diffEqSolver != NULL );
    m_diffEqSolver->setTime(time);
  }

  Scalar getTime() const
  {
    assert( m_diffEqSolver != NULL );
    return m_diffEqSolver->getTime();
  }

  void setTimeStep( Scalar dt );

  Scalar getTimeStep() const
  {
    assert( m_diffEqSolver != NULL );
    return m_diffEqSolver->getTimeStep();
  }

  const DiffEqSolver& getDiffEqSolver() const
  {
    assert( m_diffEqSolver != NULL );
    return *m_diffEqSolver;
  }

  void setDiffEqSolver( Method method );

  int ndof() const;

  /**
   * Returns all "masses" in this differential equation in a flat vector.
   *
   * \param[out] masses Vector that will contain masses. Must be the proper size.
   */
  void getMass( VecXd& masses );
  
  void setX( const VecXd& positions );
  
  /**
   * Returns all "positions" in this differential equation in a flat vector.
   *
   * \param[out] masses Vector that will contain positions. Must be the proper size.
   */
  void getX( VecXd& positions ) const;
  
  void setV( const VecXd& velocities );
  
  /**
   * Returns all "velocities" in this differential equation in a flat vector.
   *
   * \param[out] masses Vector that will contain velocities. Must be the proper size.
   */
  void getV( VecXd& velocities ) const;
  
  /**
   * This function computes the force on each degree of freedom
   * associated to the rod.
   *
   * \param[out] f The vector of accelerations on the rod.
   */
  void evaluatePDot( VecXd& f );

  /**
   * Evaluates the Jacobian of the forces on the rod.
   *
   * \param[out] J The Jacobian of the forces on the rod.
   */
  void evaluatePDotDX( Scalar scale, MatrixBase& J );

  void evaluatePDotDV( Scalar scale, MatrixBase& J );

  /**
   * This function returns the mass associated with a degree of
   * freedom of the rod.
   *
   * \param[in] i Which degree of freedom.
   * \return The mass of the degree of freedom
   */
//  Scalar getMass( int i )
//  {
//    assert( false );
//
//    //assert(i >= 0);
//    //assert(i < m_rod.ndof());
//
//    //return m_rod.getMass(i);
//    
//    return std::numeric_limits<double>::signaling_NaN();
//  }

//  Scalar getX( int i ) const
//  {
//    assert( false );
//
//    //return m_rod.getDof(i);
//
//    return std::numeric_limits<double>::signaling_NaN();
//  }
//
//  void setX( int i, Scalar x )
//  {
//    assert( false );
//
//    //m_rod.setDof(i, x);
//  }
//
//  Scalar getV( int i )
//  {
//    assert( false );
//
//    //return m_rod.getVel(i);
//
//    return std::numeric_limits<double>::signaling_NaN();
//  }
//
//  void setV( int i, Scalar v )
//  {
//    assert( false );
//
//    //m_rod.setVel(i, v);
//  }

  void increment_q( const VecXd& dq );

  void set_q(const VecXd& q);

  // TODO: Double check the indexing here. 
  // TODO: Add some sanity checks?
  void increment_qdot( const VecXd& dqd );

  void set_qdot(const VecXd& qd);

  /**
   * Adds an external force to be applied to the rod. On destruction,
   * this class will be responsible for de-allocating the memory
   * associated to the force.
   *
   * \param[in] force The external force to be applied to rods.
   */
  void addExternalForce( RodExternalForce* force )
  {
    assert( force != NULL );
    m_externalForces.push_back(force);
  }

  void addRodRodExternalForce( RodRodExternalForce* rodRodForce )
  {
    assert( rodRodForce != NULL );
    m_rodRodexternalForces.push_back(rodRodForce);
  }
  
  void addRod( ElasticRod* rod )
  {
    assert( rod != NULL );
    m_rods.push_back(rod);
  }
  
  std::vector<RodExternalForce*>& getExternalForces()
  {
    return m_externalForces;
  }

  void startStep();

  void endStep();

  void startIteration() {}

  void endIteration();

  int getMaxIterations() const
  {
    assert( m_diffEqSolver != NULL );
    return m_diffEqSolver->getMaxIterations();
  }

  void setMaxIterations( int iterations )
  {
    assert( m_diffEqSolver != NULL );
    m_diffEqSolver->setMaxIterations(iterations);
  }

  Scalar get_stol() const { return m_diffEqSolver->get_stol(); }
  void set_stol(Scalar s) { m_diffEqSolver->set_stol(s); }

  Scalar get_atol() const { return m_diffEqSolver->get_atol(); }
  void set_atol(Scalar a) { m_diffEqSolver->set_atol(a); }

  Scalar get_rtol() const { return m_diffEqSolver->get_rtol(); }
  void set_rtol(Scalar r) { m_diffEqSolver->set_rtol(r); }

  Scalar get_inftol() const { return m_diffEqSolver->get_inftol(); }
  void set_inftol(Scalar i) { m_diffEqSolver->set_inftol(i); }

  MatrixBase* createMatrix() const;

  RodBoundaryCondition* getBoundaryCondition( int rod_num )
  {
    assert( rod_num >= 0 );
    assert( rod_num < (int) m_rods.size() );
    
    return m_rods[rod_num]->getBoundaryCondition();
  }

  // TODO: Does anyone use this method? It seems kind of dangerous letting any user
  //       mess with our internal pointers :)
  void setBoundaryCondition(RodBoundaryCondition* bc)
  {
    assert( false );
    //m_boundaryCondition = bc;
  }

  void getScriptedDofs( IntArray& indices, std::vector<Scalar>& desired );

  std::vector<ElasticRod*>& getRods()
  {
    return m_rods;
  }



  void backup()
  {
    //m_backupstate.backupRod(m_rod);
    for( size_t i = 0; i < m_rods.size(); ++i ) m_backupstate[i].backupRod(*m_rods[i]);
  }

  void backupResize()
  {
    //m_backupstate.resize(m_rod);
    m_backupstate.resize(m_rods.size());
    for( size_t i = 0; i < m_rods.size(); ++i ) m_backupstate[i].resize(*m_rods[i]);
  }

  void backupRestore()
  {
    //m_backupstate.restoreRod(m_rod);
    for( size_t i = 0; i < m_rods.size(); ++i ) m_backupstate[i].restoreRod(*m_rods[i]);
  }

  void backupClear()
  {
    //m_backupstate.clear();
    for( size_t i = 0; i < m_rods.size(); ++i ) m_backupstate[i].clear();
    m_backupstate.clear();
  }


protected:

  Method m_method;
  std::vector<ElasticRod*> m_rods;
  std::vector<RodExternalForce*> m_externalForces;
  DiffEqSolver* m_diffEqSolver;

  std::vector<RodRodExternalForce*> m_rodRodexternalForces;
  
  std::vector<MinimalRodStateBackup> m_backupstate;

};

} // namespace BASim

#endif // RODTIMESTEPPER_HH
