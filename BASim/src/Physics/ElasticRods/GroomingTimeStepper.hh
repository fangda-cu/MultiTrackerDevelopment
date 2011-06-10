/**
 * \file GroomingTimeStepper.hh
 *
 * \author miklos@cs.columbia.edu
 * \date 09/03/2009
 */

#ifndef GROOMINGTIMESTEPPER_HH
#define GROOMINGTIMESTEPPER_HH

#ifdef WETA
#include "../..//Core/ObjectControllerBase.hh"
#include "RodBoundaryCondition.hh"
#include "../../Math/TimeSteppingBase.hh"
#include "../../Math/SolverUtils.hh"
#include "RodExternalForce.hh"
#include "RodExternalConservativeForce.hh"
#include "../../Math/SymplecticEuler.hh"
#include "../../Math/ImplicitEuler.hh"
#include "../../Math/StaticSolver.hh"
#include "../../Math/StaticSolver.hh"
#include "../../Physics/ElasticRods/MinimalRodStateBackup.hh"
#include "../../Util/TextLog.hh"
#else
#include "BASim/src/Core/ObjectControllerBase.hh"
#include "BASim/src/Physics/ElasticRods/RodBoundaryCondition.hh"
#include "BASim/src/Math/TimeSteppingBase.hh"
#include "BASim/src/Math/SolverUtils.hh"
#include "BASim/src/Physics/ElasticRods/RodExternalForce.hh"
#include "BASim/src/Physics/ElasticRods/RodExternalConservativeForce.hh"
#include "BASim/src/Math/SymplecticEuler.hh"
#include "BASim/src/Math/ImplicitEuler.hh"
#include "BASim/src/Math/StaticSolver.hh"
#include "BASim/src/Math/StaticSolver.hh"
#include "BASim/src/Util/TextLog.hh"
#include "BASim/src/Physics/ElasticRods/MinimalRodStateBackup.hh"
#endif

namespace BASim {

/** Class to time step a rod. */
class GroomingTimeStepper : public ObjectControllerBase
{
public:

  enum Method { SYMPL_EULER, IMPL_EULER, STATICS, NONE };

  explicit GroomingTimeStepper(ElasticRod& rod)
    : m_rod(rod)
    , m_method(NONE)
    , m_diffEqSolver(NULL)
    //, m_boundaryCondition(NULL)
    , m_backupstate()
  {
    setDiffEqSolver(SYMPL_EULER);
  }

  ~GroomingTimeStepper()
  {
    if (m_diffEqSolver != NULL) delete m_diffEqSolver;
    //if (m_boundaryCondition != NULL) delete m_boundaryCondition;
    for( int i = 0; i < (int) m_externalForces.size(); ++i )
    {
      assert( m_externalForces[i] != NULL );
      delete m_externalForces[i];
      m_externalForces[i] = NULL;
    }	
  }

  bool execute()
  {
    assert( getTimeStep() == m_rod.getTimeStep() );
    m_has_solved = m_diffEqSolver->execute();
    return m_has_solved;
  }

  bool HasSolved() { return m_has_solved; }

  void setTime(Scalar time)
  {
  //  std::cout << "Setting time in GroomingTimeStepper to be " << time  << std::endl;
    m_diffEqSolver->setTime(time);
  }

  Scalar getTime() const
  {
    double t = m_diffEqSolver->getTime();
    //std::cout << "RotTimeStepper::getTime() = " << t << std::endl;
    return m_diffEqSolver->getTime();
  }

  virtual void setTimeStep(Scalar dt)
  {
    m_diffEqSolver->setTimeStep(dt);
    m_rod.setTimeStep(dt);
  }

  virtual Scalar getTimeStep() const
  {
    return m_diffEqSolver->getTimeStep();
  }

  DiffEqSolver& getDiffEqSolver()
  {
    assert(m_diffEqSolver != NULL);

    return *m_diffEqSolver;
  }

  const DiffEqSolver& getDiffEqSolver() const
  {
    assert(m_diffEqSolver != NULL);

    return *m_diffEqSolver;
  }

  void setDiffEqSolver(Method method)
  {
    if (method == m_method) return;

    m_method = method;
    if (m_diffEqSolver != NULL) delete m_diffEqSolver;
    m_diffEqSolver = NULL;

    if (method == SYMPL_EULER)
    {
      m_diffEqSolver = new StaticSolver<GroomingTimeStepper>(*this);
    } 
    else if (method == IMPL_EULER)
    {
      //m_diffEqSolver = new ImplicitEuler<RodTimeStepper>(*this);
      m_diffEqSolver = new StaticSolver<GroomingTimeStepper>(*this);
      //m_diffEqSolver = new BacktrackingImplicitEuler<RodTimeStepper>(*this);
    }
    else if (method == STATICS) 
    {
      m_diffEqSolver = new StaticSolver<GroomingTimeStepper>(*this);
    } 
    else if (method == NONE)
    {
      m_diffEqSolver = NULL;
    }
    else
    {
      std::cout << "Unknown method specified" << std::endl;
      m_diffEqSolver = NULL;
    }

    if (m_diffEqSolver != NULL)
    {
      m_rod.setTimeStep(m_diffEqSolver->getTimeStep());
    }
  }

  int ndof() const
  {
    return m_rod.ndof();
  }

  /**
   * Returns all "masses" in this differential equation in a flat vector.
   *
   * \param[out] masses Vector that will contain masses. Must be the proper size.
   */
  void getMass( VecXd& masses )
  {
    assert( masses.size() == m_rod.ndof() );
    for( int i = 0; i < m_rod.ndof(); ++i ) masses(i) = m_rod.getMass(i);
  }

  void setX( const VecXd& positions )
  {
    assert( positions.size() == m_rod.ndof() );
    for( int i = 0; i < m_rod.ndof(); ++i ) m_rod.setDof(i, positions(i));
  }

  /**
   * Returns all "positions" in this differential equation in a flat vector.
   *
   * \param[out] masses Vector that will contain positions. Must be the proper size.
   */
  void getX( VecXd& positions ) const
  {
    assert( positions.size() == m_rod.ndof() );
    for( int i = 0; i < m_rod.ndof(); ++i ) positions(i) = m_rod.getDof(i);
  }
  
  void setV( const VecXd& velocities )
  {
    assert( velocities.size() == m_rod.ndof() );
    for( int i = 0; i < m_rod.ndof(); ++i ) m_rod.setVel(i, velocities(i));
  }
  
  /**
   * Returns all "velocities" in this differential equation in a flat vector.
   *
   * \param[out] masses Vector that will contain velocities. Must be the proper size.
   */
  void getV( VecXd& velocities ) const
  {
    assert( velocities.size() == m_rod.ndof() );
    for( int i = 0; i < m_rod.ndof(); ++i ) velocities(i) = m_rod.getVel(i);
  }

  /**
   * This function computes the force on each degree of freedom
   * associated to the rod.
   *
   * \param[out] f The vector of accelerations on the rod.
   */
  void evaluatePDot(VecXd& f)
  {
   // m_forces = f;
    // add internal forces

    // std::cout << "RodTimeStepper::evaluatePDot: rodidx = " << m_rod.globalRodIndex << "\n";
    // for (int i=0; i < m_rod.nv(); ++i)
    // {
    //   std::cout << "x[" << i << "] = " << m_rod.getVertex(i) << std::endl;
    // }

    m_rod.computeForces(f);

    //if (m_rod.viscous()) f /= m_diffEqSolver->getTimeStep();

    VecXd curr_force(f.size());

    // add external forces
    for (size_t i = 0; i < m_externalForces.size(); ++i) {
      curr_force.setZero();
      m_externalForces[i]->computeForce(m_rod, curr_force);
      f += curr_force;
      TraceStream(g_log, "") << m_externalForces[i]->getName() << " &rod = " << &m_rod << " norm = " << curr_force.norm() << '\n';
    }
  
   // m_forces = f - m_forces;
  }


  /**
   * This function computes the force on each degree of freedom
   * associated to the rod.
   *
   * \param[out] f The vector of accelerations on the rod.
   */
  void evaluateConservativeForcesEnergy(VecXd& f, Scalar& energy)
  {
    m_rod.computeConservativeForcesEnergy(f, energy);

    VecXd  curr_force(f.size());

    // add external conservative forces
    for (size_t i = 0; i < m_externalForces.size(); ++i) {
      curr_force.setZero();
      Scalar curr_energy = 0;
      if (m_externalForces[i]->isConservative())
      {
	  RodExternalConservativeForce* potentialForce = dynamic_cast<RodExternalConservativeForce*>(m_externalForces[i]);
	  potentialForce->computeForceEnergy(m_rod, curr_force, curr_energy);
	  f      += curr_force;
	  energy += curr_energy;
      }
      TraceStream(g_log, "") << m_externalForces[i]->getName() << " &rod = " << &m_rod << " potential energy = " << energy << " force norm = " << curr_force.norm() << '\n';
    }
  }


  /**
   * Evaluates the Jacobian of the forces on the rod.
   *
   * \param[out] J The Jacobian of the forces on the rod.
   */
  void evaluatePDotDX(Scalar scale, MatrixBase& J)
  {
    m_rod.computeJacobian(0, scale, J);

//     if (m_rod.viscous()) {
//       J.finalize();
//       J.scale(1.0 / m_diffEqSolver->getTimeStep());
//     }

    for (size_t i = 0; i < m_externalForces.size(); ++i) {
      m_externalForces[i]->computeForceDX(0, m_rod, scale, J);
    }
  }

  void evaluatePDotDV(Scalar scale, MatrixBase& J)
  {
    for (size_t i = 0; i < m_externalForces.size(); ++i) {
      m_externalForces[i]->computeForceDV(0, m_rod, scale, J);
    }
  }

  /**
   * This function returns the mass associated with a degree of
   * freedom of the rod.
   *
   * \param[in] i Which degree of freedom.
   * \return The mass of the degree of freedom
   */
  Scalar getMass(int i)
  {
    assert(i >= 0);
    assert(i < m_rod.ndof());

    return m_rod.getMass(i);
  }

  Scalar getX(int i) const
  {
    return m_rod.getDof(i);
  }

  void setX(int i, Scalar x)
  {
    m_rod.setDof(i, x);
  }

  Scalar getV(int i)
  {
    return m_rod.getVel(i);
  }

  void setV(int i, Scalar v)
  {
    m_rod.setVel(i, v);
  }

  void increment_q(const VecXd& dq)
  {
    for (int i = 0; i < dq.size(); ++i) {
      //setX(i, getX(i) + dq[i]);
      //setX(i, m_rod.getDof(i) + dq[i]);
      m_rod.setDof(i, m_rod.getDof(i) + dq[i]);
    }
  }
  
  void set_q(const VecXd& q)
  {
    assert( q.size() == ndof() );
    for (int i = 0; i < q.size(); ++i) 
    {
      m_rod.setDof(i,q(i));
    }
  }

  void increment_qdot(const VecXd& dqd)
  {
    // TODO: double check the indexing here :).
    if (m_rod.quasistatic()) {
      for (int i = 0; i < m_rod.nv(); ++i) {
        for (int coord = 0; coord < 3; ++coord) {
          int idx = m_rod.vertIdx(i, coord);
          //setV(idx, getV(idx) + dqd[idx]);
          //setV(idx, m_rod.getVel(idx) + dqd[idx]);
          m_rod.setVel(idx, m_rod.getVel(idx) + dqd[idx]);
        }
      }
    } else {
      for (int i = 0; i < dqd.size(); ++i) {
        //setV(i, getV(i) + dqd[i]);
        //setV(i, m_rod.getVel(i) + dqd[i]);
        m_rod.setVel(i, m_rod.getVel(i) + dqd[i]);
      }
    }
  }
  
  void set_qdot(const VecXd& qd)
  {
    assert( qd.size() == ndof() );
    if (m_rod.quasistatic()) 
    {
      for (int i = 0; i < m_rod.nv(); ++i) 
      {
        for (int coord = 0; coord < 3; ++coord) 
        {
          int idx = m_rod.vertIdx(i, coord);
          m_rod.setVel(idx, qd(idx));
        }
      }
    }
    else
    {
      for (int i = 0; i < qd.size(); ++i)
      {
        m_rod.setVel(i, qd(i));
      }
    }
  }

  /**
   * Adds an external force to be applied to the rod. On destruction,
   * this class will be responsible for de-allocating the memory
   * associated to the force.
   *
   * \param[in] force The external force to be applied to rods.
   */
  virtual void addExternalForce(RodExternalForce* force)
  {
    m_externalForces.push_back(force);
  }

  std::vector<RodExternalForce*>& getExternalForces()
  {
    return m_externalForces;
  }

  void startStep()
  {
    m_rod.viscousUpdate();
  }

  void endStep()
  {
    m_rod.updateProperties();
  }

  void startIteration() {}

  void endIteration()
  {
    m_rod.updateProperties();
  }

  int getMaxIterations() const
  {
    return m_diffEqSolver->getMaxIterations();
  }

  void setMaxIterations(int iterations)
  {
    m_diffEqSolver->setMaxIterations(iterations);
  }

  virtual Scalar get_stol() const { return m_diffEqSolver->get_stol(); }
  virtual void set_stol(Scalar s) { m_diffEqSolver->set_stol(s); }

  virtual Scalar get_atol() const { return m_diffEqSolver->get_atol(); }
  virtual void set_atol(Scalar a) { m_diffEqSolver->set_atol(a); }

  virtual Scalar get_rtol() const { return m_diffEqSolver->get_rtol(); }
  virtual void set_rtol(Scalar r) { m_diffEqSolver->set_rtol(r); }

  virtual Scalar get_inftol() const { return m_diffEqSolver->get_inftol(); }
  virtual void set_inftol(Scalar i) { m_diffEqSolver->set_inftol(i); }

  MatrixBase* createMatrix() const
  {
    SolverUtils* s = SolverUtils::instance();
    return s->createBandMatrix(m_rod.ndof(), m_rod.ndof(), 10, 10);
  }

  RodBoundaryCondition* getBoundaryCondition()
  {
    //if (m_boundaryCondition == NULL) {
    //  m_boundaryCondition = new RodBoundaryCondition(m_rod);
    //}
    //return m_boundaryCondition;
    return m_rod.getBoundaryCondition();
  }

  // TODO: Does anyone use this method? It seems kind of dangerous letting any user
  //       mess with our internal pointers :)
  void setBoundaryCondition(RodBoundaryCondition* bc)
  {
    assert( false );
    //m_boundaryCondition = bc;
  }

  void getScriptedDofs(IntArray& indices, std::vector<Scalar>& desired)
  {
    //if (m_boundaryCondition == NULL) {
    //  indices.resize(0);
    //  desired.resize(0);
    //  return;
    //}

    const RodBoundaryCondition::BCList& verts
        = m_rod.getBoundaryCondition()->scriptedVertices();
      //= m_boundaryCondition->scriptedVertices();
    const RodBoundaryCondition::BCList& edges
        = m_rod.getBoundaryCondition()->scriptedEdges();
      //= m_boundaryCondition->scriptedEdges();

    int nb = 3 * (int)verts.size() + (int)edges.size(); // # of scripted dofs 
    indices.resize(nb);
    desired.resize(nb);

    for (size_t i = 0; i < verts.size(); ++i) {
      for (int k = 0; k < 3; ++k) {
        indices[3 * i + k] = m_rod.vertIdx(verts[i], k);
	// std::cout << "RodTimeStepper is calling RodBoundaryCondition at getTime() = " << getTime() << std::endl;
        desired[3 * i + k]
          = m_rod.getBoundaryCondition()->getDesiredVertexPosition(verts[i], getTime())[k];
          //= m_boundaryCondition->getDesiredVertexPosition(verts[i])[k];
      }
    }

    for (size_t i = 0; i < edges.size(); ++i) {
      indices[3 * verts.size() + i] = m_rod.edgeIdx(edges[i]);
      desired[3 * verts.size() + i]
        = m_rod.getBoundaryCondition()->getDesiredEdgeAngle(edges[i], getTime());
        //= m_boundaryCondition->getDesiredEdgeAngle(edges[i]);
    }
  }

  ElasticRod* getRod()
  {
    return &m_rod;
  }

  void backup()
  {
    m_backupstate.backupRod(m_rod);
  }
  
  void backupResize()
  {
    m_backupstate.resize(m_rod);
  }
  
  void backupRestore()
  {
    m_backupstate.restoreRod(m_rod);
  }
  
  void backupClear()
  {
    m_backupstate.clear();
  }

protected:

  ElasticRod& m_rod;
  std::vector<RodExternalForce*> m_externalForces;
  Method m_method;
  DiffEqSolver* m_diffEqSolver;
  //RodBoundaryCondition* m_boundaryCondition;
  
  MinimalRodStateBackup m_backupstate;

  bool m_has_solved;
};

} // namespace BASim

#endif // RODTIMESTEPPER_HH