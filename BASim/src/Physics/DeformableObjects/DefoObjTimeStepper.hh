/**
 * \file DefoObjTimeStepper.hh
 *
 * \author batty@cs.columbia.edu
 * \date April 18, 2011
 */

#ifndef DEFOOBJTIMESTEPPER_H
#define DEFOOBJTIMESTEPPER_H

#include "BASim/src/Core/ObjectControllerBase.hh"

#include "BASim/src/Math/TimeSteppingBase.hh"
#include "BASim/src/Math/SolverUtils.hh"
#include "BASim/src/Math/SymplecticEuler.hh"
#include "BASim/src/Math/ImplicitEuler.hh"
#include "BASim/src/Math/SymmetricImplicitEuler.hh"
//#include "BASim/src/Math/StaticsSolver.hh"
#include "BASim/src/Physics/DeformableObjects/DeformableObject.hh"

namespace BASim {


/** Class to time step a deformable object. */
class DefoObjTimeStepper : public ObjectControllerBase
{
public:

  enum Method { SYMPL_EULER, IMPL_EULER, STATICS, NONE };

  DefoObjTimeStepper(DeformableObject& obj)
    : m_obj(obj)
    , m_method(NONE)
    , m_diffEqSolver(NULL)
  {
    setDiffEqSolver(SYMPL_EULER);
  }

  ~DefoObjTimeStepper()
  {
    if (m_diffEqSolver != NULL) delete m_diffEqSolver;
    
   /* for( int i = 0; i < (int) m_externalForces.size(); ++i )
    {
      assert( m_externalForces[i] != NULL );
      delete m_externalForces[i];
      m_externalForces[i] = NULL;
    }	*/
  }

  bool execute()
  {
    assert( getTimeStep() == m_obj.getTimeStep() );
    setTime(getTime() + getTimeStep());
    return m_diffEqSolver->execute();
  }

  void setTime(Scalar time)
  {
    m_diffEqSolver->setTime(time);
  }

  Scalar getTime() const
  {
    return m_diffEqSolver->getTime();
  }

  virtual void setTimeStep(Scalar dt)
  {
    m_diffEqSolver->setTimeStep(dt);
    m_obj.setTimeStep(dt);
  }

  Scalar getTimeStep() const
  {
    return m_diffEqSolver->getTimeStep();
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
      m_diffEqSolver = new SymplecticEuler<DefoObjTimeStepper>(*this);
    } 
    else if (method == IMPL_EULER) 
    {
      //m_diffEqSolver = new ImplicitEuler<DefoObjTimeStepper>(*this);
      m_diffEqSolver = new SymmetricImplicitEuler<DefoObjTimeStepper>(*this);
      //m_diffEqSolver = new BacktrackingImplicitEuler<RodTimeStepper>(*this);
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
      m_obj.setTimeStep(m_diffEqSolver->getTimeStep());
    }
  }

  int ndof() const
  {
    return m_obj.ndof();
  }

  /**
   * Returns all "masses" in this differential equation in a flat vector.
   *
   * \param[out] masses Vector that will contain masses. Must be the proper size.
   */
  void getMass( VecXd& masses )
  {
    assert( masses.size() == m_obj.ndof() );
    for( int i = 0; i < m_obj.ndof(); ++i ) {
      masses(i) = m_obj.getMass(i);
    }
  }

  void setX( const VecXd& positions )
  {
    assert( positions.size() == m_obj.ndof() );
    for( int i = 0; i < m_obj.ndof(); ++i ) m_obj.setDof(i, positions(i));
  }

  /**
   * Returns all "positions" in this differential equation in a flat vector.
   *
   * \param[out] masses Vector that will contain positions. Must be the proper size.
   */
  void getX( VecXd& positions ) const
  {
    assert( positions.size() == m_obj.ndof() );
    for( int i = 0; i < m_obj.ndof(); ++i ) positions(i) = m_obj.getDof(i);
  }
  
  void setV( const VecXd& velocities )
  {
    assert( velocities.size() == m_obj.ndof() );
    for( int i = 0; i < m_obj.ndof(); ++i ) m_obj.setVel(i, velocities(i));
  }
  
  /**
   * Returns all "velocities" in this differential equation in a flat vector.
   *
   * \param[out] masses Vector that will contain velocities. Must be the proper size.
   */
  void getV( VecXd& velocities ) const
  {
    assert( velocities.size() == m_obj.ndof() );
    for( int i = 0; i < m_obj.ndof(); ++i ) velocities(i) = m_obj.getVel(i);
  }

  /**
   * This function computes the force on each degree of freedom
   * associated to the rod.
   *
   * \param[out] f The vector of accelerations on the rod.
   */
  void evaluatePDot(VecXd& f)
  {
    // add internal forces
    m_obj.computeForces(f);

    //if (m_rod.viscous()) f /= m_diffEqSolver->getTimeStep();

    // add external forces
    /*for (size_t i = 0; i < m_externalForces.size(); ++i) {
      m_externalForces[i]->computeForce(m_rod, f);
    }*/
  }

  /**
   * Evaluates the Jacobian of the forces on the rod.
   *
   * \param[out] J The Jacobian of the forces on the rod.
   */
  void evaluatePDotDX(Scalar scale, MatrixBase& J)
  {
    m_obj.computeJacobian(scale, J);

//     if (m_rod.viscous()) {
//       J.finalize();
//       J.scale(1.0 / m_diffEqSolver->getTimeStep());
//     }

    /*for (size_t i = 0; i < m_externalForces.size(); ++i) {
      m_externalForces[i]->computeForceDX(0, m_rod, scale, J);
    }*/
  }

  void evaluatePDotDV(Scalar scale, MatrixBase& J)
  {
   /* for (size_t i = 0; i < m_externalForces.size(); ++i) {
      m_externalForces[i]->computeForceDV(0, m_rod, scale, J);
    }*/
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
    assert(i < m_obj.ndof());

    return m_obj.getMass(i);
  }

  Scalar getX(int i) const
  {
    return m_obj.getDof(i);
  }

  void setX(int i, Scalar x)
  {
    m_obj.setDof(i, x);
  }

  Scalar getV(int i)
  {
    return m_obj.getVel(i);
  }

  void setV(int i, Scalar v)
  {
    m_obj.setVel(i, v);
  }

  void increment_q(const VecXd& dq)
  {
    for (int i = 0; i < dq.size(); ++i) {
      //setX(i, getX(i) + dq[i]);
      //setX(i, m_rod.getDof(i) + dq[i]);
      m_obj.setDof(i, m_obj.getDof(i) + dq[i]);
    }
  }
  
  void set_q(const VecXd& q)
  {
    assert( q.size() == ndof() );
    for (int i = 0; i < q.size(); ++i) 
    {
      m_obj.setDof(i,q(i));
    }
  }

  void increment_qdot(const VecXd& dqd)
  {
    // TODO: double check the indexing here :).
    //if (m_obj.quasistatic()) {
    //  for (int i = 0; i < m_obj.nv(); ++i) {
    //    for (int coord = 0; coord < 3; ++coord) {
    //      int idx = m_m_objrod.vertIdx(i, coord);
    //      //setV(idx, getV(idx) + dqd[idx]);
    //      //setV(idx, m_rod.getVel(idx) + dqd[idx]);
    //      m_rod.setVel(idx, m_rod.getVel(idx) + dqd[idx]);
    //    }
    //  }
    //} else {
      for (int i = 0; i < dqd.size(); ++i) {
        //setV(i, getV(i) + dqd[i]);
        //setV(i, m_rod.getVel(i) + dqd[i]);
        m_obj.setVel(i, m_obj.getVel(i) + dqd[i]);
      }
    //}
  }
  
  void set_qdot(const VecXd& qd)
  {
    assert( qd.size() == ndof() );
   /* if (m_obj.quasistatic()) 
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
    {*/
      for (int i = 0; i < qd.size(); ++i)
      {
        m_obj.setVel(i, qd(i));
      }
    //}
  }

  /**
   * Adds an external force to be applied to the rod. On destruction,
   * this class will be responsible for de-allocating the memory
   * associated to the force.
   *
   * \param[in] force The external force to be applied to rods.
   */
 /* void addExternalForce(RodExternalForce* force)
  {
    m_externalForces.push_back(force);
  }*/

 /* std::vector<RodExternalForce*>& getExternalForces()
  {
    return m_externalForces;
  }*/

  void startStep()
  {
    //m_rod.viscousUpdate();
    m_obj.setTime(m_diffEqSolver->getTime());
    m_obj.startStep();
  }

  void endStep()
  {
    m_obj.endStep();
    //m_obj.updateProperties();
  }

  void startIteration() {}

  void endIteration()
  {
    //m_obj.updateProperties();
  }

  int getMaxIterations() const
  {
    return m_diffEqSolver->getMaxIterations();
  }

  void setMaxIterations(int iterations)
  {
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


  void getScriptedDofs(IntArray& indices, std::vector<Scalar>& desired)
  {
    Scalar time = getTime();
    m_obj.getScriptedDofs(indices, desired, time);
   
  }

  DeformableObject* getObj()
  {
    return &m_obj;
  }

  void backup()
  {
    //m_backupstate.backupRod(m_rod);
  }
  
  void backupResize()
  {
    //m_backupstate.resize(m_rod);
  }
  
  void backupRestore()
  {
    //m_backupstate.restoreRod(m_rod);
  }
  
  void backupClear()
  {
    //m_backupstate.clear();
  }

protected:

  DeformableObject& m_obj;
  
  Method m_method;
  DiffEqSolver* m_diffEqSolver;
  
};

} // namespace BASim

#endif // RODTIMESTEPPER_HH
