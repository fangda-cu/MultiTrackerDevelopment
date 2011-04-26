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
#include "../../Math/SymplecticEuler.hh"
#include "../../Math/SymmetricImplicitEuler.hh"
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
  enum Method { SYMPL_EULER, IMPL_EULER, NONE };

  MultipleRodTimeStepper()
  : m_method(NONE)
  , m_rods()
  , m_externalForces()
  , m_diffEqSolver(NULL)
  , m_rodRodexternalForces()
  //, m_boundaryCondition(NULL)
  {
    setDiffEqSolver(SYMPL_EULER);
  }

  ~MultipleRodTimeStepper()
  {
    if( m_diffEqSolver != NULL )
    {
      delete m_diffEqSolver;
      m_diffEqSolver = NULL;
    }

    for( int i = 0; i < (int) m_externalForces.size(); ++i )
    {
      assert( m_externalForces[i] != NULL );
      delete m_externalForces[i];
      m_externalForces[i] = NULL;
    }
    
    for( int i = 0; i < (int) m_rodRodexternalForces.size(); ++i )
    {
      assert( m_rodRodexternalForces[i] != NULL );
      delete m_rodRodexternalForces[i];
      m_rodRodexternalForces[i] = NULL;
    }
  }

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

  void setTimeStep( Scalar dt )
  {
    assert( m_diffEqSolver != NULL );
    m_diffEqSolver->setTimeStep(dt);
    for( int i = 0; i < (int) m_rods.size(); ++i )
    {
      assert( m_rods[i] != NULL );
      m_rods[i]->setTimeStep(dt);
    }
  }

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

  void setDiffEqSolver( Method method )
  {
    if( method == m_method ) return;

    m_method = method;
    if( m_diffEqSolver != NULL ) delete m_diffEqSolver;
    m_diffEqSolver = NULL;

    switch( method ) 
    {
      case SYMPL_EULER:
      {
        assert( m_diffEqSolver == NULL );
        m_diffEqSolver = new SymplecticEuler<MultipleRodTimeStepper>(*this);
        break;
      }
      case IMPL_EULER:
      {
        assert( m_diffEqSolver == NULL );
        m_diffEqSolver = new SymmetricImplicitEuler<MultipleRodTimeStepper>(*this);
        //m_diffEqSolver = new ImplicitEuler<MultipleRodTimeStepper>(*this);
        break;
      }
      //case STATICS:
      //{
      //  assert( m_diffEqSolver == NULL );
      //  m_diffEqSolver = new StaticsSolver<MultipleRodTimeStepper>(*this);
      //  break;
      //}
      default:
      {
        std::cerr << "\033[31;1mERROR IN MULTIPLERODTIMESTEPPER:\033[m Invalid integrator specified. " << std::endl;
        exit(0);
        break;
      }
    }

    assert( m_diffEqSolver != NULL );
    for( int i = 0; i < (int) m_rods.size(); ++i )
    {
      assert( m_rods[i] != NULL );
      m_rods[i]->setTimeStep(m_diffEqSolver->getTimeStep());
    }
  }

  int ndof() const
  {
    int numdof = 0;
    for( int i = 0; i < (int) m_rods.size(); ++i )
    {
      assert( m_rods[i] != NULL );
      numdof += m_rods[i]->ndof();
    }
    return numdof;
  }

  /**
   * Returns all "masses" in this differential equation in a flat vector.
   *
   * \param[out] masses Vector that will contain masses. Must be the proper size.
   */
  void getMass( VecXd& masses )
  {
    assert( masses.size() == ndof() );
    int current_dof = 0;
    for( int i = 0; i < (int) m_rods.size(); ++i )
    {
      assert( m_rods[i] != NULL );
      for( int j = 0; j < m_rods[i]->ndof(); ++j )
      {
        masses(current_dof) = m_rods[i]->getMass(j);
        ++current_dof;
      }
    }
    assert( current_dof == masses.size() );
  }
  
  void setX( const VecXd& positions )
  {
    assert( positions.size() == ndof() );
    int current_dof = 0;
    for( int i = 0; i < (int) m_rods.size(); ++i )
    {
      assert( m_rods[i] != NULL );
      for( int j = 0; j < m_rods[i]->ndof(); ++j )
      {
        m_rods[i]->setDof(j, positions(current_dof));
        ++current_dof;
      }
    }
    assert( current_dof == positions.size() );
  }
  
  /**
   * Returns all "positions" in this differential equation in a flat vector.
   *
   * \param[out] masses Vector that will contain positions. Must be the proper size.
   */
  void getX( VecXd& positions ) const
  {
    assert( positions.size() == ndof() );
    int current_dof = 0;
    for( int i = 0; i < (int) m_rods.size(); ++i )
    {
      assert( m_rods[i] != NULL );
      for( int j = 0; j < m_rods[i]->ndof(); ++j )
      {
        positions(current_dof) = m_rods[i]->getDof(j);
        ++current_dof;
      }
    }
    assert( current_dof == positions.size() );
  }
  
  void setV( const VecXd& velocities )
  {
    assert( velocities.size() == ndof() );
    int current_dof = 0;
    for( int i = 0; i < (int) m_rods.size(); ++i )
    {
      assert( m_rods[i] != NULL );
      for( int j = 0; j < m_rods[i]->ndof(); ++j )
      {
        m_rods[i]->setVel(j, velocities(current_dof));
        ++current_dof;
      }
    }
    assert( current_dof == velocities.size() );
  }
  
  /**
   * Returns all "velocities" in this differential equation in a flat vector.
   *
   * \param[out] masses Vector that will contain velocities. Must be the proper size.
   */
  void getV( VecXd& velocities ) const
  {
    assert( velocities.size() == ndof() );
    int current_dof = 0;
    for( int i = 0; i < (int) m_rods.size(); ++i )
    {
      assert( m_rods[i] != NULL );
      for( int j = 0; j < m_rods[i]->ndof(); ++j )
      {
        velocities(current_dof) = m_rods[i]->getVel(j);
        ++current_dof;
      }
    }
    assert( current_dof == velocities.size() );
  }  
  
  /**
   * This function computes the force on each degree of freedom
   * associated to the rod.
   *
   * \param[out] f The vector of accelerations on the rod.
   */
  void evaluatePDot( VecXd& f )
  {
    // Add internal forces.
    int baseindex = 0;
    for( int i = 0; i < (int) m_rods.size(); ++i )
    {
      VecXd local_forces(m_rods[i]->ndof());
      local_forces.setZero();
      m_rods[i]->computeForces(local_forces);
      f.segment(baseindex,m_rods[i]->ndof()) += local_forces;
      baseindex += m_rods[i]->ndof();
    }
    assert( baseindex == ndof() );

    // Add external forces.
    for( int i = 0; i < (int) m_externalForces.size(); ++i )
    {
      baseindex = 0;
      for( int j = 0; j < (int) m_rods.size(); ++j )
      {
        VecXd local_forces(m_rods[j]->ndof());
        local_forces.setZero();
        m_externalForces[i]->computeForce(*m_rods[j], local_forces);
        f.segment(baseindex,m_rods[j]->ndof()) += local_forces;
        baseindex += m_rods[j]->ndof();
      }
      assert( baseindex == ndof() );
    }
    
    for( int i = 0; i < (int) m_rodRodexternalForces.size(); ++i )
    {
      m_rodRodexternalForces[i]->computeForce(f);
    }
  }

  /**
   * Evaluates the Jacobian of the forces on the rod.
   *
   * \param[out] J The Jacobian of the forces on the rod.
   */
  void evaluatePDotDX( Scalar scale, MatrixBase& J )
  {
    assert( J.rows() == ndof() );
    assert( J.rows() == J.cols() );
    // Add internal forces.
    int baseindex = 0;
    for( int i = 0; i < (int) m_rods.size(); ++i )
    {
      assert( baseindex+m_rods[i]->ndof()-1 < J.rows() );
      m_rods[i]->computeJacobian(baseindex, scale, J);
      baseindex += m_rods[i]->ndof();
    }
    assert( baseindex == ndof() );

    // Add external forces.
    for( int i = 0; i < (int) m_externalForces.size(); ++i )
    {
      baseindex = 0;
      for( int j = 0; j < (int) m_rods.size(); ++j )
      {
        assert( baseindex+m_rods[j]->ndof()-1 < J.rows() );
        m_externalForces[i]->computeForceDX(baseindex, *m_rods[j], scale, J);
        baseindex += m_rods[j]->ndof();
      }
      assert( baseindex == ndof() );
    }
    
    for( int i = 0; i < (int) m_rodRodexternalForces.size(); ++i )
    {
      m_rodRodexternalForces[i]->computeForceDX(scale,J);
    }    
  }

  void evaluatePDotDV( Scalar scale, MatrixBase& J )
  {
    assert( J.rows() == ndof() );
    assert( J.rows() == J.cols() );
    // Add external forces.
    for( int i = 0; i < (int) m_externalForces.size(); ++i )
    {
      int baseindex = 0;
      for( int j = 0; j < (int) m_rods.size(); ++j )
      {
        assert( baseindex+m_rods[j]->ndof()-1 < J.rows() );
        m_externalForces[i]->computeForceDV(baseindex, *m_rods[j], scale, J);
        baseindex += m_rods[j]->ndof();
      }
      assert( baseindex == ndof() );
    }
    
    for( int i = 0; i < (int) m_rodRodexternalForces.size(); ++i )
    {
      m_rodRodexternalForces[i]->computeForceDV(scale,J);
    }
  }

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

  void increment_q( const VecXd& dq )
  {
    assert( dq.size() == ndof() );
    int current_dof = 0;
    for( int i = 0; i < (int) m_rods.size(); ++i )
    {
      assert( m_rods[i] != NULL );
      for( int j = 0; j < m_rods[i]->ndof(); ++j )
      {
        m_rods[i]->setDof( j, m_rods[i]->getDof(j) + dq(current_dof) );
        ++current_dof;
      }
    }
    assert( current_dof == dq.size() );
  }

  void set_q(const VecXd& q)
  {
    assert( q.size() == ndof() );
    int current_dof = 0;
    for( int i = 0; i < (int) m_rods.size(); ++i )
    {
      assert( m_rods[i] != NULL );
      for( int j = 0; j < m_rods[i]->ndof(); ++j )
      {
        m_rods[i]->setDof( j, q(current_dof) );
        ++current_dof;
      }
    }
    assert( current_dof == q.size() );
  }

  // TODO: Double check the indexing here. 
  // TODO: Add some sanity checks?
  void increment_qdot( const VecXd& dqd )
  {
    assert( dqd.size() == ndof() );
    int base_idx = 0;
    for( int i = 0; i < (int) m_rods.size(); ++i )
    {
      if (m_rods[i]->quasistatic()) 
      {
        for (int j = 0; j < m_rods[i]->nv(); ++j) 
        {
          for (int coord = 0; coord < 3; ++coord) 
          {
            int idx = m_rods[i]->vertIdx(j, coord);
            //setV(idx, getV(idx) + dqd[idx]);
            //m_rod.getVel(i)
            m_rods[i]->setVel(idx, m_rods[i]->getVel(idx) + dqd[base_idx+idx]);
          }
        }
      } 
      else
      {
        for (int j = 0; j < m_rods[i]->ndof(); ++j)
        {
          //setV(i, getV(i) + dqd[i]);
          m_rods[i]->setVel(j, m_rods[i]->getVel(j) + dqd[base_idx+j]);
        }
      }
      base_idx += m_rods[i]->ndof();
    }
    assert( base_idx == dqd.size() );
  }

  void set_qdot(const VecXd& qd)
  {
    assert( qd.size() == ndof() );
    int base_idx = 0;
    for( int i = 0; i < (int) m_rods.size(); ++i )
    {
      if (m_rods[i]->quasistatic()) 
      {
        for (int j = 0; j < m_rods[i]->nv(); ++j) 
        {
          for (int coord = 0; coord < 3; ++coord) 
          {
            int idx = m_rods[i]->vertIdx(j, coord);
            m_rods[i]->setVel(idx, qd[base_idx+idx]);
          }
        }
      } 
      else
      {
        for (int j = 0; j < m_rods[i]->ndof(); ++j)
        {
          m_rods[i]->setVel(j, qd[base_idx+j]);
        }
      }
      base_idx += m_rods[i]->ndof();
    }
    assert( base_idx == qd.size() );

  
//    assert( qd.size() == ndof() );
//    if (m_rod.quasistatic()) 
//    {
//      for (int i = 0; i < m_rod.nv(); ++i) 
//      {
//        for (int coord = 0; coord < 3; ++coord) 
//        {
//          int idx = m_rod.vertIdx(i, coord);
//          m_rod.setVel(idx, qd(idx));
//        }
//      }
//    }
//    else
//    {
//      for (int i = 0; i < qd.size(); ++i)
//      {
//        m_rod.setVel(i, qd(i));
//      }
//    }
  }


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

  void startStep()
  {
    for( int i = 0; i < (int) m_rods.size(); ++i ) 
    {
      assert( m_rods[i] != NULL );
      m_rods[i]->viscousUpdate();
    }
  }

  void endStep()
  {
    for( int i = 0; i < (int) m_rods.size(); ++i ) 
    {
      assert( m_rods[i] != NULL );
      m_rods[i]->updateProperties();
    }
  }

  void startIteration() {}

  void endIteration()
  {
    for( int i = 0; i < (int) m_rods.size(); ++i )
    {
      assert( m_rods[i] != NULL );
      m_rods[i]->updateProperties();
    }
  }

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

  MatrixBase* createMatrix() const
  {
    SolverUtils* s = SolverUtils::instance();
    int parameter_to_get_rid_of = -1;
    int numdof = ndof();
    return s->createSparseMatrix(numdof, numdof, parameter_to_get_rid_of);
  }

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

  void getScriptedDofs( IntArray& indices, std::vector<Scalar>& desired )
  {
    // Compute the total number of scripted degrees of freedom
    int nb = 0;
    for( int i = 0; i < (int) m_rods.size(); ++i )
    {
      nb += 3 * (int)( m_rods[i]->getBoundaryCondition()->scriptedVertices().size() );
      nb +=     (int)( m_rods[i]->getBoundaryCondition()->scriptedEdges().size() );
    }

    indices.resize(nb);
    desired.resize(nb);
    
    int base_idx = 0;
    int base_nb = 0;
    for( int i = 0; i < (int) m_rods.size(); ++i )
    {
      const RodBoundaryCondition::BCList& verts
        = m_rods[i]->getBoundaryCondition()->scriptedVertices();
      const RodBoundaryCondition::BCList& edges
        = m_rods[i]->getBoundaryCondition()->scriptedEdges();
      
      for (int j = 0; j < (int) verts.size(); ++j) 
      {
        for (int k = 0; k < 3; ++k) 
        {
          indices[base_nb] = base_idx + m_rods[i]->vertIdx(verts[j], k);
	  std::cout << "MultipleRodTimeStepper is calling RodBoundaryCondition" << std::endl;
          desired[base_nb] = m_rods[i]->getBoundaryCondition()->getDesiredVertexPosition(verts[j], getTime())[k];
          ++base_nb;
        }
      }
      
      for (int j = 0; j < (int) edges.size(); ++j) 
      {
        indices[base_nb] = base_idx + m_rods[i]->edgeIdx(edges[j]);
        desired[base_nb] = m_rods[i]->getBoundaryCondition()->getDesiredEdgeAngle(edges[j], getTime());
        ++base_nb;
      }
      
      base_idx += m_rods[i]->ndof();
    }
    assert( base_nb == (int) indices.size() );
    assert( base_idx == ndof() );
  }

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
