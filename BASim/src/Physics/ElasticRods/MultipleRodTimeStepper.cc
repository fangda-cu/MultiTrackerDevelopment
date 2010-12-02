/**
 * \file MultipleRodTimeStepper.cc
 *
 * \author smith@cs.columbia.edu
 * \date 05/29/2010
 */

#include "MultipleRodTimeStepper.hh"
#include "../../Math/SymplecticEuler.hh"
#include "../../Math/ImplicitEuler.hh"
#include "../../Math/SymmetricImplicitEuler.hh"
#include "../../Math/StaticsSolver.hh"

namespace BASim
{

MultipleRodTimeStepper::MultipleRodTimeStepper()
  : m_method(NONE)
  , m_rods()
  , m_externalForces()
  , m_diffEqSolver(NULL)
  , m_rodRodexternalForces()
  //, m_boundaryCondition(NULL)
{
  setDiffEqSolver(SYMPL_EULER);
}

MultipleRodTimeStepper::~MultipleRodTimeStepper()
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

void MultipleRodTimeStepper::setTimeStep( Scalar dt )
{
  assert( m_diffEqSolver != NULL );
  m_diffEqSolver->setTimeStep(dt);
  for( int i = 0; i < (int) m_rods.size(); ++i )
  {
    assert( m_rods[i] != NULL );
    m_rods[i]->setTimeStep(dt);
  }
}

void MultipleRodTimeStepper::setDiffEqSolver( Method method )
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
    case SYM_IMPL_EULER:
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

int MultipleRodTimeStepper::ndof() const
{
  int numdof = 0;
  for( int i = 0; i < (int) m_rods.size(); ++i )
  {
    assert( m_rods[i] != NULL );
    numdof += m_rods[i]->ndof();
  }
  return numdof;
}

void MultipleRodTimeStepper::getMass( VecXd& masses )
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
  
void MultipleRodTimeStepper::setX( const VecXd& positions )
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
void MultipleRodTimeStepper::getX( VecXd& positions ) const
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
  
void MultipleRodTimeStepper::setV( const VecXd& velocities )
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
void MultipleRodTimeStepper::getV( VecXd& velocities ) const
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
void MultipleRodTimeStepper::evaluatePDot( VecXd& f )
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
void MultipleRodTimeStepper::evaluatePDotDX( Scalar scale, MatrixBase& J )
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

void MultipleRodTimeStepper::evaluatePDotDV( Scalar scale, MatrixBase& J )
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

void MultipleRodTimeStepper::increment_q( const VecXd& dq )
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

void MultipleRodTimeStepper::set_q(const VecXd& q)
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
void MultipleRodTimeStepper::increment_qdot( const VecXd& dqd )
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

void MultipleRodTimeStepper::set_qdot(const VecXd& qd)
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

void MultipleRodTimeStepper::startStep()
{
  for( int i = 0; i < (int) m_rods.size(); ++i ) 
  {
    assert( m_rods[i] != NULL );
    m_rods[i]->viscousUpdate();
  }
}

void MultipleRodTimeStepper::endStep()
{
  for( int i = 0; i < (int) m_rods.size(); ++i ) 
  {
    assert( m_rods[i] != NULL );
    m_rods[i]->updateProperties();
  }
}

void MultipleRodTimeStepper::endIteration()
{
  for( int i = 0; i < (int) m_rods.size(); ++i )
  {
    assert( m_rods[i] != NULL );
    m_rods[i]->updateProperties();
  }
}

MatrixBase* MultipleRodTimeStepper::createMatrix() const
{
  SolverUtils* s = SolverUtils::instance();
  int parameter_to_get_rid_of = -1;
  int numdof = ndof();
  return s->createSparseMatrix(numdof, numdof, parameter_to_get_rid_of);
}

void MultipleRodTimeStepper::getScriptedDofs( IntArray& indices, std::vector<Scalar>& desired )
{
  // Compute the total number of scripted degrees of freedom
  int nb = 0;
  for( int i = 0; i < (int) m_rods.size(); ++i )
  {
    nb += 3*(int)m_rods[i]->getBoundaryCondition()->scriptedVertices().size();
    nb +=   (int)m_rods[i]->getBoundaryCondition()->scriptedEdges().size();
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
        desired[base_nb] = m_rods[i]->getBoundaryCondition()->getDesiredVertexPosition(verts[j])[k];
        ++base_nb;
      }
    }
    
    for (int j = 0; j < (int) edges.size(); ++j) 
    {
      indices[base_nb] = base_idx + m_rods[i]->edgeIdx(edges[j]);
      desired[base_nb] = m_rods[i]->getBoundaryCondition()->getDesiredEdgeAngle(edges[j]);
      ++base_nb;
    }
    
    base_idx += m_rods[i]->ndof();
  }
  assert( base_nb == (int) indices.size() );
  assert( base_idx == ndof() );
}

} // namespace BASim
