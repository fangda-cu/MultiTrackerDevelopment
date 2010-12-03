/**
 * \file RodTimeStepper.cc
 *
 * \author miklos@cs.columbia.edu
 * \date 09/03/2009
 */

#include "RodTimeStepper.hh"
#include "../../Math/SymplecticEuler.hh"
#include "../../Math/ImplicitEuler.hh"
#include "../../Math/SymmetricImplicitEuler.hh"
#include "../../Math/StaticsSolver.hh"

namespace BASim {

RodTimeStepper::RodTimeStepper(ElasticRod& rod)
  : m_rod(rod)
  , m_method(NONE)
  , m_diffEqSolver(NULL)
  //, m_boundaryCondition(NULL)
  , m_backupstate()
{
  setDiffEqSolver(SYMPL_EULER);
}

RodTimeStepper::~RodTimeStepper()
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

void RodTimeStepper::setDiffEqSolver(Method method)
{
  if (method == m_method) return;

  m_method = method;
  if (m_diffEqSolver != NULL) delete m_diffEqSolver;
  m_diffEqSolver = NULL;

  if (method == SYMPL_EULER) {
    m_diffEqSolver = new SymplecticEuler<RodTimeStepper>(*this);

  } else if (method == IMPL_EULER) {
    m_diffEqSolver = new ImplicitEuler<RodTimeStepper>(*this);

  } else if (method == SYM_IMPL_EULER) {
    m_diffEqSolver = new SymmetricImplicitEuler<RodTimeStepper>(*this);

  } else if (method == STATICS) {
    m_diffEqSolver = new StaticsSolver<RodTimeStepper>(*this);

  } else if (method == NONE) {
    m_diffEqSolver = NULL;

  } else {
    std::cout << "Unknown method specified" << std::endl;
    m_diffEqSolver = NULL;

  }

  if (m_diffEqSolver != NULL) {
    m_rod.setTimeStep(m_diffEqSolver->getTimeStep());
  }
}

void RodTimeStepper::getScriptedDofs(IntArray& indices, std::vector<Scalar>& desired)
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
      desired[3 * i + k]
        = m_rod.getBoundaryCondition()->getDesiredVertexPosition(verts[i])[k];
        //= m_boundaryCondition->getDesiredVertexPosition(verts[i])[k];
    }
  }

  for (size_t i = 0; i < edges.size(); ++i) {
    indices[3 * verts.size() + i] = m_rod.edgeIdx(edges[i]);
    desired[3 * verts.size() + i]
      = m_rod.getBoundaryCondition()->getDesiredEdgeAngle(edges[i]);
      //= m_boundaryCondition->getDesiredEdgeAngle(edges[i]);
  }
}

} // namespace BASim
