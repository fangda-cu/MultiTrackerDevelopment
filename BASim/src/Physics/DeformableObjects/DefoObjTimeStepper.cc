#include "DefoObjTimeStepper.hh"
#include "Shells/ShellSurfaceTensionForce.hh"

namespace BASim {

   
MatrixBase* DefoObjTimeStepper::createMatrix() const
{
  SolverUtils* s = SolverUtils::instance();
  s->setMatrixType(SolverUtils::EIGEN_SPARSE_MATRIX);
  MatrixBase* newMat = s->createSparseMatrix(m_obj.ndof(), m_obj.ndof(), 30);
  return newMat;
}

void DefoObjTimeStepper::evaluateConservativeForcesEnergy(VecXd& f, Scalar& energy)
{
  m_obj.computeConservativeForcesEnergy(f, energy);
   
}

Scalar DefoObjTimeStepper::determineMaxDt(const VecXd & force)
{
    VertexProperty<Vec3d> vforce(&m_obj);
    
    for (VertexIterator vit = m_obj.vertices_begin(); vit != m_obj.vertices_end(); ++vit)
    {
        int dofbase = m_obj.getPositionDofBase(*vit);
        vforce[*vit] = force.segment<3>(dofbase);
    }
    
    Scalar global_max_dt = 1.0;
    for (EdgeIterator eit = m_obj.edges_begin(); eit != m_obj.edges_end(); ++eit)
    {
        VertexHandle v0 = m_obj.fromVertex(*eit);
        VertexHandle v1 = m_obj.toVertex(*eit);
        Vec3d edge = m_obj.getVertexPosition(v1) - m_obj.getVertexPosition(v0);
        Vec3d force = vforce[v1] - vforce[v0];
        Scalar approaching_velocity = -force.dot(edge) / edge.squaredNorm();
        
        Scalar max_dt = 0;
        if (approaching_velocity >= 0)
            max_dt = std::numeric_limits<Scalar>::infinity();
        else
            max_dt = 1 / approaching_velocity * 0.8;    // allow the edge to be shrinked by 80% in one time step at most
        
        if (max_dt < global_max_dt)
            global_max_dt = max_dt;
    }
    
    return global_max_dt;
}

}

