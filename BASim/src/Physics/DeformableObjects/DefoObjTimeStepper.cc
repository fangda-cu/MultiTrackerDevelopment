#include "DefoObjTimeStepper.hh"


namespace BASim {

   
MatrixBase* DefoObjTimeStepper::createMatrix() const
{
  SolverUtils* s = SolverUtils::instance();
  s->setMatrixType(SolverUtils::PETSC_MATRIX);
  MatrixBase* newMat = s->createSparseMatrix(m_obj.ndof(), m_obj.ndof(), 30);
  return newMat;
}



}