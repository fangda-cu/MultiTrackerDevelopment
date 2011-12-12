#include "PetscMatrix.hh"

namespace BASim {

int PetscMatrix::setZero()
{
  int ierr = MatZeroEntries(m_M);
  CHKERRQ(ierr);
  return 0;
}

int PetscMatrix::resetNonzeros() {
   //now allow non-zero pattern edits again
   int ierr = MatSetOption(m_M, MAT_NEW_NONZERO_LOCATIONS, PETSC_TRUE);
   CHKERRQ(ierr);
   ierr = MatSetOption(m_M, MAT_NEW_NONZERO_LOCATION_ERR, PETSC_FALSE);
   CHKERRQ(ierr);
   return 0;
}

int PetscMatrix::zeroRows(const IntArray& idx, Scalar diag)
{
  if(idx.size() < 1) return 0;
  MatSetOption(m_M, MAT_KEEP_NONZERO_PATTERN, PETSC_TRUE);
  int ierr = MatZeroRows(m_M, idx.size(), &idx[0], diag);
  CHKERRQ(ierr);
  return 0;
}
/*
  int PetscMatrix::assemble()
  {
  int ierr = MatAssemblyBegin(m_M, MAT_FINAL_ASSEMBLY);
  CHKERRQ(ierr);
  ierr = MatAssemblyEnd(m_M, MAT_FINAL_ASSEMBLY);
  CHKERRQ(ierr);
  return 0;
  }
*/

int PetscMatrix::finalize()
{
   int ierr = MatAssemblyBegin(m_M, MAT_FINAL_ASSEMBLY);
   CHKERRQ(ierr);
   ierr = MatAssemblyEnd(m_M, MAT_FINAL_ASSEMBLY);
   CHKERRQ(ierr);
   return 0;
}

int PetscMatrix::finalizeNonzeros() {
   int ierr = MatAssemblyBegin(m_M, MAT_FINAL_ASSEMBLY);
   CHKERRQ(ierr);
   ierr = MatAssemblyEnd(m_M, MAT_FINAL_ASSEMBLY);
   CHKERRQ(ierr);
   ierr = MatSetOption(m_M, MAT_NEW_NONZERO_LOCATIONS, PETSC_FALSE);
   CHKERRQ(ierr);
   ierr = MatSetOption(m_M, MAT_NEW_NONZERO_LOCATION_ERR, PETSC_TRUE);
   CHKERRQ(ierr);
   return 0;
}

}