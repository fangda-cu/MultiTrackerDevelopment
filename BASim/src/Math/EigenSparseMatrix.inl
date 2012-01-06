/**
 * \file PetscMatrix.inl
 *
 * \author miklos@cs.columbia.edu
 * \date 09/08/2009
 */

inline EigenSparseMatrix::EigenSparseMatrix(int s)
  : MatrixBase(s, s), m_dynamic(s,s), m_static(s,s)
{
}

inline EigenSparseMatrix::EigenSparseMatrix(int r, int c, int nnz)
  : MatrixBase(r, c), m_dynamic(r,c), m_static(r,c)
{
  m_dynamic.reserve(nnz);
}

inline EigenSparseMatrix::EigenSparseMatrix(const EigenSparseMatrix& M)
  : MatrixBase(m_rows, m_cols), m_dynamic(M.m_dynamic), m_static(M.m_static)
{
}

inline EigenSparseMatrix::~EigenSparseMatrix()
{
}


inline Scalar EigenSparseMatrix::operator() (int r, int c) const
{
  return m_dynamic.coeff(r,c);
}

inline int EigenSparseMatrix::set(int r, int c, Scalar val)
{
  m_dynamic.coeffRef(r,c) = val;
  return 0;
}

inline int EigenSparseMatrix::add(int r, int c, Scalar val)
{
   m_dynamic.coeffRef(r,c) += val;
   return 0;
}

inline int EigenSparseMatrix::add(const IntArray& rowIdx, const IntArray& colIdx,
                            const MatXd& values)
{
  
  for(unsigned int i = 0; i < rowIdx.size(); ++i) {
    for(unsigned int j = 0; j < colIdx.size(); ++j) {
      m_dynamic.coeffRef(rowIdx[i], colIdx[j]) += values(i,j);
    }
  }

   return 0;
}

inline int EigenSparseMatrix::add(const IndexArray& rowIdx, const IndexArray& colIdx,
                            const MatXd& values)
{
  for(int i = 0; i < rowIdx.size(); ++i) {
    for(int j = 0; j < colIdx.size(); ++j) {
      m_dynamic.coeffRef(rowIdx[i], colIdx[i]) += values(i,j);
    }
  }

  return 0;
}

inline bool EigenSparseMatrix::isApproxSymmetric( Scalar eps ) const 
{
  //TODO Fill this in.
  return true;
}

inline std::string EigenSparseMatrix::name() const
{
   return "SimpleSparseMatrix";
}

