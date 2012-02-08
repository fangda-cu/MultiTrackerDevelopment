/**
 * \file EigenSparseMatrix.inl
 *
 * \author batty@cs.columbia.edu
 * \date Feb 8, 2012
 */

inline EigenSparseMatrix::EigenSparseMatrix(int s)
  : MatrixBase(s, s), m_dynamic(s,s)
{
}

inline EigenSparseMatrix::EigenSparseMatrix(int r, int c, int nnz)
  : MatrixBase(r, c), m_dynamic(r,c)
{
  //nnz is est #per row, so multiply to get a good estimate.
  m_triplets.reserve(nnz);
  m_dynamic.reserve(nnz);
}

inline EigenSparseMatrix::EigenSparseMatrix(const EigenSparseMatrix& M)
  : MatrixBase(m_rows, m_cols), m_dynamic(M.m_dynamic), m_triplets(M.m_triplets)
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
  //m_dynamic.coeffRef(r,c) = val;
  //TODO: Not implemented
  return 0;
}

inline int EigenSparseMatrix::add(int r, int c, Scalar val)
{
   //m_dynamic.coeffRef(r,c) += val;
  m_triplets.push_back( Eigen::Triplet<Scalar>(r,c,val) );
   return 0;
}

inline int EigenSparseMatrix::add(const IntArray& rowIdx, const IntArray& colIdx,
                            const MatXd& values)
{
  
  for(unsigned int i = 0; i < rowIdx.size(); ++i) {
    for(unsigned int j = 0; j < colIdx.size(); ++j) {
      m_triplets.push_back( Eigen::Triplet<Scalar>(rowIdx[i],colIdx[j],values(i,j)) );
      //m_dynamic.coeffRef(rowIdx[i], colIdx[j]) += values(i,j);
    }
  }

   return 0;
}

inline int EigenSparseMatrix::add(const IndexArray& rowIdx, const IndexArray& colIdx,
                            const MatXd& values)
{
  for(int i = 0; i < rowIdx.size(); ++i) {
    for(int j = 0; j < colIdx.size(); ++j) {
      m_triplets.push_back( Eigen::Triplet<Scalar>(rowIdx[i],colIdx[j],values(i,j)) );
      //m_dynamic.coeffRef(rowIdx[i], colIdx[j]) += values(i,j);
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

