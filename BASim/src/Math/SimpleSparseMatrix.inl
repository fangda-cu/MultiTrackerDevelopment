/**
 * \file PetscMatrix.inl
 *
 * \author miklos@cs.columbia.edu
 * \date 09/08/2009
 */

inline SimpleSparseMatrix::SimpleSparseMatrix(int s)
  : MatrixBase(s, s)

  , m_data(s, std::map<int,Scalar>())
{
}

inline SimpleSparseMatrix::SimpleSparseMatrix(int r, int c, int nnz)
  : MatrixBase(r, c), m_data(r, std::map<int,Scalar>())
{
}

inline SimpleSparseMatrix::SimpleSparseMatrix(const SimpleSparseMatrix& M)
  : MatrixBase(0, 0), m_data(M.m_data)
{
}

inline SimpleSparseMatrix::~SimpleSparseMatrix()
{
}


inline Scalar SimpleSparseMatrix::operator() (int r, int c) const
{
   std::map<int,Scalar>::const_iterator it = m_data[r].find(c);
   Scalar val = 0;
   if(it != m_data[r].end())
      val = it->second;
  return (Scalar) val;
}

inline int SimpleSparseMatrix::set(int r, int c, Scalar val)
{
  m_data[r][c] = val;
  return 0;
}

inline int SimpleSparseMatrix::add(int r, int c, Scalar val)
{
   m_data[r][c] += val;
   return 0;
}

inline int SimpleSparseMatrix::add(const IntArray& rowIdx, const IntArray& colIdx,
                            const MatXd& values)
{
  
  for(unsigned int i = 0; i < rowIdx.size(); ++i)
    for(unsigned int j = 0; j < colIdx.size(); ++j)
      m_data[rowIdx[i]][colIdx[j]] += values(i,j);

   return 0;
}

inline int SimpleSparseMatrix::add(const IndexArray& rowIdx, const IndexArray& colIdx,
                            const MatXd& values)
{
  for(int i = 0; i < rowIdx.size(); ++i)
    for(int j = 0; j < colIdx.size(); ++j)
      m_data[rowIdx[i]][colIdx[j]] += values(i,j);

  return 0;
}

inline bool SimpleSparseMatrix::isApproxSymmetric( Scalar eps ) const 
{
  //TODO Fill this in.
  return true;
}

inline std::string SimpleSparseMatrix::name() const
{
   return "SimpleSparseMatrix";
}

