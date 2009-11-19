/**
 * \file MKLMatrix.hh
 *
 * \author acoull@wetafx.co.nz
 * \date 11/10/2009
 */

#ifndef MKL_MATRIX_HH
#define MKL_MATRIX_HH

#include <vector>

using namespace std;

namespace BASim {

class MKLMatrix : public MatrixBase
{
public:
  MKLMatrix(int s);
  MKLMatrix(int r, int c, int nnz = 1);
  MKLMatrix(const MKLMatrix& M);
  ~MKLMatrix();

  Scalar operator() (int i, int j) const;
  int set(int r, int c, Scalar val);
  int add(int r, int c, Scalar val);
  int add(const IntArray& rowIdx, const IntArray& colIdx, const MatXd& values);
  int add(const IndexArray& rowIdx, const IndexArray& colIdx,
          const MatXd& values);
  int scale(Scalar val);
  int setZero();
  int zeroRows(const IntArray& idx, Scalar diag = 1.0);
  int multiply(VecXd& y, Scalar s, const VecXd& x);
    
  int finalize();
  void transpose();
  
  //const Mat& getMKLMatrix() const { return m_M; }
  //Mat& getMKLMatrix() { return m_M; }

  size_t numRows() const { return m_rows; }
  size_t numColumns() const { return m_cols; }
      
  vector<Scalar>& getData() { return m_data; }
  vector<Scalar> copyOfData() const { return m_data; }
  
protected:

  vector<Scalar> m_data;
  //Mat m_M;       // the underlying Petsc matrix
  //Vec m_v, m_w;  // m_v is size rows(), m_w is size cols()
};

#include "MKLMatrix.inl"

} // namespace BASim

  
  /*
    
    MKLMatrix();
    MKLMatrix( size_t i_size );
    MKLMatrix( size_t i_rows, size_t i_columns );
    MKLMatrix( const MKLMatrix& i_M );
    ~MKLMatrix();
   
    size_t numRows() const { return m_numRows; }
   size_t numColumns() const { return m_numColumns; }
   vector<double>& data() { return m_data; }
   vector<double> copyOfData() const { return m_data; }
   
   void resize( size_t i_numRows, size_t i_numColumns );

   int add( uint i_row, uint i_column, Real i_value );

   void transpose();

   // C = A * B
   static int multiply( MKLMatrix& i_A,  MKLMatrix& i_B, MKLMatrix& o_C);
   
   // b = r*M*a
   int multiply(Real r, const Vector<Real>& a, Vector<Real>& b);
   
   // C = P^T * A
   static int PtA( MKLMatrix& P,  MKLMatrix& A, MKLMatrix& C);

    void set(uint o_row, uint o_column, Real value);

   Real get( uint i_row, uint i_column ) const;

   int clear();
   
   void print();
    void printMatlabFile(std::string filename);
   
   void getSuperDiagonal( vector<double>& d );
   void getDiagonal( vector<double>& d );
   void getSubDiagonal( vector<double>& d );
   
   // =================== Not Implemented below this line ===================
   
   int diagonalScale(Vector<Real> &a);

   int shift(Real a);

   // b = P^T * a
   static int Pta(const MKLMatrix& P, const Vector<Real>& a, Vector<Real>& b);
   

   int zeroRows(std::vector<int>& idx, Real diag = 1.0);

   int finalize() const;

protected:
    size_t m_numRows;
    size_t m_numColumns;
    vector<double> m_data;
};


/*
class MKLMatrix 
{
public:
    MKLMatrix();
    MKLMatrix( size_t i_rows, size_t i_columns );
    ~MKLMatrix(); 
    
    size_t rows() const { return m_rows; }
    size_t columns() const { return m_columns; }
    vector<double>& data() { return m_data; }
    
    double& operator()( size_t row, size_t col );
    friend MKLMatrix operator*( const MKLMatrix& m1, const MKLMatrix &m2 );
    
    friend std::ostream& operator<< (std::ostream& os, MKLMatrix& M)
    {
        os << "[";
        for ( size_t i = 0; i < M.rows(); i++ )
            for ( size_t j = 0; j < M.columns(); j++)
            {
                os << M( i, j ) << ( j < M.columns()-1 ? " " : (i < M.rows()-1 ? "\n" : "]") );
            }
            
        return os;
    }
    
private:
    size_t m_rows;
    size_t m_columns;
    vector<double> m_data;
};*/

#endif

