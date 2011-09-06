#include "IncidenceMatrix.hh"

#include "BAsim/src/io/SerializationUtils.hh"

namespace BASim {

int signum(int val) {
  return (val >= 0 ? 1 : -1);

}

IncidenceMatrix::IncidenceMatrix() : 
   n_rows(0), n_cols(0), m_indices(0)
{
}

IncidenceMatrix::IncidenceMatrix(unsigned int rows, unsigned int cols) : 
   n_rows(rows), n_cols(cols), 
   m_indices(rows, std::vector<int>())
{
}

unsigned int IncidenceMatrix::getNumEntriesInRow(unsigned int row) const {
 assert(row < n_rows);
 return m_indices[row].size();
}

int IncidenceMatrix::getValueByIndex(unsigned int i, unsigned int index_in_row) const {
 assert(i < n_rows);
 assert(index_in_row < m_indices[i].size());

 return signum(m_indices[i][index_in_row]);
}

unsigned int IncidenceMatrix::getColByIndex(unsigned int i, unsigned int index_in_row) const {
  assert(i < n_rows);
  assert(index_in_row < m_indices[i].size());

  return abs(m_indices[i][index_in_row]) - 1;
}

void IncidenceMatrix::set(unsigned int i, unsigned int j, int new_val) {
   assert(i < n_rows && j < n_cols);
   if(new_val == 0) {
     zero(i,j);
     return;
   }

   assert(new_val == 1 || new_val == -1);
  
   int colShift = j+1;
   for(unsigned int k=0; k<m_indices[i].size(); ++k){
      if(abs(m_indices[i][k])==colShift){
         m_indices[i][k]=signum(new_val)*colShift; //change its sign to match
         return;
      }
      else if(abs(m_indices[i][k])>(int)j){
         insert_vec(m_indices[i], (int)k, signum(new_val)*colShift); //set sign to match
         return;
      }
   }
   m_indices[i].push_back(signum(new_val)*colShift);
}

int IncidenceMatrix::get(unsigned int i, unsigned int j) const {
   assert(i < n_rows && j < n_cols);
  
   int colShift = j+1;
   for(unsigned int k=0; k<m_indices[i].size(); ++k){
      if(abs(m_indices[i][k])==(int)colShift){
         return signum(m_indices[i][k]);
      }
   }
   return 0;
}

void IncidenceMatrix::zero(unsigned int i, unsigned int j) {
   assert(i<n_rows && j < n_cols);

   int colShift = j+1;
   for(unsigned int k=0; k<m_indices[i].size(); ++k){
      if(abs(m_indices[i][k])==colShift){
         remove_vec(m_indices[i], k);
         return;
      }
      else if(abs(m_indices[i][k]) > colShift)
         return;
   }
}

bool IncidenceMatrix::exists(unsigned int i, unsigned int j) const {
   assert(i<n_rows && j < n_cols);
  
   int colShift = j+1;
   for(unsigned int k=0; k<m_indices[i].size(); ++k){
      if(abs(m_indices[i][k])==colShift){
         return true;
      }
   }
   return false;
}

void IncidenceMatrix::addRows(unsigned int rows) {
   n_rows += rows;
   m_indices.resize(n_rows);
}

void IncidenceMatrix::addCols(unsigned int cols) {
   n_cols += cols;
}

void IncidenceMatrix::zeroRow( unsigned int i )
{
   m_indices[i].clear();
}

void IncidenceMatrix::printMatrix() const 
{
   printf("Dimensions (%d,%d):\n", n_rows, n_cols);
   for(unsigned int row = 0; row < m_indices.size(); ++row) {
      printf("%d: ", row);
      for(unsigned int i = 0; i < m_indices[row].size(); ++i)
        printf("%c%d ", m_indices[row][i] > 0?'+':'-', abs(m_indices[row][i])-1);
      printf("\n");
   }
}

void IncidenceMatrix::zeroAll()
{
  for(unsigned int i = 0; i < getNumRows(); ++i)
    zeroRow(i);
}

void IncidenceMatrix::serialize(std::ofstream& of, const IncidenceMatrix& val) {
   assert( of.is_open() );
   
   serializeVal(of,val.n_rows);
   serializeVal(of,val.n_cols);
   serializeVectorVectorInt(of, val.m_indices);
}

void IncidenceMatrix::load(std::ifstream& ifs, IncidenceMatrix& val) {
   assert( ifs.is_open() );

   int rows, cols;
   loadVal(ifs,rows);
   loadVal(ifs,cols);
   val.n_rows = rows;
   val.n_cols = cols;
   loadVectorVectorInt(ifs, val.m_indices);
   
}


}