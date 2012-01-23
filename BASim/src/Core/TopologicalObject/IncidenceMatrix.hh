/**
* \file IncidenceMatrix.hh
*
* \author batty@cs.columbia.edu
* \date 04/07/2011 
*/

#ifndef INCIDENCEMATRIX_HH
#define INCIDENCEMATRIX_HH

#include <cassert>
#include <vector>
#include <fstream>

namespace BASim {

  //A simple std::vector-based sparse compressed row incidence matrix 
  //to store the topology of our simplex mesh structure.
  //It needs to be resize-able in order to add/delete simplices.
  class IncidenceMatrix {

    //We could also swap this for a more highly tuned dynamic sparse matrix structure,
    //if there's an appropriate one.

  private:

    //Matrix dimensions
    unsigned int n_rows, n_cols;

    //For each row, a list of all column indices (sorted). The sign indicates whether the value is intended to be: +1 or -1.
    //NOTE: We shift all column indices up by 1, so the zero'th column is enabled to have a sign!
    std::vector< std::vector<int> > m_indices; 

  public:
    IncidenceMatrix();
    IncidenceMatrix(unsigned int rows, unsigned int cols);

    //Matrix dimensions
    unsigned int getNumRows() const {return n_rows;}
    unsigned int getNumCols() const {return n_cols;}
    void addRows(unsigned int rows);
    void addCols(unsigned int cols);

    //Regular accessors
    void set(unsigned int i, unsigned int j, int new_val);
    int  get(unsigned int i, unsigned int j) const;
    bool exists(unsigned int i, unsigned int j) const;
    void zero(unsigned int i, unsigned int j);
    void zeroRow(unsigned int i);
    void zeroAll();

    //Constant-time access within the row, by row-index rather than column number
    unsigned int getNumEntriesInRow(unsigned int row) const;
    unsigned int getColByIndex(unsigned int i, unsigned int index_in_row) const;
    int getValueByIndex(unsigned int i, unsigned int index_in_row) const;
    
    unsigned int getIndexByCol(unsigned int i, unsigned int col) const;

    void setByIndex(unsigned int i, unsigned int index_in_row, unsigned int col, int value);
    
    void cycleRow(unsigned int i);
    
    //Debugging
    void printMatrix() const;
    
    //Serialization
    static void serialize(std::ofstream& of, const IncidenceMatrix& val);
    static void load(std::ifstream& ifs, IncidenceMatrix& val);
  };

} // namespace BASim


#endif
