// ---------------------------------------------------------
//
//  accelerationgrid.cpp
//  Tyson Brochu 2008
//  
//  A grid-based collision test culling structure.
//
// ---------------------------------------------------------

// ---------------------------------------------------------
// Includes
// ---------------------------------------------------------

#include <accelerationgrid.h>

#include <array3.h>
#include <limits>
#include <util.h>
#include <vec.h>
#include <vector>
#include <wallclocktime.h>

// ---------------------------------------------------------
// Global externs
// ---------------------------------------------------------

// ---------------------------------------------------------
// Local constants, typedefs, macros
// ---------------------------------------------------------

// ---------------------------------------------------------
// Static function definitions
// ---------------------------------------------------------

// ---------------------------------------------------------
// Member function definitions
// ---------------------------------------------------------

// --------------------------------------------------------
///
/// Default constructor
///
// --------------------------------------------------------

namespace ElTopo {

AccelerationGrid::AccelerationGrid() :
m_cells(0,0,0),
m_elementidxs(0),
m_elementxmins(0),
m_elementxmaxs(0),
m_elementquery(0),
m_lastquery(0),
m_gridxmin(0,0,0),
m_gridxmax(0,0,0),
m_cellsize(0,0,0),
m_invcellsize(0,0,0),
m_elementcount(0)
{
    Vec3st dims(1,1,1);
    Vec3d xmin(0,0,0), xmax(1,1,1);
    set(dims, xmin, xmax);
}

// --------------------------------------------------------
///
/// Calls assignment operator, which does a deep copy.
///
// --------------------------------------------------------

AccelerationGrid::AccelerationGrid(AccelerationGrid& other) :
m_cells(0,0,0),
m_elementidxs(0),
m_elementxmins(0),
m_elementxmaxs(0),
m_elementquery(0),
m_lastquery(0),
m_gridxmin(0,0,0),
m_gridxmax(0,0,0),
m_cellsize(0,0,0),
m_invcellsize(0,0,0),
m_elementcount(0)
{
    
    // Call assignment operator
    *this = other;
    
}

// --------------------------------------------------------
///
/// Deep copy.
///
// --------------------------------------------------------

AccelerationGrid& AccelerationGrid::operator=( const AccelerationGrid& other)
{
    m_cells.resize( other.m_cells.ni, other.m_cells.nj, other.m_cells.nk, 0 );
    for ( size_t i = 0; i < m_cells.a.size(); ++i )
    {
        if (other.m_cells.a[i])
        {
            m_cells.a[i] = new std::vector<size_t>();
            *(m_cells.a[i]) = *(other.m_cells.a[i]);
        }
    }
    
    m_elementcount = other.m_elementcount;
    m_elementidxs = other.m_elementidxs;
    m_elementxmins = other.m_elementxmins;
    m_elementxmaxs = other.m_elementxmaxs;
    m_elementquery = other.m_elementquery;
    m_lastquery = other.m_lastquery;
    m_gridxmin = other.m_gridxmin;
    m_gridxmax = other.m_gridxmax;
    m_cellsize = other.m_cellsize;
    m_invcellsize = other.m_invcellsize;   
    
    return *this;
}

// --------------------------------------------------------
///
/// Destructor: clear all grids
///
// --------------------------------------------------------

AccelerationGrid::~AccelerationGrid()
{
    clear();
}

// --------------------------------------------------------
///
/// Define the grid, given the extents of the domain and the number of desired voxels along each dimension
///
// --------------------------------------------------------

void AccelerationGrid::set( const Vec3st& dims, const Vec3d& xmin, const Vec3d& xmax )
{
    m_gridxmin = xmin;
    m_gridxmax = xmax;
    
    for(unsigned int i = 0; i < 3; i++)
    {
        m_cellsize[i] = (m_gridxmax[i]-m_gridxmin[i])/dims[i];
        m_invcellsize[i] = 1.0 / m_cellsize[i];
    }
    
    clear();
    
    m_cells.resize(dims[0], dims[1], dims[2]);    
    for(size_t i = 0; i < m_cells.a.size(); i++)
    {
        m_cells.a[i] = 0;
    }   
}

// --------------------------------------------------------
///
/// Generate a set of voxel indices from a pair of AABB extents
///
// --------------------------------------------------------

void AccelerationGrid::boundstoindices(const Vec3d& xmin, const Vec3d& xmax, Vec3i& xmini, Vec3i& xmaxi)
{
    
    xmini[0] = (int) std::floor((xmin[0] - m_gridxmin[0]) * m_invcellsize[0]);
    xmini[1] = (int) std::floor((xmin[1] - m_gridxmin[1]) * m_invcellsize[1]);
    xmini[2] = (int) std::floor((xmin[2] - m_gridxmin[2]) * m_invcellsize[2]);
    
    xmaxi[0] = (int) std::floor((xmax[0] - m_gridxmin[0]) * m_invcellsize[0]);
    xmaxi[1] = (int) std::floor((xmax[1] - m_gridxmin[1]) * m_invcellsize[1]);
    xmaxi[2] = (int) std::floor((xmax[2] - m_gridxmin[2]) * m_invcellsize[2]);
    
    if(xmini[0] < 0) xmini[0] = 0;
    if(xmini[1] < 0) xmini[1] = 0;
    if(xmini[2] < 0) xmini[2] = 0;
    
    if(xmaxi[0] < 0) xmaxi[0] = 0;
    if(xmaxi[1] < 0) xmaxi[1] = 0;
    if(xmaxi[2] < 0) xmaxi[2] = 0;
    
    assert( m_cells.ni < INT_MAX );
    assert( m_cells.nj < INT_MAX );
    assert( m_cells.nk < INT_MAX );
    
    if(xmaxi[0] >= (int)m_cells.ni) xmaxi[0] = (int)m_cells.ni-1;
    if(xmaxi[1] >= (int)m_cells.nj) xmaxi[1] = (int)m_cells.nj-1;
    if(xmaxi[2] >= (int)m_cells.nk) xmaxi[2] = (int)m_cells.nk-1;
    
    if(xmini[0] >= (int)m_cells.ni) xmini[0] = (int)m_cells.ni-1;
    if(xmini[1] >= (int)m_cells.nj) xmini[1] = (int)m_cells.nj-1;
    if(xmini[2] >= (int)m_cells.nk) xmini[2] = (int)m_cells.nk-1;
    
}

// --------------------------------------------------------
///
/// Add an object with the specified index and AABB to the grid
///
// --------------------------------------------------------

void AccelerationGrid::add_element(size_t idx, const Vec3d& xmin, const Vec3d& xmax)
{
    
   //if(m_elementxmins.capacity() <= idx) {
   //   //double the capacities (does resize do this behind the scenes, or only increase by one?)
   //   m_elementidxs.reserve(m_elementcount*2);
   //   m_elementxmins.reserve(m_elementcount*2);
   //   m_elementxmaxs.reserve(m_elementcount*2);
   //   m_elementquery.reserve(m_elementcount*2);
   //}

    if(m_elementcount <= idx)
    {
        if(m_elementidxs.size() <= idx) {
           m_elementidxs.resize(idx+1); //only ever grow m_elementidxs. but we won't clear it, since we don't want to have to reallocate its contained vectors.
           m_elementidxs[idx].reserve(10); //reserve some space in the vector
        }
        m_elementidxs[idx].clear(); //make sure the ``new" one is clear, if it already existed.

        m_elementxmins.resize(idx+1);
        m_elementxmaxs.resize(idx+1);
        m_elementquery.resize(idx+1);
        m_elementcount = idx+1;
    }
    
    
    m_elementxmins[idx] = xmin;
    m_elementxmaxs[idx] = xmax;
    m_elementquery[idx] = 0;
        
    Vec3i xmini, xmaxi;
    boundstoindices(xmin, xmax, xmini, xmaxi);
    
    Vec3st cur_index;
    for(cur_index[2] = xmini[2]; cur_index[2] <= xmaxi[2]; cur_index[2]++)
    {
        for(cur_index[1] = xmini[1]; cur_index[1] <= xmaxi[1]; cur_index[1]++)
        {
           for(cur_index[0] = xmini[0]; cur_index[0] <= xmaxi[0]; cur_index[0]++)
           {
                std::vector<size_t>*& cell = m_cells(cur_index[0], cur_index[1], cur_index[2]);
                if(!cell) {
                    cell = new std::vector<size_t>();
                    cell->reserve(10);
                }
                
                cell->push_back(idx);
                m_elementidxs[idx].push_back(cur_index);
            }
        }
    }
}

// --------------------------------------------------------
///
/// Remove an object with the specified index from the grid
///
// --------------------------------------------------------

void AccelerationGrid::remove_element(size_t idx)
{
    
    if ( idx >= m_elementcount ) { return; }
    
    for(size_t c = 0; c < m_elementidxs[idx].size(); c++)
    {
        Vec3st cellcoords = m_elementidxs[idx][c];
        std::vector<size_t>* cell = m_cells(cellcoords[0], cellcoords[1], cellcoords[2]);
        
        std::vector<size_t>::iterator it = cell->begin();
        while(*it != idx)
        {
            it++;
        }
        
        cell->erase(it);
    }
    
    m_elementidxs[idx].clear();
    
}


// --------------------------------------------------------
///
/// Reset the specified object's AABB
///
// --------------------------------------------------------
bool boxes_overlap(Vec3i low_0, Vec3i high_0, Vec3i low_1, Vec3i high_1) {
   if (high_0[0] < low_1[0])  return false; // a is left of b
   if (low_0[0]  > high_1[0]) return false; // a is right of b
   if (high_0[1] < low_1[1]) return false; // a is above b
   if (low_0[1]  > high_1[1]) return false; // a is below b
   if (high_0[2] < low_1[2]) return false; // a is in front of b
   if (low_0[2]  > high_1[2]) return false; // a is behind b
   return true; // boxes overlap
}

void AccelerationGrid::update_element(size_t idx, const Vec3d& xmin, const Vec3d& xmax)
{
   
   remove_element(idx);
   add_element(idx, xmin, xmax);
   return;
   /*
   //do an add and remove simultaneously (before we had done them separately, slower)

   bool is_new = false;
   if(m_elementcount <= idx)
   {
      if(m_elementidxs.size() <= idx) {
         m_elementidxs.resize(idx+1); //only ever grow m_elementidxs. but we won't clear it, since we don't want to have to reallocate its contained vectors.
         m_elementidxs[idx].reserve(30); //reserve some space in the vector
      }
      m_elementidxs[idx].clear(); //make sure the ``new" one is clear, if it already existed.

      m_elementxmins.resize(idx+1);
      m_elementxmaxs.resize(idx+1);
      m_elementquery.resize(idx+1);
      m_elementcount = idx+1;
      
      //initialize to nothingness (inverted box)
      m_elementxmins[idx] = Vec3d(-1,-1,-1);
      m_elementxmaxs[idx] = Vec3d(-2,-2,-2);
      m_elementquery[idx] = 0;
      is_new = true;
   }
   
   
   Vec3d xmin_old(0,0,0), xmax_old(0,0,0); 
   Vec3i xmini_new, xmaxi_new;
   Vec3i xmini_old, xmaxi_old;

   if(!is_new) {
      xmin_old = m_elementxmins[idx];
      xmax_old = m_elementxmaxs[idx];
   
      //get old and new index bounds
      boundstoindices(xmin,     xmax,     xmini_new, xmaxi_new);
      boundstoindices(xmin_old, xmax_old, xmini_old, xmaxi_old);
   }

   //try to do something smarter if the element has only moved slightly
   if(!is_new && boxes_overlap(xmini_old, xmaxi_old, xmini_new, xmaxi_new)) {
      
      //determine union of the two boxes
      Vec3i total_min = min_union(xmini_old, xmini_new);
      Vec3i total_max = max_union(xmaxi_old, xmaxi_new);
      
      m_elementxmins[idx] = xmin;
      m_elementxmaxs[idx] = xmax;
      m_elementquery[idx] = 0;

      //iterate over all the cells, old and new, updating as needed 
      Vec3st cur_index;
      for(cur_index[2] = total_min[2]; cur_index[2] <= total_max[2]; cur_index[2]++)
      {
         bool in_new_z = cur_index[2] >= xmini_new[2] && cur_index[2] <= xmaxi_new[2];
         bool in_old_z = cur_index[2] >= xmini_old[2] && cur_index[2] <= xmaxi_old[2];
         if(!in_new_z && !in_old_z) continue;

         for(cur_index[1] = total_min[1]; cur_index[1] <= total_max[1]; cur_index[1]++)
         {
            bool in_new_y = cur_index[1] >= xmini_new[1] && cur_index[1] <= xmaxi_new[1];
            bool in_old_y = cur_index[1] >= xmini_old[1] && cur_index[1] <= xmaxi_old[1];
            if(!in_new_y && !in_old_y) continue;

            for(cur_index[0] = total_min[0]; cur_index[0] <= total_max[0]; cur_index[0]++)
            {
               bool in_new_x = cur_index[0] >= xmini_new[0] && cur_index[0] <= xmaxi_new[0];
               bool in_old_x = cur_index[0] >= xmini_old[0] && cur_index[0] <= xmaxi_old[0];
               if(!in_new_x && !in_old_x) continue;

               bool in_new = in_new_x && in_new_y && in_new_z;
               bool in_old = in_old_x && in_old_y && in_old_z;
               
               assert(in_new || in_old);

               if(in_new) //in new set 
               {
                  if(!in_old) { //not in old, we need to add it
                     std::vector<size_t>*& cell = m_cells(cur_index[0], cur_index[1], cur_index[2]);
                     if(!cell) {
                        cell = new std::vector<size_t>();
                        cell->reserve(10);
                     }

                     cell->push_back(idx);
                     m_elementidxs[idx].push_back(cur_index);
                  }
                  //else: in both new and old, so we don't need to change anything!
               }
               else if(in_old) {//not in new set, but it is in the old set; must delete it
                  std::vector<size_t>*& cell = m_cells(cur_index[0], cur_index[1], cur_index[2]);

                  //erase the index of the element in the cell
                  std::vector<size_t>::iterator it = cell->begin();
                  while(*it != idx)
                     it++;
                  cell->erase(it);

                  //erase the index of the *cell* in the *element* -> is this pricy?
                  std::vector<Vec3st>::iterator it2 = m_elementidxs[idx].begin();
                  while(*it2 != cur_index) 
                     it2++;
                  m_elementidxs[idx].erase(it2);

               }
            }
         }
      }
   }
   else { //the old and new regions don't overlap, so do it the dumb way

      if(!is_new) {
         //erase all the old ones
         for(size_t c = 0; c < m_elementidxs[idx].size(); c++)
         {
            Vec3st cellcoords = m_elementidxs[idx][c];
            std::vector<size_t>*& cell = m_cells(cellcoords[0], cellcoords[1], cellcoords[2]);

            std::vector<size_t>::iterator it = cell->begin();
            while(*it != idx)
            {
               it++;
            }

            cell->erase(it);
         }
      }

      //update the data
      m_elementidxs[idx].clear();
      m_elementxmins[idx] = xmin;
      m_elementxmaxs[idx] = xmax;
      m_elementquery[idx] = 0;

      Vec3i xmini, xmaxi;
      boundstoindices(xmin, xmax, xmini, xmaxi);

      //now add the geometry back into the acceleration structure
      Vec3st cur_index;
      for(cur_index[2] = xmini[2]; cur_index[2] <= xmaxi[2]; cur_index[2]++)
      {
         for(cur_index[1] = xmini[1]; cur_index[1] <= xmaxi[1]; cur_index[1]++)
         {
            for(cur_index[0] = xmini[0]; cur_index[0] <= xmaxi[0]; cur_index[0]++)
            {
               std::vector<size_t>*& cell = m_cells(cur_index[0], cur_index[1], cur_index[2]);
               if(!cell) {
                  cell = new std::vector<size_t>();
                  cell->reserve(10);
               }

               cell->push_back(idx);
               m_elementidxs[idx].push_back(cur_index);
            }
         }
      }
   }
   */

}

// --------------------------------------------------------
///
/// Remove all elements from the grid
///
// --------------------------------------------------------

void AccelerationGrid::clear()
{
    for(size_t i = 0; i < m_cells.a.size(); i++)
    {
        std::vector<size_t>*& cell = m_cells.a[i];  
        if(cell)
        {
            delete cell;
            cell = 0;
           //cell->clear(); //don't clear the memory, since we likely want to reuse it.
        }
    }
    
    //clear the entries for each element, but don't clear the overall vector itself, since we
    //don't want to reallocate memory for the individual vectors.
    //m_elementidxs.clear();
    for(int i = 0; i < m_elementidxs.size(); ++i) { 
       m_elementidxs[i].clear();
    }

    m_elementxmins.clear();
    m_elementxmaxs.clear();
    m_elementquery.clear();
    m_lastquery = 0;

    m_elementcount = 0;
    
}

// --------------------------------------------------------
///
/// Return the set of elements which have AABBs overlapping the query AABB.
///
// --------------------------------------------------------

void AccelerationGrid::find_overlapping_elements( const Vec3d& xmin, const Vec3d& xmax, std::vector<size_t>& results ) 
{
    if(m_lastquery == std::numeric_limits<unsigned int>::max())
    {
        std::vector<unsigned int>::iterator iter = m_elementquery.begin();
        for( ; iter != m_elementquery.end(); ++iter )
        {
            *iter = 0;
        }
        m_lastquery = 0;
    }
    
    ++m_lastquery;
    
    Vec3i xmini, xmaxi;
    boundstoindices(xmin, xmax, xmini, xmaxi);
    
    for(int k = xmini[2]; k <= xmaxi[2]; ++k)
    {
      for(int j = xmini[1]; j <= xmaxi[1]; ++j)
      {
         for(int i = xmini[0]; i <= xmaxi[0]; ++i)
         {

                std::vector<size_t>* cell = m_cells(i, j, k);
                
                if(cell)
                {
                    for( std::vector<size_t>::const_iterator citer = cell->begin(); citer != cell->end(); ++citer)
                    {
                        size_t oidx = *citer;
                        
                        // Check if the object has already been found during this query
                        
                        if(m_elementquery[oidx] < m_lastquery)
                        {
                            
                            // Object has not been found.  Set m_elementquery so that it will not be tested again during this query.
                            
                            m_elementquery[oidx] = m_lastquery;
                            
                            const Vec3d& oxmin = m_elementxmins[oidx];
                            const Vec3d& oxmax = m_elementxmaxs[oidx];
                            
                            if( (xmin[0] <= oxmax[0] && xmin[1] <= oxmax[1] && xmin[2] <= oxmax[2]) &&
                               (xmax[0] >= oxmin[0] && xmax[1] >= oxmin[1] && xmax[2] >= oxmin[2]) )
                            {
                                results.push_back(oidx);
                            }
                            
                        }
                    }
                }
                
            }
        }
    }
}

}

