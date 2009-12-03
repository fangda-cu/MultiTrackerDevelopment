// UniformGrid.hh
//
// This class is (lovingly) ripped of from code on Robert Bridson's
// homepage. Modified slightly to use our data types, except for the Array
// type
//


#ifndef UNIFORMGRID_H
#define UNIFORMGRID_H

#include <tr1/unordered_map>
#include"CandidateCollision.hh"
//#include <BASim/src/Physics/ElasticRods/ElasticRod.hh>
#include "BASim/src/Physics/ElasticRods/Stencil.hh"

#include "BASim/src/Physics/DegreeOfFreedom.hh"
#include "BASim/src/Physics/PhysObject.hh"
#include "BASim/src/Physics/ElasticRods/ElasticRod.hh"


namespace BASim {

typedef tr1::unordered_map<int, std::pair<int, ElasticRod *> > EdgeIndexMap;
typedef tr1::unordered_map<int, std::pair<int, ElasticRod *> >::iterator EdgeIndexMapIterator;
    
class UniformGrid
{
public:
    UniformGrid();
    ~UniformGrid();

    void getProximities(std::vector<ElasticRod *> &rods, Collisions &collisions);
    void getContinuousTimeCollisions(std::vector<ElasticRod*> &rods, Real dt, Collisions &collisions);

    void addElement(unsigned int idx, Vec3d &xmin, Vec3d &xmax);
//    void updateElement(unsigned int idx, Vec3d &xmin, Vec3d &xmax);
//    void removeElement(unsigned int idx);

    void set(Vec3i &dims, Vec3d &xmin, Vec3d &xmax);

    void clear();
    
    void findOverlappingElements(Vec3d &xmin, Vec3d &xmax, std::vector<uint> &idxs);

    void writeFile(std::string filename);
    bool readFile(std::string filename);

    void getGridDims(Vec3d &xmin, Vec3d &xmax, Vec3d &dX)
    { xmin = gridxmin; xmax = gridxmax; dX = cellsize;}
    
    inline void addScalertoVec3d( Vec3d& x, Scalar s )
    {
        for (unsigned int i=0; i<3; ++i)
            x[i] += s;
    }
    
    inline Vec3d divideVec3dByVec3d( const Vec3d& x, const Vec3d& y )
    {
        Vec3d o;
        for (unsigned int i=0; i<3; ++i)
            o[i] = x[i] / y[i];
        return o;
    }
    
    inline Vec3i floor(const Vec3d &a)
    { 
        Vec3i rounded;
        for (unsigned int i=0; i<3; ++i)
            rounded[i] = (int)std::floor(a[i]);
    
        return rounded; 
    }
    
    bool isLessThanOrEqual(const Vec3d &a, const Vec3d &b)
    { 
       for(unsigned int i = 0; i < 3; i++)
          if(a[i] > b[i])
             return false;
       return true;
    }
    
    bool isGreaterThanOrEqual(const Vec3d &a, const Vec3d &b)
    { 
       for(unsigned int i = 0; i < 3; i++)
          if(a[i] < b[i])
             return false;
       return true;
    }


protected:
    void boundsToIndices(Vec3d &xmin, Vec3d &xmax, Vec3i &xmini, Vec3i &xmaxi);

protected:
    bridson::Array<std::vector<unsigned int>* > cells;
//    std::vector<std::vector<Vec3ui> > elementidxs;
    std::vector<Vec3d> elementxmins;
    std::vector<Vec3d> elementxmaxs;
    std::vector<unsigned int> elementquery;
    unsigned int lastquery;
    Vec3d gridxmin, gridxmax;
    Vec3d cellsize;

private:

};

}

#endif

