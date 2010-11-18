/* 
 *
 * See http://twiki.wetafx.co.nz/Weta/VolumetricHairCollisions for more documentation
 *
 */



#ifndef __VOLUMETRIC_COLLISIONS_CPU_HH__
#define __VOLUMETRIC_COLLISIONS_CPU_HH__

// Maya #defines this I think and Physbam has an MPI class called Status that
// gets screwed up on compilation
#undef Status 

#include "VolumetricCollisions.hh"

#include <vector>
#include <iostream>
#include <cmath>

#include <Geometry/SEGMENTED_CURVE.h>
#include <Grids/GRID_3D.h>
#include <Incompressible_Flows/INCOMPRESSIBLE_UNIFORM.h>
#include <Matrices_And_Vectors/VECTOR.h>
#include <Grids/UNIFORM_GRID_ITERATOR_FACE_3D.h>
#include <Grids/UNIFORM_GRID_ITERATOR_CELL_3D.h>

using namespace std;

class VolumetricCollisionsCPU : public VolumetricCollisions
{

    typedef PhysBAM::VECTOR<double,3> TV;
    typedef PhysBAM::GRID_3D<double> GRID_3D;
    typedef typename GRID_3D::CELL_ITERATOR pbCellIterator;
    typedef typename GRID_3D::FACE_ITERATOR pbFaceIterator;

public:
    VolumetricCollisionsCPU();
    explicit VolumetricCollisionsCPU( RodDataMap &rodDataMap );
    ~VolumetricCollisionsCPU();

    void respondVolumetricCollisions( RodDataMap &rodDataMap, Real targetEdgeDensity, 
                                      Real volumetricRadius, Real gridDx, 
                                      BASim::Vec3d separationCondition,  double flip,
                                      double slip, CollisionMeshDataHashMap &collisionMeshes );

    void draw( bool displayGrid, Real displayGridVelocitiesMultiplier, Real maxDisplayDensity,
               bool displayPsiN, bool displayPsiD );

private:
    void reinitialiseRods( RodDataMap &rodDataMap );
    void deleteRods();

    ElasticRod*  initialiseRodMap( RodDataMap& i_rodDataMap );
    ElasticRod* nextRod();

    bool _initialised;
    GRID_3D _grid;
    PhysBAM::SEGMENTED_CURVE<TV>* _hairCurve;
    PhysBAM::PROJECTION_UNIFORM<GRID_3D>* _projection;

    PhysBAM::FACE_ARRAYS_3D<double> _faceWeights, _initialVelocities, _deltas;
    PhysBAM::ARRAYS_3D<double> _cellWeights;
    PhysBAM::ARRAYS_3D<TV> _cellVelocities;
    PhysBAM::ARRAYS_3D<bool> _separatingCells;

    int m_stepNumber;

    RodDataMap* m_rodDataMap;
    RodDataMapIterator m_rodDataMapIterator;    
    int m_rodIndex;
    double m_flip;
    double m_slip;
};

#endif //__VOLUMETRIC_COLLISIONS_CPU_HH__
