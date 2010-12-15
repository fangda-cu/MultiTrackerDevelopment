// CollisionMeshData.hh
//

#ifndef COLLISIONMESHDATA_H
#define COLLISIONMESHDATA_H

#include "UniformGrid.hh"
#include "../Math/Math.hh"
#include <string>
#include <tr1/unordered_map>

#include "../Core/Definitions.hh"

//#include "BoundingVolumeTree.hh"
#include "LevelSet.hh"
//#include "AdaptiveLevelSet.hh"

namespace BASim {
    
class CollisionMeshData : public CollisionObject
{
public:
    CollisionMeshData();
    virtual ~CollisionMeshData();

    void initialize();

    void clearAll();

    void reset(vector<Vec3d>& points);
    void update(vector<Vec3d>& points, std::string filename, int currFrame);

    void interpolate(Real percentage);

    void draw();

//       BoundingVolumeTree& getTree() { return _bvTree; }

    // CollisionObject interface
    //
    
    virtual Positions& getPositions()
    { return currPositions; }
    
    virtual Velocities& getVelocities()
    { return velocities; }
    
    virtual Indices& getEdgeIndices()
    { return _edgeIndices; }

    virtual Indices& getTriangleIndices()
    { return triangleIndices; }

    virtual Real getMassInverse(uint vertex)
    { return 0.0; }

    virtual Real getFrictionCoefficient()
    { return _friction; }

    virtual Real getThickness()
    { return _thickness; }

    void setLevelsetDx(Real levelsetDx){ _levelsetDx = levelsetDx; }
    
    void setFPS(int fps){_fps = fps;}

    void setThickness(Real thickness) { _thickness = thickness; }

    void setFriction(Real friction) { _friction = friction; }

    void setSeparationStrength(Real ss) { _separationStrength = ss; }

    Real getSeparationStrength() const { return _separationStrength; }

    void setDamping(Real dmp) { _damping = dmp; }

    Real getDamping() const { return _damping; }

    void setCoefficientOfRestitution(Real cor) { _coefficientOfRestitution = cor; }

    Real getCoefficientOfRestitution() const { return _coefficientOfRestitution; }

    void setFullCollisions(bool fc) { _fullCollisions = fc; }

    bool getFullCollisions() const { return _fullCollisions; }

    // Level set functions
    //

    void sizeLevelSet( Vec3d &origin, bridson::Vec3ui &dims, Real &dx, Real length[3]);
    void buildLevelSet();
    void buildLevelSet(int currFrame);
    Real getLevelSetValue(Vec3d& x, Vec3d& v);
    void getGradient(Vec3d &x, Vec3f &grad);

    Vec3d& vertex(uint i);
    const Vec3d& vertex(uint i) const;

    Vec3d& velocity(uint i);
    const Vec3d& velocity(uint i) const;
    
    Vec3d& triangleVertex(uint faceID, uint i);
    const Vec3d& triangleVertex(uint faceID, uint i) const;

    Vec3d& triangleVelocity(uint faceID, uint i);
    const Vec3d& triangleVelocity(uint faceID, uint i) const;

    bool initialized()
    { return _initialized; }

    std::vector<uint> triangleIndices;

    // The previous positions from the last timestep and
    // the positions at the current timestep
    //
    Positions prevPositions;
    Positions currPositions;

    // The positions from the last frame and the positions
    // from the next frame
    //
    Positions oldPositions;
    Positions newPositions;

    vector<Positions> allPositions;

    // The velocities that take us from oldPositions to newPositions
    // in (1 / fps) time
    //
    Velocities velocities;
    
    // Maximum velocity per frame
    float m_maxVelocityMag;

    int _nbrTriangles;
    
    Indices _triIndices;
    Indices _edgeIndices;
    
    Vec3Indices _tri;
    Vec3Indices _triangleEdgeIndices;

//    BoundingVolumeTree _bvTree;
    UniformGrid _grid;

    Real _thickness;
    Real _friction;
    Real _separationStrength;
    Real _coefficientOfRestitution;
    Real _damping;

    bool _fullCollisions;
    bool _initialized;

    int _fps;

    // Level set data
    //
    Real _levelsetDx;
    Real _percent;
    std::vector< Vec3d > _x;
    std::vector< Vec3d > _v;
     
    LevelSet *_phiPrevious;
    LevelSet *_phiCurrent;
//    AdaptiveLevelSet *_phiPrevious;
//    AdaptiveLevelSet *_phiCurrent;

    bool recordToFile;
    int recordFrames;
    std::string recordFilename;
    
    void writeMeshesToFile();
    
protected:
    void updateGrid(vector<Vec3d>& points, std::string filename="");

};

}

#endif

