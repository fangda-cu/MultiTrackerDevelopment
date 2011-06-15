// LevelSet.hh
//

#ifndef _LEVELSET_H_
#define _LEVELSET_H_

#include "../Core/EigenIncludes.hh"

#include <vector>
#include <string>
#include <fstream>

/*#include "common/Defines.hh"

#include "math/Vector.hh"

#include "common/vec.hh"
#include "common/array3.hh"*/

#include "Bridson/Vector.hh"
#include "Bridson/vec.hh"
#include "Bridson/array3.hh"
#include "Bridson/util.hh"

namespace BASim {

class LevelSet
{
public:
    LevelSet();
    ~LevelSet();

//    void buildLevelSet(const std::vector<bridson::Vec3ui> &triangles,
//                       const std::vector<bridson::Vec3f>  &x,
//                       const std::vector<bridson::Vec3f>  &v,
//                       const bridson::Vec3f &origin, Real length[3], int nbrTriangles, bool bound=true);
/*
    void buildLevelSet(const std::vector<bridson::Vec3ui> &triangles,
                       const std::vector<uint> &triIndices,
                       const std::vector<bridson::Vec3f>  &x,
                       const std::vector<bridson::Vec3f>  &v,
                       const bridson::Vec3f &origin, Real length[3],
                       Real dx, int nx, int ny, int nz, int nbrTriangles);
*/
    void buildLevelSet(const Vec3Indices &triangles,
                       const Indices &triIndices,
                       const std::vector<bridson::Vec3f>  &x,
                       const std::vector<bridson::Vec3f>  &v,
                       const bridson::Vec3f &origin, Real length[3],
                       Real dx, int nx, int ny, int nz, int nbrTriangles,
                       Eigen::Matrix4f& transformMatrix);

    Real getLevelSetValue(Vec3<Real> x);
    Real getLevelSetValueVelocity(Vec3<Real> &x, Vec3<Real> &v);
    void getGradient(Vec3<Real> &x, Vec3<Real> &grad);

    bridson::Vec3f& getOrigin() { return _origin; }

    Real getGridSize() { return _dx; }

    int getNbrX() { return _phi.ni; }
    int getNbrY() { return _phi.nj; }
    int getNbrZ() { return _phi.nk; }

    bridson::Array3f& getPhi() { return _phi; }
    bridson::Array3<bridson::Vec3f, bridson::Array1<bridson::Vec3f> >& getPhiVel() { return _phiVel; }

    void draw();

    void writeFile(std::fstream &levelSetFile);
    void loadFile(std::fstream &levelSetFile);

    bool isInitialized() { return _initialized; }

    // This stores the transformation matrix for the mesh the level set is created from
    void setTransformationMatrix( Eigen::Matrix4f& i_matrix );

protected:
    bridson::Array3f _phi;
    bridson::Array3<bridson::Vec3f, bridson::Array1<bridson::Vec3f> > _phiVel;

    bridson::Vec3f _origin;
    Real _dx;

    bool _initialized;

protected:
    void buildLevelSet(const Vec3Indices &triangles,
                       const Indices &triIndices,
                       const std::vector<bridson::Vec3f>  &x,
                       const std::vector<bridson::Vec3f>  &v,
                       const bridson::Vec3f &origin,
                       float dx,
                       int ni, int nj, int nk,
                       bridson::Array3f &phi,
                       bridson::Array3<bridson::Vec3f, bridson::Array1<bridson::Vec3f> > &phiVel);

    float point_triangle_distance(const bridson::Vec3f &p,
                                  const bridson::Vec3f &a, const bridson::Vec3f &b, const bridson::Vec3f &c,
                                  float &t1, float &t2, float &t3);

    void check_neighbour(const Vec3Indices &tri,
                         const std::vector<bridson::Vec3f> &x,
                         const std::vector<bridson::Vec3f> &v,
                         bridson::Array3f &phi,
                         bridson::Array3<bridson::Vec3f, bridson::Array1<bridson::Vec3f> > &phi_vel,
                         bridson::Array3i &closest_tri,
                         const bridson::Vec3f &gx,
                         int i0, int j0, int k0, int i1, int j1, int k1);

    void sweep(const Vec3Indices &tri,
               const std::vector<bridson::Vec3f> &x,
               const std::vector<bridson::Vec3f> &v,
               bridson::Array3f &phi,
               bridson::Array3<bridson::Vec3f, bridson::Array1<bridson::Vec3f> > &phi_vel,
               bridson::Array3i &closest_tri,
               const bridson::Vec3f &origin,
               float dx,
               int di, int dj, int dk);

    int orientation(double x1, double y1, double x2, double y2, double &twice_signed_area);

    bool point_in_triangle_2d(double x0, double y0, 
                              double x1, double y1, double x2, double y2, double x3, double y3,
                              double& a, double& b, double& c);

private:
    // The mesh may not be at the origin when it is created. So remember where it was so that
    // later when we're asked to sample the level set we can remove the initial transform
    // before working out the true place to ask.
    Eigen::Matrix4f m_transformMatrixAtCreation;
    Eigen::Matrix4f m_currentTransformMatrix;
    
    // Debug data to check whether the values returned work for meshes that move or
    // were created with a transform on them
    std::vector< Eigen::Vector4f > m_realRequestPositions;
    std::vector< Eigen::Vector4f > m_transformedRequestPositions;
    std::vector< Eigen::Vector4f > m_grad;
    std::vector< Eigen::Vector4f > m_gradPosition;
};

}

#endif

