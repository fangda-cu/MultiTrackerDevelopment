// LevelSet.hh
//

#ifndef _LEVELSET_H_
#define _LEVELSET_H_

#include "CollisionObject.hh"

#include <vector>
#include <string>
#include <fstream>

#include "BridsonArray3.hh"
#include "BridsonVec.hh"

#include "RodneyVector.hh"


namespace BASim {

//#include "common/vec.hh"
// 
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
                       const std::vector<BASim::Vec3d>  &x,
                       const std::vector<BASim::Vec3d>  &v,
                       const BASim::Vec3d &origin, Real length[3],
                       Real dx, int nx, int ny, int nz, int nbrTriangles);



    Real getLevelSetValue(Vec3d x);
    Real getLevelSetValueVelocity(BASim::Vec3d &x, BASim::Vec3d &v);
    void getGradient(Vec3<Real> &x, Vec3<Real> &grad);

    Vec3d& getOrigin() { return _origin; }

    Real getGridSize() { return _dx; }

    int getNbrX() { return _phi.ni; }
    int getNbrY() { return _phi.nj; }
    int getNbrZ() { return _phi.nk; }

    bridson::Array3f& getPhi() { return _phi; }
    bridson::Array3<BASim::Vec3d, bridson::Array1<BASim::Vec3d> >& getPhiVel() { return _phiVel; }

    void draw();

    void writeFile(std::fstream &levelSetFile);
    void loadFile(std::fstream &levelSetFile);

    bool isInitialized() { return _initialized; }

protected:
    bridson::Array3f _phi;
    bridson::Array3<BASim::Vec3d, bridson::Array1<BASim::Vec3d> > _phiVel;

    Vec3d _origin;
    Real _dx;

    bool _initialized;

protected:
    void buildLevelSet(const Vec3Indices &triangles,
                       const Indices &triIndices,
                       const std::vector<BASim::Vec3d>  &x,
                       const std::vector<BASim::Vec3d>  &v,
                       const Vec3d &origin,
                       float dx,
                       int ni, int nj, int nk,
                       bridson::Array3f &phi,
                       bridson::Array3<BASim::Vec3d, bridson::Array1<BASim::Vec3d> > &phiVel);

    float point_triangle_distance(const BASim::Vec3d &p,
                                  const BASim::Vec3d &a, const BASim::Vec3d &b, const BASim::Vec3d &c,
                                  float &t1, float &t2, float &t3);

    void check_neighbour(const Vec3Indices &tri,
                         const std::vector<BASim::Vec3d> &x,
                         const std::vector<BASim::Vec3d> &v,
                         bridson::Array3f &phi,
                         bridson::Array3<BASim::Vec3d, bridson::Array1<BASim::Vec3d> > &phi_vel,
                         bridson::Array3i &closest_tri,
                         const BASim::Vec3d &gx,
                         int i0, int j0, int k0, int i1, int j1, int k1);

    void sweep(const Vec3Indices &tri,
                             const std::vector<Vec3d> &x,
                             const std::vector<Vec3d> &v,
                             bridson::Array3f &phi,
                             bridson::Array3<Vec3d, bridson::Array1<Vec3d> > &phi_vel,
                             bridson::Array3i &closest_tri,
                             const Vec3d &origin,
                             float dx,
                             int di, int dj, int dk);

    int orientation(double x1, double y1, double x2, double y2, double &twice_signed_area);

    bool point_in_triangle_2d(double x0, double y0, 
                              double x1, double y1, double x2, double y2, double x3, double y3,
                              double& a, double& b, double& c);

private:

};

}

#endif

