/**
 * \file ShellStickyRepulsionForce.hh
 *
 * \author batty@cs.columbia.edu
 * \date Nov 22, 2011
 */

#ifndef SHELLSTICKYREPULSIONFORCE_HH
#define SHELLSTICKYREPULSIONFORCE_HH

#include "BASim/src/Physics/DeformableObjects/Shells/ElasticShellForce.hh"

#include "BASim/src/Math/ADT/adreal.h"
#include "BASim/src/Math/ADT/advec.h"
#include <set>

//A set of springs, each between a tri and a vert.

namespace BASim {


typedef Scalar Real;
typedef CVec3T<Real> Vector3d;

const int NumRepulsionDof = 12;
const int NumRepulsionVerts = 4;

class ShellStickyRepulsionForce : public ElasticShellForce {

public:

  ShellStickyRepulsionForce (ElasticShell& shell, const std::string& name = "ShellStickyRepulsionForce", Scalar timestep = 1.0);
  virtual ~ShellStickyRepulsionForce () {}

  std::string getName() const;

  Scalar globalEnergy() const;
  void globalForce(VecXd& force) const;
  void globalJacobian(Scalar scale, MatrixBase& Jacobian) const;
  
  void addSpring(const FaceHandle& fh, const VertexHandle& vh, const Vec3d& baryCoords, Scalar stiffness, Scalar damping, Scalar restlen);
  void clearSprings();
  void getSpringLists(std::vector<VertexHandle> &verts, std::vector<FaceHandle>& tris, std::vector<Vec3d>& barycoords);
  bool springExists(const FaceHandle& f, const VertexHandle& v);

  void clearSprings(VertexHandle& v);
  void clearSprings(FaceHandle& f);

protected:

  bool gatherDOFs(const FaceHandle& fh, const VertexHandle& vh, std::vector<Vec3d>& deformed, std::vector<Vec3d>& undeformed_damp, std::vector<int>& indices) const;

  Scalar elementEnergy(const std::vector<Vec3d>& deformed, const Vec3d& baryCoords, Scalar strength, Scalar restlen) const;
  void elementForce(const std::vector<Vec3d>& deformed, const Vec3d& baryCoords, Scalar strength, Scalar restlen,
                    Eigen::Matrix<Scalar, NumRepulsionDof , 1>& force) const;
  void elementJacobian(const std::vector<Vec3d>& deformed, const Vec3d& baryCoords, Scalar strength, Scalar restlen,
                       Eigen::Matrix<Scalar, NumRepulsionDof , NumRepulsionDof >& J) const;
  
  //List of springs
  std::vector<FaceHandle> m_faces;
  std::vector<VertexHandle> m_vertices;
  std::vector<Vec3d> m_barycoords;
  std::vector<Scalar> m_stiffnesses;
  std::vector<Scalar> m_damping;
  std::vector<Scalar> m_restlen;

  //for fast checking to see if a given spring already exists
  std::set< std::pair<int,int> > m_springset; 
  

  Scalar m_timestep; //for damping/viscosity
};




}


#endif //SHELLSTICKYREPULSIONFORCE_HH
