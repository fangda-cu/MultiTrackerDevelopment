/**
 * \file RodForce.hh
 *
 * \author miklos@cs.columbia.edu
 * \date 09/01/2009
 */

#ifndef RODFORCE_HH
#define RODFORCE_HH

#ifdef WETA
#include "Stencil.hh"
#endif
#include "ElasticRod.hh"

namespace BASim
{

/** Base class for a force that acts on rods. */
class RodForce
{
public:

  typedef ElasticRod::vertex_handle vertex_handle;
  typedef ElasticRod::edge_handle   edge_handle;

  explicit RodForce(ElasticRod& rod, const std::string& name = "RodForce");
  virtual ~RodForce() {}

  std::string getName() const;

  void computeKb(Vec3d& kb, const Vec3d& x0, const Vec3d& x1, const Vec3d& x2);
  Vec3d computeKb(const Vec3d& x0, const Vec3d& x1, const Vec3d& x2);

  void computeDkb(Mat3dArray& Dkb, const Vec3d& x0, const Vec3d& x1,
                  const Vec3d& x2);
  Mat3dArray computeDkb(const Vec3d& x0, const Vec3d& x1, const Vec3d& x2);

  virtual Scalar globalEnergy() = 0;
  virtual void globalForce(VecXd& force) = 0;
  virtual void globalJacobian(int baseidx, Scalar scale, MatrixBase& Jacobian) = 0;
  
  // TODO(sainsley) : RENAME THESE METHODS TO INDICATE THEY ARE ACCUMULATORS
  // TODO(sainsley) : remove default without getting compile errors for non-sym forces
  virtual void globalForceEnergy(VecXd& force, Scalar& energy)
  {
    // more efficient implementations are strongly encouraged!!!
    globalForce(force);
    energy = globalEnergy();
  }
  virtual void globalJacobianForceEnergy(int baseidx, Scalar scale, MatrixBase& Jacobian, 
					 VecXd& force, Scalar& energy) {

    globalJacobian(baseidx, scale, Jacobian);
    globalForce(force);
    energy = globalEnergy();
  }

  virtual void updateProperties() {}
  virtual void updateStiffness() {}
  virtual void updateUndeformedStrain() {}
  virtual void updateReferenceDomain() {}
  
  virtual void setReferenceLengths(std::vector<Scalar>& vals) {}

  virtual void verifyProperties() {}

  virtual void updatePlasticity(Scalar maxKappa) {}
  
  virtual void updateUndeformedConfiguration(std::vector<Scalar>& vals) {}
  
  virtual void globalReverseJacobian(MatrixBase& Jacobian) {}
  virtual void updateReverseUndeformedStrain(const VecXd& e) {}
  
  // 'Attaches' any links to properties within children. Used for serialization.
  virtual void reattatchProperties() {};

  bool viscous() const { return m_viscous; }
  void setViscous(bool v) { m_viscous = v; updateStiffness(); }

protected:

    ElasticRod& m_rod;
    std::string m_name;
    bool m_viscous;

    static Mat2d J;
    static Mat2d Jt;
};

template<class Stencil>
class RodForceT: public RodForce
{
public:

    typedef typename Stencil::iterator iterator;

    explicit RodForceT(ElasticRod& rod, const std::string& name = "RodForce") :
        RodForce(rod, name), m_stencil(rod)
    {
    }

protected:

    Stencil m_stencil;
};

} // namespace BASim

#endif // RODFORCE_HH
