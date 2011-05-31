/**
 * \file RodMayaForces.hh
 *
 * \author acoull@wetafx.co.nz
 * \date 22/04/2011
 */

#ifndef RODMAYAFORCES_HH
#define RODMAYAFORCES_HH

namespace BASim {

/** External forces coming from Maya fields that get applied to the rod. */
class RodMayaForces : public RodExternalForce
{
public:

  /**
   * Constructor for the force, does nothing as the force is just passed through from Maya
   *
   **/
  explicit RodMayaForces( ElasticRod* i_elasticRod ) : m_elasticRod( i_elasticRod )
  {
    m_name = "maya";    
    resetExternalMayaForces();
  }

  /**
   * Computes the force from Maya for the given rod and adds it
   * to the given vector.
   *
   * \param[in] rod The rod to compute the gravity force for.
   * \param[out] force Vector storing the forces on the rod.
   */

  // Is it safe to store the rod from the constructor like I'm doing? 
  // Why is there a rod passed in here, is it because these external force
  // classes are supposed to be independant of rods? This one definitely 
  // isn't! - Alasdair
  void computeForce(const ElasticRod& rod, VecXd& force)
  {
    for (int v = 0; v < rod.nv(); ++v) 
    {
      for (int coord = 0; coord < 3; ++coord) 
      {
        force( rod.vertIdx( v, coord ) ) += m_forcesFromMaya[ v ][ coord ];
      }
    }
  }

  void computeForceDX(int baseidx, const ElasticRod& rod, Scalar scale, MatrixBase& J) {}

  void computeForceDV(int baseidx, const ElasticRod& rod, Scalar scale, MatrixBase& J) {}

  /** The forces come directly from Maya and we store them here then pass them to the stepper
   * in the computeForce() function when we are asked to do so.
   **/
  void addExternalForceToVertex( const size_t i_vertexIndex, const BASim::Vec3d i_force )
  {
    if ( m_forcesFromMaya.size() > i_vertexIndex )
    {
      m_forcesFromMaya[ i_vertexIndex ] += i_force;
    }
  }
    
  void resetExternalMayaForces()
  {
    if ( m_elasticRod )
    {
      m_forcesFromMaya.resize( m_elasticRod->nv() );
      for ( size_t v = 0; v < m_elasticRod->nv(); ++v )
      {
        m_forcesFromMaya[ v ] = BASim::Vec3d( 0.0, 0.0, 0.0 );
      }
    }
  }   

protected:

  Vec3d m_gravity;
  
  std::vector< Vec3d > m_forcesFromMaya;
  ElasticRod* m_elasticRod;
};

} // namespace BASim

#endif // RODGRAVITY_HH
