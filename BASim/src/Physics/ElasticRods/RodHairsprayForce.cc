#include "RodHairsprayForce.hh"

namespace BASim {

RodHairsprayForce::RodHairsprayForce(ElasticRod &rod, vector<Scalar>& ks, vector<Vec3d>& curvePositions)
  : m_rod(rod), m_curvePositions(curvePositions)
{
  m_ks.resize(ks.size());
  m_ds.resize(ks.size());
  for (size_t i=0; i<ks.size(); ++i)
  {
      m_ks[i] = ks[i];
      m_ds[i] = 2.0 * std::sqrt(m_ks[i] * rod.getVertexMass(i));
  }
}

RodHairsprayForce::~RodHairsprayForce()
{
}

void RodHairsprayForce::computeForce(const ElasticRod& rod, VecXd& F)
{
  // We get a rod passed in but we are assuming it's the same rod we were given when
  // initialising the class. If not we're screwed. This quirk comes because we are
  // creating the hairspray force as an external force like gravity when it is
  // actually attached to a specific rod.

  for (int i=0; i<rod.nv(); ++i)
  {
    if (m_ks[i] == 0 || m_rod.vertFixed(i))
      continue;
    else
      cerr << "Calculating vertex " << i << "\n";
    
    Vec3d f = m_ks[i] * (m_curvePositions[i] - m_rod.getVertex(i)) - m_ds[i] * m_rod.getVelocity(i);

    if (!rod.vertFixed(i))
      for (int k=0; k<3; ++k)
        F[rod.vertIdx(i, k)] += f[k]; 
  }
}

}
