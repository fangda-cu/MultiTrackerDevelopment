// RodLevelSetForce.cc
//

#include "RodLevelSetForce.hh"
#ifdef WETA
#include "../../Math/Math.hh"
#else
#include "BASim/src/Math/Math.hh"
#endif

#include "../../Util/TextLog.hh"

#include "../../Collisions/LevelSet.hh"


using namespace BASim;

RodLevelSetForce::RodLevelSetForce( LevelSet* levelSet ) 
{
    m_name = "level set penalty";

    m_levelSet = levelSet;
    assert( m_levelSet );
    assert( m_levelSet->isInitialized() );
}

RodLevelSetForce::~RodLevelSetForce()
{
}


void RodLevelSetForce::computeForce(const ElasticRod& rod, VecXd& force) const
{
    Scalar dummyEnergy;
    computeForceEnergy(rod, force, dummyEnergy);
}


void RodLevelSetForce::computeForceDX(int baseindex, const ElasticRod& rod, Scalar scale, MatrixBase& J) const
{
    // explicit implementation (for now)
}


Scalar RodLevelSetForce::computeEnergy(const ElasticRod& rod) const
{
    Scalar energy;
    VecXd dummyforce;
    
    computeForceEnergy(rod, dummyforce, energy);
    return energy;
}

 
void RodLevelSetForce::computeForceEnergy(const ElasticRod& rod, VecXd& force, Scalar& energy) const
{
    Scalar stiffness = 1;

    for (int vidx = 0; vidx < (int) rod.nv(); vidx++)
    {
        Vec3d v0 = rod.getVertex(vidx);

        Vec3<Real> v0_otherMathLibrary(v0[0], v0[1], v0[2]);

	Scalar d = m_levelSet->getLevelSetValue( v0_otherMathLibrary );

        Vec3<Real> dgrad_otherMathLibrary;
	m_levelSet->getGradient(v0_otherMathLibrary, dgrad_otherMathLibrary);

	Vec3d dgrad( dgrad_otherMathLibrary[0], dgrad_otherMathLibrary[1], dgrad_otherMathLibrary[2] );

        if (d >= 0) continue;

	Scalar curr_energy = 0.5 * stiffness * d*d;

        Vec3d curr_force = stiffness * d * dgrad; 

        for (int i = 0; i < 3; ++i)
        {
            force[rod.vertIdx(vidx, i)] += curr_force[i];
        }

        TraceStream(g_log, "RodLevelSetForce::computeForceEnergy") << "vidx = " << vidx << " d = " << d << " grad d = " << dgrad << " energy = " << curr_energy << " force = " << curr_force << "\n";

    }
}
