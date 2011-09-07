// RodLevelSetForce.cc
//

#include "RodLevelSetForce.hh"
#ifdef WETA
#include "../../Math/Math.hh"
#else
#include "BASim/src/Math/Math.hh"
#endif
#include "BASim/src/Math/MatrixBase.hh"

#include "../../Util/TextLog.hh"

#include "../../Collisions/LevelSet.hh"


using namespace BASim;

RodLevelSetForce::RodLevelSetForce( LevelSet* levelSet ) 
{
    m_name = "level set penalty";

    m_stiffness = 1;

    m_thickness = 0.1;

    m_subsamples = 0;

    m_levelSet = levelSet;
    assert( m_levelSet );
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

    assert( m_levelSet->isInitialized() );

    for (int vidx = 1; vidx < (int) rod.nv(); vidx++)
    {

        Vec3d v0 = rod.getVertex(vidx);

        Scalar alpha = 0.0;
        computeLocalForceDX( baseindex, rod, v0, vidx, alpha, scale, J );

        if ( vidx == rod.nv() - 1 )
        	continue;

        Vec3d sample_v = v0;
        Vec3d v1 = rod.getVertex(vidx+1);

        for ( int i = 0; i < m_subsamples; i++ )
        {
        	sample_v += ( v1 - v0 ) / ( m_subsamples + 1 );
        	alpha += 1.0 / ( m_subsamples + 1 );
        	computeLocalForceDX( baseindex, rod, sample_v, vidx, alpha, scale, J );
        }
    }
}

void RodLevelSetForce::computeLocalForceDX(int baseindex, const ElasticRod& rod, Vec3d& v0, int vidx,
		const Scalar alpha, Scalar scale, MatrixBase& J) const
{
	IntArray indices(3);
	MatXd localJ(3, 3);

	Vec3<Real> v0_otherMathLibrary(v0[0], v0[1], v0[2]);

	Scalar d = m_levelSet->getLevelSetValue( v0_otherMathLibrary );

	if (d >= m_thickness) return;

	Vec3<Real> dgrad_otherMathLibrary;
	m_levelSet->getGradient(v0_otherMathLibrary, dgrad_otherMathLibrary);

	Vec3d dgrad( dgrad_otherMathLibrary[0], dgrad_otherMathLibrary[1], dgrad_otherMathLibrary[2] );

	localJ.setZero();

	for (int i = 0; i < 3; ++i)
	{
		for (int j = 0; j < 3; ++j)
		{
			Scalar val = - scale * m_stiffness * dgrad[i] * dgrad[j];
			localJ(i,j) = val;
		}
	}

	// distribute local jacobian between this vertex and the next
	for (int i = 0; i < 3; ++i)
	{
		indices[i] = baseindex + rod.vertIdx(vidx, i);
	}
	J.add(indices, indices, (1.0-alpha)*localJ);

	if ( alpha == 0.0 || vidx == rod.nv() - 1 ) return;

	for (int i = 0; i < 3; ++i)
	{
		indices[i] = baseindex + rod.vertIdx(vidx + 1, i);
	}
	J.add(indices, indices, alpha*localJ);

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
    assert( m_levelSet->isInitialized() );

    // NOTE (sainsley) : we shouldn't check the fixed vertex as it may already be in the scalp
    for ( int vidx = 1; vidx < (int) rod.nv(); vidx++ )
    {
        Vec3d v0 = rod.getVertex(vidx);

        Scalar alpha = 0.0;
        computeLocalForceEnergy( rod, v0, vidx, alpha, force, energy );

        if ( vidx == rod.nv() - 1 )
        	continue;

        Vec3d sample_v = v0;
        Vec3d v1 = rod.getVertex(vidx+1);

        for ( int i = 0; i < m_subsamples; i++ )
        {
        	sample_v += ( v1 - v0 ) / ( m_subsamples + 1 );
        	alpha += 1.0 / ( m_subsamples + 1 );
        	computeLocalForceEnergy( rod, sample_v, vidx, alpha, force, energy );
        }

        // TODO (sainsley) : add subchecks here
        /*Vec3<Real> v0_otherMathLibrary(v0[0], v0[1], v0[2]);

        Scalar d = m_levelSet->getLevelSetValue( v0_otherMathLibrary );

        if (d >= m_thickness) continue;

        Vec3<Real> dgrad_otherMathLibrary;
        m_levelSet->getGradient(v0_otherMathLibrary, dgrad_otherMathLibrary);

        Vec3d dgrad( dgrad_otherMathLibrary[0], dgrad_otherMathLibrary[1], dgrad_otherMathLibrary[2] );

        Scalar curr_energy;
        Vec3d curr_force;

        curr_energy = 0.5 * m_stiffness * (m_thickness-d)*(m_thickness-d);
        curr_force = m_stiffness * (m_thickness-d) * dgrad;

        for (int i = 0; i < 3; ++i)
        {
            force[rod.vertIdx(vidx, i)] += curr_force[i];
        }

        energy += curr_energy;

        TraceStream(g_log, "RodLevelSetForce::computeForceEnergy") << "vidx = " << vidx << " stiff = " << m_stiffness << " d = " << d << " grad d = " << dgrad << " energy = " << curr_energy << " force = " << curr_force << "\n";
        */
    }
}

void RodLevelSetForce::computeLocalForceEnergy(const ElasticRod& rod, Vec3d& v0, int vidx, const Scalar alpha, VecXd& force, Scalar& energy) const
{
	Vec3<Real> v0_otherMathLibrary(v0[0], v0[1], v0[2]);

	Scalar d = m_levelSet->getLevelSetValue( v0_otherMathLibrary );

	if (d >= m_thickness) return;

	Vec3<Real> dgrad_otherMathLibrary;
	m_levelSet->getGradient(v0_otherMathLibrary, dgrad_otherMathLibrary);

	Vec3d dgrad( dgrad_otherMathLibrary[0], dgrad_otherMathLibrary[1], dgrad_otherMathLibrary[2] );

	Scalar curr_energy;
	Vec3d curr_force;

	curr_energy = 0.5 * m_stiffness * (m_thickness-d) * (m_thickness-d);
	curr_force = m_stiffness * (m_thickness-d) * dgrad;

	for (int i = 0; i < 3; ++i)
	{
		// distribute force along edge
		force[rod.vertIdx(vidx, i)] += (1.0-alpha)*curr_force[i];
		if ( vidx < rod.nv() - 1 )
		{
			force[rod.vertIdx(vidx+1, i)] += alpha*curr_force[i];
		}
	}

	energy += curr_energy;

	TraceStream(g_log, "RodLevelSetForce::computeForceEnergy") << "vidx = " << vidx << " stiff = " << m_stiffness << " d = " << d << " grad d = " << dgrad << " energy = " << curr_energy << " force = " << curr_force << "\n";

}
