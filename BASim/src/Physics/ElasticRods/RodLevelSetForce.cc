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

    m_stiffness = 1;

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
    MatXd localJ(3, 3);
    IntArray indices(3);

    assert( m_levelSet->isInitialized() );

    for (int vidx = 1; vidx < (int) rod.nv(); vidx++)
    {
        Vec3d v0 = rod.getVertex(vidx);

        Vec3<Real> v0_otherMathLibrary(v0[0], v0[1], v0[2]);

        Scalar d = m_levelSet->getLevelSetValue( v0_otherMathLibrary );
        Scalar tol = 0.1;
	    // TODO: grab this distance tolerance from the implicit settings in maya
        if (d >= tol) continue;
        //if (d >= 0) continue;

        Vec3<Real> dgrad_otherMathLibrary;
        m_levelSet->getGradient(v0_otherMathLibrary, dgrad_otherMathLibrary);

        Vec3d dgrad( dgrad_otherMathLibrary[0], dgrad_otherMathLibrary[1], dgrad_otherMathLibrary[2] );

        Scalar curr_energy = 0.5 * m_stiffness * (tol-d)*(tol-d); //0.5 * m_stiffness * d*d;

        localJ.setZero();

        for (int i = 0; i < 3; ++i)
        {
        	for (int j = 0; j < 3; ++j)
        	{
        		Scalar val = - scale * m_stiffness * dgrad[i] * dgrad[j];
        		localJ(i,j) = val;
        	}
        }

        for (int i = 0; i < 3; ++i)
        {
            indices[i] = baseindex + rod.vertIdx(vidx, i);
        }
        J.add(indices, indices, localJ);
    }
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
    for (int vidx = 1; vidx < (int) rod.nv(); vidx++)
    {
        Vec3d v0 = rod.getVertex(vidx);

        Vec3<Real> v0_otherMathLibrary(v0[0], v0[1], v0[2]);

        Scalar d = m_levelSet->getLevelSetValue( v0_otherMathLibrary );

        Scalar tol = 0.1;
        // TODO: grab this distance tolerance from the implicit settings in maya
        if (d >= tol) continue;

        Vec3<Real> dgrad_otherMathLibrary;
        m_levelSet->getGradient(v0_otherMathLibrary, dgrad_otherMathLibrary);

        Vec3d dgrad( dgrad_otherMathLibrary[0], dgrad_otherMathLibrary[1], dgrad_otherMathLibrary[2] );

        //Scalar d_orig = d;
        // todo factor in tolerance here for penalty spring compression
        Scalar curr_energy;
        Vec3d curr_force;
        //if ( d >= 0 )
        //{
        	//d = fabs(d);
        curr_energy = 0.5 * m_stiffness * (tol-d)*(tol-d);
        curr_force = m_stiffness * (tol-d) * dgrad;
        //}
        //else
        //{
        	 //d = fabs(d);
        	// curr_energy = 0.5 * m_stiffness*m_stiffness * (tol-d)*(tol-d);
        	 //curr_force = m_stiffness*m_stiffness * (tol-d) * dgrad;
        //}

        // todo factor in tolerance here for penalty spring compression
        //Vec3d curr_force = m_stiffness * fabs(d-tol) * dgrad;

        for (int i = 0; i < 3; ++i)
        {
            force[rod.vertIdx(vidx, i)] += curr_force[i];
        }

        energy += curr_energy;

        TraceStream(g_log, "RodLevelSetForce::computeForceEnergy") << "vidx = " << vidx << " stiff = " << m_stiffness << " d = " << d << " grad d = " << dgrad << " energy = " << curr_energy << " force = " << curr_force << "\n";
    }
}
