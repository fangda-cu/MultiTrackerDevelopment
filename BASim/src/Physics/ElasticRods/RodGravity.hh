/**
 * \file RodGravity.hh
 *
 * \author miklos@cs.columbia.edu
 * \date 09/04/2009
 */

#ifndef RODGRAVITY_HH
#define RODGRAVITY_HH

namespace BASim
{

/** External force that applies gravity to a rod. */
class RodGravity: public RodExternalConservativeForce
{
public:

    /**
     * Constructor for creating the gravity force.
     *
     * \param[in] gravity The acceleration due to gravity to be applied to rods.
     */

    explicit RodGravity(const Vec3d& gravity) :
        m_gravity(gravity)
    {
        m_name = "gravity";
    }

    /**
     * Computes the force due to gravity for the given rod and adds it
     * to the given vector.
     *
     * \param[in] rod The rod to compute the gravity force for.
     * \param[out] force Vector storing the forces on the rod.
     */

    Scalar computeEnergy(const ElasticRod& rod) const
    {
        Scalar energy = 0;

        for (int i = 0; i < rod.nv(); ++i)
        {
            energy += rod.getVertexMass(i) * m_gravity.dot( rod.getVertex(i));
        }

        return energy;
    }

    void computeForce(const ElasticRod& rod, VecXd& force) const
    {
        VecXd force1 = force;

        for (int i = 0; i < rod.nv(); ++i)
        {
            Vec3d f = rod.getVertexMass(i) * m_gravity;
            for (int coord = 0; coord < 3; ++coord)
            {
                force(rod.vertIdx(i, coord)) += f(coord);
            }
        }

        //    std::cout << "gravity\n" << force - force1 <<"\n";

    }


    void computeForceEnergy(const ElasticRod& rod, VecXd& force, Scalar& energy) const
    {
        VecXd force1 = force;

        for (int i = 0; i < rod.nv(); ++i)
        {
            energy += rod.getVertexMass(i) * m_gravity.dot( rod.getVertex(i));

            Vec3d f = rod.getVertexMass(i) * m_gravity;
            for (int coord = 0; coord < 3; ++coord)
            {
                force(rod.vertIdx(i, coord)) += f(coord);
            }
        }

        //    std::cout << "gravity\n" << force - force1 <<"\n";

    }


    void computeForceDX(int baseidx, const ElasticRod& rod, Scalar scale, MatrixBase& J) const
    {
    }

    void computeForceDV(int baseidx, const ElasticRod& rod, Scalar scale, MatrixBase& J) const
    {
    }

protected:

    Vec3d m_gravity;
};

} // namespace BASim

#endif // RODGRAVITY_HH
