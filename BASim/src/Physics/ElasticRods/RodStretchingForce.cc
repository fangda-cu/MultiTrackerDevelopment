/**
 * \file RodStretchingForce.cc
 *
 * \author miklos@cs.columbia.edu
 * \date 09/01/2009
 */

#include "RodStretchingForce.hh"
#include "../../Math/Math.hh"
#include "BASIM/src/Math/MatrixBase.hh"

#ifdef TEST_ROD_STRETCHING
#include "BASim/src/Physics/ElasticRods/Tests/RodStretchingTest.hh"
#endif // TEST_ROD_STRETCHING
using namespace std;

namespace BASim
{

RodStretchingForce::RodStretchingForce(ElasticRod& rod, bool vscs, bool runinit) :
    RodForceT<EdgeStencil> (rod, "RodStretchingForce")
{
    if (runinit)
    {
        if (!vscs)
        {
            m_rod.add_property(m_ks, "stretching stiffness");
            m_rod.add_property(m_refLength, "stretching ref length");
        }
        else
        {
            m_rod.add_property(m_ks, "viscous stretching stiffness");
            m_rod.add_property(m_refLength, "viscous stretching ref length");
        }

        setViscous(vscs);

        updateStiffness();
        updateUndeformedStrain();
    }
    else
    {
        m_viscous = vscs;
    }
}

void RodStretchingForce::updateUndeformedConfiguration(std::vector<Scalar>& vals)
{
    /*std::cout << "stretch origin\n";

     {
     iterator end = m_stencil.end();
     for (m_stencil = m_stencil.begin(); m_stencil != end; ++m_stencil) {
     edge_handle& eh = m_stencil.handle();
     std::cout << getRefLength(eh) << "\n";
     }
     }*/

    //std::cout << "stretch new\n";

    int i = 0;
    iterator end = m_stencil.end();
    for (m_stencil = m_stencil.begin(); m_stencil != end; ++m_stencil, ++i)
    {
        if (i == 0)
            continue;
        edge_handle& eh = m_stencil.handle();
        setRefLength(eh, vals[i - 1]);
        //    std::cout << getRefLength(eh) << "\n";
    }
}

void RodStretchingForce::setReferenceLengths(std::vector<Scalar>& vals)
{
    updateUndeformedConfiguration(vals);
}

void RodStretchingForce::gatherDofs(SpringDofStruct& dofs, const edge_handle& eh)
{
    dofs.x[0] = m_rod.getFromVertex(eh);
    dofs.x[1] = m_rod.getToVertex(eh);
    dofs.edge = m_rod.getEdge(eh);
    dofs.tangent = m_rod.getTangent(eh);
    dofs.currLength = m_rod.getEdgeLength(eh);
    dofs.restLength = getRefLength(eh);
    dofs.stiffness = getKs(eh);

    // std::cout << "RodStretchingForce: strain = " << dofs.currLength/dofs.restLength << " l = " << dofs.currLength << " lrest = " << dofs.restLength << " x0 = " << dofs.x[0] << " x1 = " << dofs.x[1] << std::endl;
}

Scalar RodStretchingForce::globalEnergy()
{
    Scalar energy = 0;

    iterator end = m_stencil.end();
    for (m_stencil = m_stencil.begin(); m_stencil != end; ++m_stencil)
    {
        edge_handle& eh = m_stencil.handle();
        Scalar localEnergy = elementEnergy(eh);
        energy += localEnergy;

#ifdef TEST_ROD_STRETCHING
        testEnergy(localEnergy, eh);
#endif // TEST_ROD_STRETCHING
    }
    return energy;
}

Scalar RodStretchingForce::elementEnergy(const edge_handle& eh)
{
    Scalar ks = getKs(eh);
    if (ks == 0.0)
        return 0;

    Scalar refLength = getRefLength(eh);
    Scalar len = m_rod.getEdgeLength(eh);

    return ks / 2.0 * square(len / refLength - 1.0) * refLength;
}

void RodStretchingForce::globalForce(VecXd& force)
{
    ElementForce localForce;

    const iterator end = m_stencil.end();
    for (m_stencil = m_stencil.begin(); m_stencil != end; ++m_stencil)
    {
        elementForce(localForce, m_stencil.handle());

        const int fi = m_stencil.firstIndex();
        for (int i = 0; i < 3; ++i)
            force(fi + i) += localForce(i);
        for (int i = 3; i < 6; ++i)
            force(fi + 1 + i) += localForce(i);

#ifdef TEST_ROD_STRETCHING
        testForce(localForce, eh);
#endif // TEST_ROD_STRETCHING
    }

#ifdef TEST_ROD_STRETCHING
    globalEnergy();
#endif // TEST_ROD_STRETCHING
    //std::cout << (viscous() ? "VISCOUS " : "") << "STRETCHING FORCE\n";
    //std::cout << "norm = " << (force - force1).norm() << " force = " << force-force1 << "\n\n";
}

void RodStretchingForce::elementForce(ElementForce& force, const SpringDofStruct& dofs)
{
    Vec3d f = dofs.stiffness * (dofs.currLength / dofs.restLength - 1.0) * dofs.tangent;
    force.segment<3> (0) = f;
    force.segment<3> (3) = -f;
}

void RodStretchingForce::elementForce(ElementForce& force, const edge_handle& eh)
{
    const Scalar ks = getKs(eh);
    if (ks == 0.0)
        return;

    const Scalar refLength = getRefLength(eh);
    const Scalar len = m_rod.getEdgeLength(eh);
    const Vec3d& tangent = m_rod.getTangent(eh);

    const Vec3d f = ks * (len / refLength - 1.0) * tangent;
    force.segment<3> (0) = f;
    force.segment<3> (3) = -f;
}

void RodStretchingForce::globalJacobian(int baseidx, Scalar scale, MatrixBase& Jacobian)
{
    ElementJacobian localJ;
    const iterator end = m_stencil.end();
    for (m_stencil = m_stencil.begin(); m_stencil != end; ++m_stencil)
    {
        edge_handle& eh = m_stencil.handle();
        localJ.setZero();
        elementJacobian(localJ, eh);

#ifdef TEST_ROD_STRETCHING
        testJacobian(localJ, eh);
#endif // TEST_ROD_STRETCHING
        localJ *= scale;
        Jacobian.edgeStencilAdd(m_stencil.firstIndex() + baseidx, localJ);
    }

    assert(isSymmetric(Jacobian));
}

void RodStretchingForce::globalForceEnergy(VecXd& force, Scalar& energy)
{
    IndexArray indices;
    ElementForce localForce;

    iterator end = m_stencil.end();

    for (m_stencil = m_stencil.begin(); m_stencil != end; ++m_stencil)
    {
        edge_handle& eh = m_stencil.handle();
        elementForceEnergy(localForce, energy, eh);
        m_stencil.indices(indices);
        for (int i = 0; i < indices.size(); ++i)
            force(indices(i)) += localForce(i);
    }

}

void RodStretchingForce::elementForceEnergy(ElementForce& force, Scalar& energy, const edge_handle& eh)
{
    // NOTE(sainsley): there is no efficiency gain in this method. There is overlap in the
    // elementForce(... edge_handle) method and the elementEnergy method, but the former doesn't seem
    // to be in use.

    SpringDofStruct dofs;
    gatherDofs(dofs, eh);

    // set force
    Vec3d f = dofs.stiffness * (dofs.currLength / dofs.restLength - 1.0) * dofs.tangent;
    force.segment<3> (0) = f;
    force.segment<3> (3) = -f;

    Scalar ks = getKs(eh);
    if (ks == 0.0)
        return;

    // accumulate energy
    Scalar refLength = getRefLength(eh);
    Scalar len = m_rod.getEdgeLength(eh);
    energy += ks / 2.0 * square(len / refLength - 1.0) * refLength;
}

void RodStretchingForce::globalJacobianForceEnergy(int baseidx, Scalar scale, MatrixBase& Jacobian, VecXd& force,
        Scalar& energy)
{
    IndexArray indices;
    ElementJacobian localJ;
    ElementForce localForce;
    MatXd adder;

    iterator end = m_stencil.end();
    for (m_stencil = m_stencil.begin(); m_stencil != end; ++m_stencil)
    {
        edge_handle& eh = m_stencil.handle();
        localJ.setZero();
        elementJacobianForceEnergy(localJ, localForce, energy, eh);
        adder = localJ;
        adder *= scale;
        m_stencil.indices(indices);
        for (int i = 0; i < (int) indices.size(); ++i)
        {
            force(indices(i)) += localForce(i);
            indices(i) += baseidx;
        }
        Jacobian.add(indices, indices, adder);
    }
}

void RodStretchingForce::elementJacobianForceEnergy(ElementJacobian& Jacobian, ElementForce& force, Scalar& energy,
        const edge_handle& eh)
{
    SpringDofStruct dofs;
    gatherDofs(dofs, eh);

    // set force
    Vec3d f = dofs.stiffness * (dofs.currLength / dofs.restLength - 1.0) * dofs.tangent;
    force.segment<3> (0) = f;
    force.segment<3> (3) = -f;

    Scalar ks = getKs(eh);
    if (ks == 0.0)
        return;

    const Vec3d& e = m_rod.getEdge(eh);
    Scalar len = m_rod.getEdgeLength(eh);
    Scalar refLength = getRefLength(eh);
    Mat3d M = ks * ((1.0 / refLength - 1.0 / len) * Mat3d::Identity() + 1.0 / len * outerProd(e, e) / square(len));

    // set jacobian
    for (int i = 0; i < 3; ++i)
    {
        for (int j = 0; j < 3; ++j)
        {
            Jacobian(i, j) = Jacobian(3 + i, 3 + j) = -M(i, j);
            Jacobian(3 + i, j) = Jacobian(i, 3 + j) = M(i, j);
        }
    }

    // accumulate energy
    energy += ks / 2.0 * square(len / refLength - 1.0) * refLength;
}

void RodStretchingForce::elementJacobian(ElementJacobian& Jacobian, const edge_handle& eh)
{
    Scalar ks = getKs(eh);
    if (ks == 0.0)
        return;

    const Vec3d& e = m_rod.getEdge(eh);
    Scalar len = m_rod.getEdgeLength(eh);
    Scalar refLength = getRefLength(eh);
    Mat3d M = ks * ((1.0 / refLength - 1.0 / len) * Mat3d::Identity() + 1.0 / len * outerProd(e, e) / square(len));

    for (int i = 0; i < 3; ++i)
    {
        for (int j = 0; j < 3; ++j)
        {
            Jacobian(i, j) = Jacobian(3 + i, 3 + j) = -M(i, j);
            Jacobian(3 + i, j) = Jacobian(i, 3 + j) = M(i, j);
        }
    }

    assert(isSymmetric(Jacobian));
}

void RodStretchingForce::globalReverseJacobian(MatrixBase& J)
{
    if (viscous())
        return; // don't need to include viscous force - zero

    IndexArray indices;

    iterator end = m_stencil.end();
    m_stencil = m_stencil.begin();

    ++m_stencil; // skip first edge

    unsigned int eid = 1;
    for (; m_stencil != end; ++m_stencil, ++eid)
    {
        edge_handle& eh = m_stencil.handle();

        Scalar ks = getKs(eh);
        if (ks == 0.0)
            return;

        const Vec3d& e = m_rod.getEdge(eh);
        //Scalar len = m_rod.getEdgeLength(eh);
        Scalar refLength = getRefLength(eh);

        Vec3d fdx = -ks / (refLength * refLength) * e;

        m_stencil.indices(indices);

        //    std::cout << "str id  " << indices << "\n" << fdx << "\n";

        if (eid > 1)
        {
            for (int i = 0; i < 3; i++)
            {
                J.add(indices[i] - 7, (eid - 1) * 4 + 3, fdx(i));
            }
        }

        //    J.print();

        for (int i = 3; i < 6; i++)
        {
            //      cout << indices[i] - 7 << " " << (eid-1) * 4 + 3 << " "<< -fdx(i-3) <<"\n";
            J.add(indices[i] - 7, (eid - 1) * 4 + 3, -fdx(i - 3));
        }

        //    J.print();
    }
}

// we know index
void RodStretchingForce::updateReverseUndeformedStrain(const VecXd& e)
{
    if (viscous())
        return;

    m_stencil = m_stencil.begin();
    iterator end = m_stencil.end();

    ++m_stencil;
    unsigned int eid = 1; // eh.id() ?

    for (; m_stencil != end; ++m_stencil, ++eid)
    {
        edge_handle& eh = m_stencil.handle();
        setRefLength(eh, e((eid - 1) * 4 + 3));
    }
}

const Scalar& RodStretchingForce::getKs(const edge_handle& eh) const
{
    return m_rod.property(m_ks)[eh];
}

void RodStretchingForce::setKs(const edge_handle& eh, const Scalar& ks)
{
    m_rod.property(m_ks)[eh] = ks;
}

const Scalar& RodStretchingForce::getRefLength(const edge_handle& eh) const
{
    return m_rod.property(m_refLength)[eh];
}

void RodStretchingForce::setRefLength(const edge_handle& eh, const Scalar& length)
{
    m_rod.property(m_refLength)[eh] = length;
}

void RodStretchingForce::updateStiffness()
{
    Scalar E = m_rod.getYoungsModulus();
    if (viscous())
    {
        E = 3 * m_rod.getViscosity() / m_rod.getTimeStep();
    }

    iterator end = m_stencil.end();
    for (m_stencil = m_stencil.begin(); m_stencil != end; ++m_stencil)
    {
        edge_handle& eh = m_stencil.handle();
        Scalar a = m_rod.radiusA(eh);
        Scalar b = m_rod.radiusB(eh);
        setKs(eh, E * M_PI * a * b);

        //std::cout << E * M_PI * a * b << "\n";

    }
}

void RodStretchingForce::updateUndeformedStrain()
{
    iterator end = m_stencil.end();
    for (m_stencil = m_stencil.begin(); m_stencil != end; ++m_stencil)
    {
        edge_handle& eh = m_stencil.handle();
        setRefLength(eh, m_rod.getEdgeLength(eh));
    }
}

#ifdef TEST_ROD_STRETCHING

void RodStretchingForce::testEnergy(const Scalar& energy,
        const edge_handle& eh) const
{
    Scalar ks = getKs(eh);
    Scalar refLength = getRefLength(eh);
    Scalar len = m_rod.getEdgeLength(eh);
    const Vec3d& x0 = m_rod.getFromVertex(eh);
    const Vec3d& x1 = m_rod.getToVertex(eh);
    Scalar mathEnergy;
    rodStretchingEnergyTest(mathEnergy, energy, x0, x1, ks, refLength);
}

void RodStretchingForce::testForce(const ElementForce& force,
        const edge_handle& eh) const
{
    Scalar ks = getKs(eh);
    Scalar refLength = getRefLength(eh);
    Scalar len = m_rod.getEdgeLength(eh);
    const Vec3d& tangent = m_rod.getTangent(eh);
    const Vec3d& x0 = m_rod.getFromVertex(eh);
    const Vec3d& x1 = m_rod.getToVertex(eh);
    ElementForce mathForce;
    rodStretchingForceTest(mathForce, force, x0, x1, ks, refLength);
}

void RodStretchingForce::testJacobian(const ElementJacobian& Jacobian,
        const edge_handle& eh) const
{
    const Vec3d& x0 = m_rod.getFromVertex(eh);
    const Vec3d& x1 = m_rod.getToVertex(eh);
    Scalar ks = getKs(eh);
    Scalar refLength = getRefLength(eh);
    ElementJacobian mathJ;
    rodStretchingJacobianTest(mathJ, Jacobian, x0, x1, ks, refLength);
}

#endif // TEST_ROD_STRETCHING
} // namespace BASim
