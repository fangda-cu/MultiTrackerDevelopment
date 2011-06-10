/*
 * BARodStepper_CollisionResponse.cc
 *
 *  Created on: 19/05/2011
 *
 */

#include <typeinfo>
#include "BARodStepper.hh"
#include "../../Threads/MultithreadedStepper.hh"
#include "../../Core/Timer.hh"
#include "../../Collisions/Collision.hh"

#include <iostream.h>
#ifdef HAVE_OPENMP
#include <omp.h>
#endif

namespace BASim
{

/**
 * Inelastic response
 */
void BARodStepper::exertVertexImpulse(const Vec3d& I, const double& m, const int& idx, VecXd& v)
{
    assert(m > 0.0);
    assert(idx >= 0);
    assert(idx < getNumVerts());
    assert(v.size() == getNumDof());

    v.segment<3> (3 * idx) += I / m;
}

void BARodStepper::exertEdgeImpulse(const Vec3d& I, const double& m0, const double& m1, const double& alpha, const int& idx0,
        const int& idx1, VecXd& v)
{
    assert(m0 > 0.0);
    assert(m1 > 0.0);
    assert(alpha >= 0.0);
    assert(alpha <= 1.0);
    assert(idx0 >= 0);
    assert(idx0 < getNumVerts());
    assert(idx1 >= 0);
    assert(idx1 < getNumVerts());
    assert(v.size() == getNumDof());

    v.segment<3> (3 * idx0) += ((1 - alpha) / m0) * I;
    v.segment<3> (3 * idx1) += (alpha / m1) * I;
}

void BARodStepper::exertFaceImpulse(const Vec3d& I, const double& m0, const double& m1, const double& m2, const double& u,
        const double& v, const double& w, const int& idx0, const int& idx1, const int& idx2, VecXd& vel)
{
    assert(m0 > 0.0);
    assert(m1 > 0.0);
    assert(m2 > 0.0);
    //assert( u >= 0.0 ); assert( u <= 1.0 );
    //assert( v >= 0.0 ); assert( v <= 1.0 );
    //assert( w >= 0.0 ); assert( w <= 1.0 );
    assert(approxEq(u + v + w, 1.0));
    assert(idx0 >= 0);
    assert(idx0 < getNumVerts());
    assert(idx1 >= 0);
    assert(idx1 < getNumVerts());
    assert(idx2 >= 0);
    assert(idx2 < getNumVerts());
    assert(vel.size() == getNumDof());

    vel.segment<3> (3 * idx0) += u * I / m0;
    vel.segment<3> (3 * idx1) += v * I / m1;
    vel.segment<3> (3 * idx2) += w * I / m2;
}

void BARodStepper::exertInelasticImpulse(EdgeEdgeCTCollision& eecol)
{
    assert(eecol.e0_v0 >= 0);
    assert(eecol.e0_v0 < getNumVerts());
    assert(eecol.e0_v1 >= 0);
    assert(eecol.e0_v1 < getNumVerts());
    assert(eecol.e1_v0 >= 0);
    assert(eecol.e1_v0 < getNumVerts());
    assert(eecol.e1_v1 >= 0);
    assert(eecol.e1_v1 < getNumVerts());

    // Compute the relative velocity of the edges at the collision point
    //   Vec3d relvel = m_geodata.computeRelativeVelocity(clssn.e0_v0, clssn.e0_v1, clssn.e1_v0, clssn.e1_v1, clssn.s, clssn.t);
    //   double magrelvel = relvel.dot(clssn.n);
    //    double magrelvel = eecol.computeRelativeVelocity(m_geodata);

    // TraceStream(g_log, "") << "BARodStepper:exertInelasticImpulse: pre-impulse e-e relative velocity = "
    //         << eecol.computeRelativeVelocity() << '\n';

    if (eecol.GetCachedRelativeVelocity() >= 0.0)
    {
        DebugStream(g_log, "")
                << "\033[31;1mWARNING IN BRIDSON STEPPER:\033[m Relative velocity computed \
                  incorrectly before applying edge-edge inelastic impulse (bug in normal \
                  computation?). Magnitude of relative velocity: "
                << eecol.computeRelativeVelocity() << ". Rod: " << getContainingRod(eecol.e0_v0) << '\n';
    }
    assert(eecol.GetCachedRelativeVelocity() < 0.0);

    // Add some extra "kick" to relative velocity to account for FPA errors
    //eecol.ApplyRelativeVelocityKick();
    Vec3d I = eecol.computeInelasticImpulse();
    //computeEdgeEdgeInelasticImpulse(m_masses[eecol.e0_v0], m_masses[eecol.e0_v1], m_masses[eecol.e1_v0],
    //m_masses[eecol.e1_v1], eecol.s, eecol.t, eecol.GetRelativeVelocity(), eecol.n);

    // TraceStream(g_log, "") << "BARodStepper::exertInelasticEdgeEdgeImpulse: (" << eecol.e0_v0 << "-" << eecol.e0_v1 << ", " << eecol.e1_v0
    //         << "-" << eecol.e1_v1 << ") s=" << eecol.s << " t=" << eecol.t << " I=" << I << '\n';

    exertEdgeImpulse(-I, m_masses[eecol.e0_v0], m_masses[eecol.e0_v1], eecol.s, eecol.e0_v0, eecol.e0_v1, m_vnphalf);
    exertEdgeImpulse(I, m_masses[eecol.e1_v0], m_masses[eecol.e1_v1], eecol.t, eecol.e1_v0, eecol.e1_v1, m_vnphalf);

    // TraceStream(g_log, "") << "BARodStepper:exertInelasticImpulse: post-impulse e-e relative velocity = "
    //         << eecol.computeRelativeVelocity() << '\n';

    assert(eecol.computeRelativeVelocity() >= 0);
}

void BARodStepper::exertInelasticImpulse(VertexFaceCTCollision& vfcol)
{
    assert(vfcol.v0 >= 0);
    assert(vfcol.v0 < getNumVerts());
    assert(vfcol.f0 >= 0);
    assert(vfcol.f0 < getNumVerts());
    assert(vfcol.f1 >= 0);
    assert(vfcol.f1 < getNumVerts());
    assert(vfcol.f2 >= 0);
    assert(vfcol.f2 < getNumVerts());

    // Compute the relative velocity of the edges at the collision point
    //   double magrelvel = vfcol.computeRelativeVelocity(m_geodata);
    assert(vfcol.GetCachedRelativeVelocity() < 0.0);

    // TraceStream(g_log, "") << "BARodStepper:exertInelasticImpulse: pre-impulse v-f relative velocity = "
    //         << vfcol.GetCachedRelativeVelocity() << '\n';

    // Add some extra "kick" to relative velocity to account for FPA errors
    //vfcol.ApplyRelativeVelocityKick();
    Vec3d I = vfcol.computeInelasticImpulse();
    // computeVertexFaceInelasticImpulse(m_masses[vfcol.v0], m_masses[vfcol.f0], m_masses[vfcol.f1], m_masses[vfcol.f2],
    // vfcol.u, vfcol.v, vfcol.w, eecol.GetRelativeVelocity(), vfcol.n);

    // TraceStream(g_log, "") << "BARodStepper::exertInelasticImpulse: (" << vfcol.v0 << ", " << vfcol.f0 << "-" << vfcol.f1 << "-"
    //         << vfcol.f2 << ") (u,v,w)=(" << vfcol.u << "," << vfcol.v << "," << vfcol.w << ") I=" << I << '\n';

    exertFaceImpulse(-I, m_masses[vfcol.f0], m_masses[vfcol.f1], m_masses[vfcol.f2], vfcol.u, vfcol.v, vfcol.w, vfcol.f0,
            vfcol.f1, vfcol.f2, m_vnphalf);
    exertVertexImpulse(I, m_masses[vfcol.v0], vfcol.v0, m_vnphalf);

    // TraceStream(g_log, "") << "BARodStepper:exertInelasticImpulse: post-impulse v-f relative velocity = "
    //         << vfcol.computeRelativeVelocity() << '\n';

    assert(vfcol.computeRelativeVelocity() >= 0);
}

/**
 * Compliant inelastic response
 */
bool BARodStepper::executeIterativeInelasticImpulseResponse(std::vector<bool>& failed_collisions_rods,
        std::vector<bool>& stretching_rods)
{
    bool all_rods_collisions_ok = true;

    // Check whether the solver left some rods stretched
    checkLengths(stretching_rods);

    // Detect continuous time collisions
    std::list<Collision*> collisions_list;
    TraceStream(g_log, "") << "Detecting collisions...\n";
    m_collision_detector->getCollisions(collisions_list, ContinuousTime);
    TraceStream(g_log, "") << "Initial potential collisions: " << m_collision_detector->m_potential_collisions << "\n";

    // Iteratively apply inelastic impulses
    for (int itr = 0; !collisions_list.empty() && itr < m_perf_param.m_collision.m_max_iterations; ++itr)
    {
        TraceStream(g_log, "") << "CTcollision response iteration " << itr << '\n';
        TraceStream(g_log, "") << "Detected " << collisions_list.size() << " continuous time collisions (potential: "
                << m_collision_detector->m_potential_collisions << ")\n";

        // Just sort the collision times to maintain some rough sense of causality
        collisions_list.sort(CompareTimes);

        if (m_perf_param.m_skipRodRodCollisions)
        {
            // Keep only one collision per rod
            std::set<int> already_collided_rods;
            for (std::list<Collision*>::iterator col_it = collisions_list.begin(); col_it != collisions_list.end(); col_it++)
            {
                CTCollision* collision = dynamic_cast<CTCollision*> (*col_it);
                if (collision)
                {
                    // Only keep collisions involving a rod we see for the first time
                    int colliding_rod = getContainingRod(collision->GetRodVertex());
                    if (already_collided_rods.find(colliding_rod) == already_collided_rods.end())
                        already_collided_rods.insert(colliding_rod);
                    else
                    {
                        delete collision;
                        collisions_list.erase(col_it--);
                    }
                }
            }
            TraceStream(g_log, "") << "of which " << collisions_list.size() << " are on different rods\n";

            // Now apply response to the remaining collisions
            for (std::list<Collision*>::iterator col_it = collisions_list.begin(); col_it != collisions_list.end(); col_it++)
            {
                const CTCollision* collision = dynamic_cast<CTCollision*> (*col_it);
                assert(collision != NULL);

                // So, which rod was it? Here we are assuming that this is not a rod-rod collision
                const int collidingRodIdx = getContainingRod(collision->GetRodVertex());
                ElasticRod* const collidingRod = m_rods[collidingRodIdx];
                RodSelectionType oneRodList;
                oneRodList.push_back(collidingRodIdx);

                // Save pre-impulse velocities
                VecXd velBackup(3 * collidingRod->nv());
                for (int i = 0; i < 3 * collidingRod->nv(); i++)
                    velBackup[i] = m_vnphalf[m_base_dof_indices[collidingRodIdx] + i];

                // Position the rod at end of time step so the compliance works on the "right" configuration
                restorePositions(m_xn + m_dt * m_vnphalf, oneRodList);
                collidingRod->updateProperties();

                exertCompliantInelasticImpulse(collision);

                if (m_perf_param.m_enable_explosion_detection)
                {
                    // Prepare the rod for explosion checking
                    restorePositions(m_xn + m_dt * m_vnphalf, oneRodList);
                    collidingRod->updateProperties();
                    collidingRod->computeForces(*m_endForces[collidingRodIdx]);

                    // Line search on the impulse size: if full impulse causes explosion, halve it and so on.
                    int splitCounter = 5;
                    double splitFactor = 1.0;
                    while (hadExplosion(collidingRodIdx))
                    {
                        splitFactor *= 0.5;
                        splitCounter--;
                        TraceStream(g_log, "") << "Downsizing impulse for rod " << collidingRodIdx << " to " << splitFactor
                                * 100.0 << "%\n";
                        // Collision is marked as failed if the original impulse was exploding
                        failed_collisions_rods[collidingRodIdx] = true;
                        // Interpolate velocity between pre-impulse (velBackup) and (resized) post-impulse (m_vnphalf)
                        for (int v = 0; v < collidingRod->nv(); ++v)
                        {
                            m_vnphalf.segment<3> (m_base_dof_indices[collidingRodIdx] + 3 * v) = 0.5 * m_vnphalf.segment<3> (
                                    m_base_dof_indices[collidingRodIdx] + 3 * v) + 0.5 * velBackup.segment<3> (3 * v);
                            m_collision_immune[m_base_vtx_indices[collidingRodIdx] + v] = true;
                        }
                        if (splitCounter == 0)
                            break;

                        // Prepare the rod again for explosion checking
                        restorePositions(m_xn + m_dt * m_vnphalf, oneRodList);
                        collidingRod->updateProperties();
                        collidingRod->computeForces(*m_endForces[collidingRodIdx]);
                    }
                }

                // Test for stretching and if so, mark the rod immune
                stretching_rods[collidingRodIdx] = !checkLength(collidingRodIdx);
                if (stretching_rods[collidingRodIdx])
                {
                    // Declare the rod collision-immune for the rest of the time step
                    for (int j = 0; j < collidingRod->nv(); ++j)
                        m_collision_immune[m_base_vtx_indices[collidingRodIdx] + j] = true;
                }

                delete collision;
                collisions_list.erase(col_it--);
            }

            // Detect remaining collisions (including at the end of the last iteration, so we know what failed)
            TraceStream(g_log, "") << "Detecting collisions...\n";
            m_collision_detector->getCollisions(collisions_list, ContinuousTime, false); // No need to update the mesh bvh bounding boxes.
        }
        else
            while (!collisions_list.empty())
            {
                CTCollision* collision = dynamic_cast<CTCollision*> (collisions_list.front());
                collisions_list.pop_front();
                if (collision)
                    exertCompliantInelasticImpulse(collision);
                delete collision;
                m_collision_detector->updateContinuousTimeCollisions();
            }
    }

    // If collisions were still detected before exiting the loop
    if (!collisions_list.empty())
    {
        all_rods_collisions_ok = false;

        TraceStream(g_log, "") << "Remains " << collisions_list.size() << " unresolved collisions (potential: "
                << m_collision_detector->m_potential_collisions << ")\n";

        // Just in case we haven't emptied the collisions but exited when itr == m_num_inlstc_itrns
        for (std::list<Collision*>::iterator col = collisions_list.begin(); col != collisions_list.end(); col++)
        {
            // Let's see which collisions are remaining.
            EdgeEdgeCTCollision* eecol = dynamic_cast<EdgeEdgeCTCollision*> (*col);
            VertexFaceCTCollision* vfcol = dynamic_cast<VertexFaceCTCollision*> (*col);
            if (vfcol)
            {
                assert(isRodVertex(vfcol->v0));
                failed_collisions_rods[getContainingRod(vfcol->v0)] = true;
            }
            if (eecol)
            {
                if (isRodVertex(eecol->e0_v0) && isRodVertex(eecol->e0_v1))
                {
                    assert(getContainingRod(eecol->e0_v0) == getContainingRod(eecol->e0_v1));
                    failed_collisions_rods[getContainingRod(eecol->e0_v0)] = true;
                }
                if (isRodVertex(eecol->e1_v0) && isRodVertex(eecol->e1_v1))
                {
                    assert(getContainingRod(eecol->e1_v0) == getContainingRod(eecol->e1_v1));
                    failed_collisions_rods[getContainingRod(eecol->e1_v0)] = true;
                }
            }
            delete *col;
        }
    }

    return all_rods_collisions_ok;
}

void BARodStepper::computeCompliantLHS(MatrixBase* lhs, int rodidx)
{
    assert(lhs != NULL);
    assert(rodidx >= 0);
    assert(rodidx < m_number_of_rods);

    //double emphasizeStaticEquilibium = 0;
    //InfoStream(g_log,"") << "BARodStepper::computeCompliantLHS: emphasizing static equilibrium = " << emphasizeStaticEquilibium;

    // lhs = -h^2*dF/dx
    lhs->setZero();
    //TraceStream(g_log, "") << "WARNING: COMPLIANCE is disabled!" << '\n';
    //m_rods[rodidx]->computeJacobian(0, -m_dt * m_dt - emphasizeStaticEquilibium, *lhs);
    m_rods[rodidx]->computeJacobian(0, -m_dt * m_dt, *lhs);
    lhs->finalize();

    // lhs = M - h^2*dF/dx
    for (int i = 0; i < m_rods[rodidx]->ndof(); ++i)
        lhs->add(i, i, m_rods[rodidx]->getMass(i));
    lhs->finalize();

    assert(lhs->isApproxSymmetric(1.0e-6));

    // Set rows/cols of twist DOFs to identity
    // TODO: Handle this better
    EPropHandle<int> edgeIdx;
    m_rods[rodidx]->property_handle(edgeIdx, "edge index");
    std::vector<int> dofs_to_clear;
    for (ElasticRod::edge_iter eit = m_rods[rodidx]->edges_begin(); eit != m_rods[rodidx]->edges_end(); ++eit)
    {
        int twist_dof = m_rods[rodidx]->property(edgeIdx)[*eit];
        dofs_to_clear.push_back(twist_dof);
    }

    lhs->zeroRows(dofs_to_clear, 1.0);
    lhs->finalize();
    lhs->zeroCols(dofs_to_clear, 1.0);
    lhs->finalize();

    assert(lhs->isApproxSymmetric(1.0e-6));

    // TEMP
    //  ElasticRod::RodForces& forces = m_rods[rodidx]->getForces();
    //  ElasticRod::RodForces::iterator fIt;
    //  for (fIt = forces.begin(); fIt != forces.end(); ++fIt)
    //    //(*fIt)->globalJacobian(baseidx, scale, J);
    //    TraceStream(g_log, "") << "Force: " << (*fIt)->getName() << '\n';
    // END TEMP
}

void BARodStepper::exertCompliantInelasticImpulse(const CTCollision* cllsn)
{
    assert(cllsn->isAnalysed());
    const EdgeEdgeCTCollision* eecol = dynamic_cast<const EdgeEdgeCTCollision*> (cllsn);
    const VertexFaceCTCollision* vfcol = dynamic_cast<const VertexFaceCTCollision*> (cllsn);
    //TraceStream(g_log, "") << "BARodStepper:exertCompliantInelasticImpulse: pre-impulse e-e relative velocity = " << cllsn->GetRelativeVelocity() << '\n';

    if (eecol)
    {
        exertCompliantInelasticEdgeEdgeImpulse(*eecol);
        //exertInelasticImpulse(cllsn.getEdgeEdge());
    }
    else if (vfcol)
    {
        exertCompliantInelasticVertexFaceImpulse(*vfcol);
        //exertCompliantInelasticVertexFaceImpulse(cllsn.getVertexFace());
        //exertInelasticImpulse(cllsn.getVertexFace());
    }
    //TraceStream(g_log, "") << "BARodStepper:exertCompliantInelasticImpulse: post-impulse e-e relative velocity = " << cllsn->GetRelativeVelocity() << '\n';
}

void BARodStepper::exertCompliantInelasticVertexFaceImpulse(const VertexFaceCTCollision& vfcol)
{
    //   TraceStream(g_log, "") << "Vertex-face compliant inelastic impulse" << '\n';
    //   TraceStream(g_log, "") << vfcol << '\n';

    // TraceStream(g_log, "") << "BARodStepper::exertCompliantInelasticImpulse: (" << vfcol.v0 << ", " << vfcol.f0 << "-" << vfcol.f1 << "-"
    //         << vfcol.f2 << ") (u,v,w)=(" << vfcol.u << "," << vfcol.v << "," << vfcol.w << ")" << '\n';

    // For now, assume vertex is free and entire face is fixed
    assert(!m_geodata.isVertexFixed(vfcol.v0));
    assert(YATriangle(vfcol.f0, vfcol.f1, vfcol.f2).IsFixed(m_geodata));

    // Determine which rod the free vertex belongs to
    int rodidx = getContainingRod(vfcol.v0);

    // Ensure the free vertex belongs to a rod
    assert(rodidx >= 0);
    assert(rodidx < (int) m_number_of_rods);

    // If the rod has not solved properly, no need to compute its collision response
    if (!m_steppers[rodidx]->HasSolved())
    {// NB: this is redundant because we declared the rod collision-immune already.
        DebugStream(g_log, "") << "WARNING: attempt to do vertex-face collision with non-dependable rod\n";
        return;
    }

    // Determine which vertex of the rod the free vertex is
    assert(m_base_dof_indices[rodidx] % 3 == 0);
    int v0 = vfcol.v0 - m_base_dof_indices[rodidx] / 3;
    assert(v0 >= 0);
    assert(v0 < m_rods[rodidx]->nv());

    TraceStream(g_log, "") << "CTcollision: vertex " << v0 << " of rod " << rodidx << '\n';

    // Compute the relative velocity of the collision
    //double magrelvel = vfcol.computeRelativeVelocity(m_geodata);
    assert(vfcol.computeRelativeVelocity() < 0.0);

    // Get storage for lhs of linear system, get a solver
    LinearSystemSolver* lss = m_solver_collection.getLinearSystemSolver(m_rods[rodidx]->ndof());
    MatrixBase* lhs = lss->m_lhs;
    assert(lhs != NULL);
    LinearSolverBase* solver = lss->m_solver;
    assert(solver != NULL);

    // Compute M - h^2*dF/dx
    computeCompliantLHS(lhs, rodidx);

    // Compute the 'base' index of the vertex
    int base0 = m_rods[rodidx]->vertIdx(v0, 0);
    assert(base0 % 4 == 0);

    int ndof = m_rods[rodidx]->ndof();

    int rodbase = m_base_dof_indices[rodidx];

    // Compute the 'global' normal
    std::vector<VecXd> normal;

    normal.push_back(VecXd(ndof));
    normal.back().setZero();
    normal.back().segment<3> (base0) = vfcol.GetNormal();

    // Compute three normals per scripted vertex
    const std::vector<int>& scriptedverts = m_rods[rodidx]->getBoundaryCondition()->scriptedVertices();
    for (size_t i = 0; i < scriptedverts.size(); ++i)
    {
        for (int j = 0; j < 3; ++j)
        {
            normal.push_back(VecXd(ndof));
            normal.back().setZero();
            // Zero the velocity of this particular dof of this particular vertex
            normal.back()(m_rods[rodidx]->vertIdx(scriptedverts[i], j)) = 1.0;
        }
    }
    int numconstraints = (int) (normal.size());
    int nvdof = 3 * m_rods[rodidx]->nv();

    //TraceStream(g_log, "") << "nc: " << numconstraints << '\n';
    //for( size_t i = 0; i < n.size(); ++i ) TraceStream(g_log, "") << n[i] << '\n';

    // Determine the desired values for each constraint
    VecXd desired_values(numconstraints);
#ifndef NDEBUG
    desired_values.setConstant(std::numeric_limits<double>::signaling_NaN());
#endif
    // First constraint, the collision constraint, is 0
    desired_values(0) = 0.0;
    // The 'scripted vertex' constraints are the respective components of the velocity
    int curdof = 1;
    for (size_t i = 0; i < scriptedverts.size(); ++i)
    {
        assert(scriptedverts[i] >= 0);
        assert((int) scriptedverts[i] < m_rods[rodidx]->nv());
        desired_values(curdof++) = m_vnphalf(rodbase + 3 * scriptedverts[i] + 0);
        desired_values(curdof++) = m_vnphalf(rodbase + 3 * scriptedverts[i] + 1);
        desired_values(curdof++) = m_vnphalf(rodbase + 3 * scriptedverts[i] + 2);
    }
    assert(curdof == numconstraints);
    assert(desired_values.size() == numconstraints);

    //TraceStream(g_log, "") << "Constraint values: " << desired_values << '\n';

    // Currently, all fixed vertex constraints normalized
#ifndef NDEBUG
    for (int i = 0; i < numconstraints; ++i)
        assert(approxEq(normal[i].norm(), 1.0, 1.0e-9));
#endif

    std::vector<VecXd> ntilde;
    for (int i = 0; i < numconstraints; ++i)
    {
        ntilde.push_back(VecXd(ndof));
        ntilde.back().setZero();
        int status = solver->solve(ntilde.back(), normal[i]);
        if (status < 0)
            DebugStream(g_log, "")
                    << "\033[31;1mWARNING IN IMPLICITEULER:\033[m Problem during linear solve detected. exertCompliantInelasticVertexFaceImpulse. Time: "
                    << m_t << '\n';
    }

    // Ensure the edge degrees of freedom experience no impulse
#ifndef NDEBUG
    for (ElasticRod::edge_iter eit = m_rods[rodidx]->edges_begin(); eit != m_rods[rodidx]->edges_end(); ++eit)
    {
        for (int i = 0; i < numconstraints; ++i)
            assert(ntilde[i](m_rods[rodidx]->edgeIdx(*eit)) == 0.0);
    }
#endif

    // Vectors restriced to just vertex DoFs. TODO: do all of this in place in the vectors I already have
    std::vector<VecXd> posnn;
    for (int i = 0; i < numconstraints; ++i)
        posnn.push_back(VecXd(nvdof));
#ifndef NDEBUG
    for (int i = 0; i < numconstraints; ++i)
        posnn[i].setConstant(std::numeric_limits<double>::signaling_NaN());
#endif

    std::vector<VecXd> posnntilde;
    for (int i = 0; i < numconstraints; ++i)
        posnntilde.push_back(VecXd(nvdof));
#ifndef NDEBUG
    for (int i = 0; i < numconstraints; ++i)
        posnntilde[i].setConstant(std::numeric_limits<double>::signaling_NaN());
#endif

    for (int i = 0; i < numconstraints; ++i)
    {
        int theidx = 0;
        for (ElasticRod::vertex_iter vit = m_rods[rodidx]->vertices_begin(); vit != m_rods[rodidx]->vertices_end(); ++vit)
        {
            int vert_dof = m_rods[rodidx]->vertIdx(*vit, 0);
            posnn[i].segment<3> (3 * theidx) = normal[i].segment<3> (vert_dof);
            ++theidx;
        }
    }

    for (int i = 0; i < numconstraints; ++i)
    {
        int theidx = 0;
        for (ElasticRod::vertex_iter vit = m_rods[rodidx]->vertices_begin(); vit != m_rods[rodidx]->vertices_end(); ++vit)
        {
            int vert_dof = m_rods[rodidx]->vertIdx(*vit, 0);
            posnntilde[i].segment<3> (3 * theidx) = ntilde[i].segment<3> (vert_dof);
            ++theidx;
        }
    }

    // Compute the Lagrange multipliers
    // TODO: This matrix is symmetric, exploit that fact to avoid computations
    MatXd lglhs(numconstraints, numconstraints);
#ifndef NDEBUG
    lglhs.setConstant(std::numeric_limits<double>::signaling_NaN());
#endif
    for (int i = 0; i < numconstraints; ++i)
        for (int j = 0; j < numconstraints; ++j)
            lglhs(i, j) = posnn[i].dot(posnntilde[j]);
    assert(approxSymmetric(lglhs, 1.0e-6));

    Eigen::VectorXd lgrhs(numconstraints);
#ifndef NDEBUG
    lgrhs.setConstant(std::numeric_limits<double>::signaling_NaN());
#endif
    lgrhs(0) = vfcol.computeRelativeVelocity(); //posnN0.dot(m_vnphalf.segment(m_base_dof_indices[rodidx],posnN0.size()));
    assert(lgrhs(0) < 0.0);
    for (int i = 1; i < numconstraints; ++i)
        lgrhs(i) = posnn[i].dot(m_vnphalf.segment(rodbase, nvdof));
    lgrhs *= -1.0;

    //TraceStream(g_log, "") << lgrhs << '\n';

    lgrhs += desired_values;

    //TraceStream(g_log, "") << lgrhs << '\n';
    //TraceStream(g_log, "") << desired_values << '\n';

    Eigen::VectorXd alpha(numconstraints); // = lglhs.inverse()*lgrhs;
    alpha = lglhs.lu().solve(lgrhs);

    assert(alpha(0) >= 0.0);

    double magrelvel = vfcol.computeRelativeVelocity();

    // TraceStream(g_log, "") << "BARodStepper::exertCompliantInelasticVertexFaceImpulse: pre-impulse view of collision: " << vfcol
    //         << '\n';
    //TraceStream(g_log, "") << "BARodStepper::exertCompliantInelasticVertexFaceImpulse: pre-impulse velocities: " << m_vnphalf.segment(rodbase, nvdof) << '\n';

    for (int i = 0; i < numconstraints; ++i)
        m_vnphalf.segment(rodbase, nvdof) += alpha(i) * posnntilde[i];

    // TraceStream(g_log, "") << "BARodStepper::exertCompliantInelasticVertexFaceImpulse: post-impulse view of collision: " << vfcol
    //         << '\n';
    //TraceStream(g_log, "") << "BARodStepper::exertCompliantInelasticVertexFaceImpulse: post-impulse velocities: " << m_vnphalf.segment(rodbase, nvdof) << '\n';

    // Ensure that the impulse eliminated the realtive velocity
    //#ifndef NDEBUG
    double postmagrelvel = vfcol.computeRelativeVelocity();
    // TraceStream(g_log, "") << "BARodStepper::exertCompliantInelasticVertexFaceImpulse: relative velocity pre-impulse = " << magrelvel
    //         << " post-impulse = " << postmagrelvel << '\n';
    // Ensure the inelastic impulse decreased the realtive velocity
    assert(fabs(postmagrelvel) <= fabs(magrelvel));
    // Should add some 'extra kick' to ensure collision gets killed, but for now just be content with small velocity
    assert(fabs(postmagrelvel) < 1.0e-9);
    //assert( postmagrelvel < 0.0 );
    //#endif

    applyInextensibilityVelocityFilter(rodidx);

    // Ensure the 'scripted vertices' achieved the desired velocity
#ifndef NDEBUG
    curdof = 1;
    for (size_t i = 0; i < scriptedverts.size(); ++i)
    {
        assert(scriptedverts[i] >= 0);
        assert((int) scriptedverts[i] < m_rods[rodidx]->nv());
        //TraceStream(g_log, "") << m_vnphalf(rodbase+3*scriptedverts[i]+0) << "   " << desired_values(curdof) << '\n';
        assert(approxEq(m_vnphalf(rodbase + 3 * scriptedverts[i] + 0), desired_values(curdof), 1.0e-6));
        ++curdof;
        //TraceStream(g_log, "") << m_vnphalf(rodbase+3*scriptedverts[i]+1) << "   " << desired_values(curdof) << '\n';
        assert(approxEq(m_vnphalf(rodbase + 3 * scriptedverts[i] + 1), desired_values(curdof), 1.0e-6));
        ++curdof;
        //TraceStream(g_log, "") << m_vnphalf(rodbase+3*scriptedverts[i]+2) << "   " << desired_values(curdof) << '\n';
        assert(approxEq(m_vnphalf(rodbase + 3 * scriptedverts[i] + 2), desired_values(curdof), 1.0e-6));
        ++curdof;
    }
#endif
}

void BARodStepper::exertCompliantInelasticEdgeEdgeImpulse(const EdgeEdgeCTCollision& eecol)
{
    //   TraceStream(g_log, "") << "Edge-edge compliant inelastic impulse" << '\n';
    //   TraceStream(g_log, "") << eecol << '\n';

    // Determine if either edge is totally fixed
    bool rod0fixed = YAEdge(eecol.e0_v0, eecol.e0_v1).IsFixed(m_geodata);
    bool rod1fixed = YAEdge(eecol.e1_v0, eecol.e1_v1).IsFixed(m_geodata);

    // Skip collisions in which both edges are fixed
    if (rod0fixed && rod1fixed)
        return;
    // Neither rod is totally fixed

    else if (!rod0fixed && !rod1fixed)
        exertCompliantInelasticEdgeEdgeImpulseBothFree(eecol);
    // One rod is totally fixed

    else
        exertCompliantInelasticEdgeEdgeImpulseOneFixed(eecol);
}

void BARodStepper::exertCompliantInelasticEdgeEdgeImpulseBothFree(const EdgeEdgeCTCollision& eecol)
{
    // Ensure both edges have two free vertices
    //assert( isEntireEdgeFree(eecol.e0_v0,eecol.e0_v1) );
    //assert( isEntireEdgeFree(eecol.e1_v0,eecol.e1_v1) );

    assert(!YAEdge(eecol.e0_v0, eecol.e0_v1).IsFixed(m_geodata) && !YAEdge(eecol.e1_v0, eecol.e1_v1).IsFixed(m_geodata));

    // TraceStream(g_log, "") << "BARodStepper::exertCompliantInelasticEdgeEdgeImpulseBothFree: (x[" << eecol.e0_v0 << "]="
    //         << m_geodata.GetPoint(eecol.e0_v0) << " - x[" << eecol.e0_v1 << "]=" << m_geodata.GetPoint(eecol.e0_v1) << ",   x["
    //         << eecol.e1_v0 << "]=" << m_geodata.GetPoint(eecol.e1_v0) << " - x[" << eecol.e1_v1 << "]=" << m_geodata.GetPoint(
    //         eecol.e1_v1) << ")" << " s=" << eecol.s << " t=" << eecol.t << '\n';

    // Determine which rod each edge belongs to
    int rod0 = getContainingRod(eecol.e0_v0);
    assert(rod0 == getContainingRod(eecol.e0_v1));
    int rod1 = getContainingRod(eecol.e1_v0);
    assert(rod1 == getContainingRod(eecol.e1_v1));

    // Don't do self-collisions, for now
    assert(rod0 != rod1);

    if (!m_steppers[rod0]->HasSolved() || !m_steppers[rod1]->HasSolved())
    { // This should never happen because rod-rod collisions preclude selective adaptivity.
        DebugStream(g_log, "") << "WARNING: attempt to do rod-rod collision on non-dependable rods";
        return;
    }

    TraceStream(g_log, "") << "CTcollision: rod " << rod0 << " vs. rod " << rod1 << '\n';

    // Compute the relative velocity, which must be negative
    //    Vec3d relvel = m_geodata.computeRelativeVelocity(eecol.e0_v0, eecol.e0_v1, eecol.e1_v0, eecol.e1_v1, eecol.s, eecol.t);
    //    double magrelvel = relvel.dot(eecol.n);
    //   double magrelvel = eecol.computeRelativeVelocity(m_geodata);
    assert(eecol.computeRelativeVelocity() < 0.0);

    // Determine the total number of degrees of freedom in each rod
    int rod0ndof = m_rods[rod0]->ndof();
    int rod1ndof = m_rods[rod1]->ndof();

    // Determine the number of vertices in each rod
    int rod0nv = m_rods[rod0]->nv();
    int rod1nv = m_rods[rod1]->nv();

    // Determine the number of vertices in each rod
    int rod0nvdof = 3 * rod0nv;
    int rod1nvdof = 3 * rod1nv;

    //TraceStream(g_log, "") << "Ndof 0: " << rod0ndof << '\n';
    //TraceStream(g_log, "") << "Ndof 1: " << rod1ndof << '\n';
    //TraceStream(g_log, "") << "NV 0: " << rod0nv << '\n';
    //TraceStream(g_log, "") << "NV 1: " << rod1nv << '\n';
    //TraceStream(g_log, "") << "NVdof 0: " << rod0nvdof << '\n';
    //TraceStream(g_log, "") << "NVdof 1: " << rod1nvdof << '\n';

    // Determine where in the 'global' vertex dof pool each rod begins
    int rod0base = m_base_dof_indices[rod0];
    assert(rod0base % 3 == 0);
    int rod1base = m_base_dof_indices[rod1];
    assert(rod1base % 3 == 0);

    //TraceStream(g_log, "") << "rod0base: " << rod0base << '\n';
    //TraceStream(g_log, "") << "rod1base: " << rod1base << '\n';

    // Vertex number in 'position N'
    int i0 = eecol.e0_v0 - rod0base / 3;
    assert(i0 >= 0);
    assert(i0 < m_rods[rod0]->nv());
    int i1 = eecol.e0_v1 - rod0base / 3;
    assert(i1 >= 0);
    assert(i1 < m_rods[rod0]->nv());
    int j0 = eecol.e1_v0 - rod1base / 3;
    assert(j0 >= 0);
    assert(j0 < m_rods[rod1]->nv());
    int j1 = eecol.e1_v1 - rod1base / 3;
    assert(j1 >= 0);
    assert(j1 < m_rods[rod1]->nv());

    // Base indices within each rod's storage
    int base_i0 = m_rods[rod0]->vertIdx(i0, 0);
    assert(base_i0 >= 0);
    assert(base_i0 < rod0ndof);
    int base_i1 = m_rods[rod0]->vertIdx(i1, 0);
    assert(base_i1 >= 0);
    assert(base_i1 < rod0ndof);
    int base_j0 = m_rods[rod1]->vertIdx(j0, 0);
    assert(base_j0 >= 0);
    assert(base_j0 < rod1ndof);
    int base_j1 = m_rods[rod1]->vertIdx(j1, 0);
    assert(base_j1 >= 0);
    assert(base_j1 < rod1ndof);

    // Constraint normals for rod 0
    std::vector<VecXd> n0;
    // Collision constraint
    n0.push_back(VecXd(rod0ndof));
    n0.back().setZero();
    n0.back().segment<3> (base_i0) = -(1.0 - eecol.s) * eecol.GetNormal();
    n0.back().segment<3> (base_i1) = -eecol.s * eecol.GetNormal();
    // 'Scripted vertex' constraints
    const std::vector<int>& scriptedverts0 = m_rods[rod0]->getBoundaryCondition()->scriptedVertices();
    for (size_t i = 0; i < scriptedverts0.size(); ++i)
    {
        for (int j = 0; j < 3; ++j)
        {
            n0.push_back(VecXd(rod0ndof));
            n0.back().setZero();
            // Zero the velocity of this particular dof of this particular vertex
            n0.back()(m_rods[rod0]->vertIdx(scriptedverts0[i], j)) = 1.0;
        }
    }
    // Determine the desired value of each constraint
    int nc0 = (int) (n0.size());
    VecXd cval0(nc0);
#ifndef NDEBUG
    cval0.setConstant(std::numeric_limits<double>::signaling_NaN());
#endif
    // First constraint, the collision constraint, is 0
    cval0(0) = 0.0;
    // The 'scripted vertex' constraints are the respective components of the velocity
    int curdof = 1;
    for (size_t i = 0; i < scriptedverts0.size(); ++i)
    {
        assert(scriptedverts0[i] >= 0);
        assert((int) scriptedverts0[i] < m_rods[rod0]->nv());
        cval0(curdof++) = m_vnphalf(rod0base + 3 * scriptedverts0[i] + 0);
        cval0(curdof++) = m_vnphalf(rod0base + 3 * scriptedverts0[i] + 1);
        cval0(curdof++) = m_vnphalf(rod0base + 3 * scriptedverts0[i] + 2);
    }
    assert(curdof == nc0);
    assert(cval0.size() == nc0);
    // Ensure collision constraint adds up to actual normal
#ifndef NDEBUG
    Vec3d testn0 = n0[0].segment<3> (base_i0) + n0[0].segment<3> (base_i1);
    Vec3d actln0 = -eecol.GetNormal();
    assert(approxEq(testn0, actln0, 1.0e-6));
#endif
    // Currently, all fixed vertex constraints normalized
#ifndef NDEBUG
    for (int i = 1; i < nc0; ++i)
        assert(approxEq(n0[i].norm(), 1.0, 1.0e-9));
#endif
    // Ensure the edge degrees of freedom experience no impulse
#ifndef NDEBUG
    for (int i = 0; i < nc0; ++i)
        for (ElasticRod::edge_iter eit = m_rods[rod0]->edges_begin(); eit != m_rods[rod0]->edges_end(); ++eit)
            assert(n0[i](m_rods[rod0]->edgeIdx(*eit)) == 0.0);
#endif

    //for( int i = 0; i < nc0; ++i ) TraceStream(g_log, "") << "n0:    " << n0[i] << '\n';
    //for( int i = 0; i < nc0; ++i ) TraceStream(g_log, "") << "cval0: " << cval0[i] << '\n';

    // Constraint normals for rod 1
    std::vector<VecXd> n1;
    // Collision constraint
    n1.push_back(VecXd(rod1ndof));
    n1.back().setZero();
    n1.back().segment<3> (base_j0) = (1.0 - eecol.t) * eecol.GetNormal();
    n1.back().segment<3> (base_j1) = eecol.t * eecol.GetNormal();
    // 'Scripted vertex' constraints
    const std::vector<int>& scriptedverts1 = m_rods[rod1]->getBoundaryCondition()->scriptedVertices();
    for (size_t i = 0; i < scriptedverts1.size(); ++i)
    {
        for (int j = 0; j < 3; ++j)
        {
            n1.push_back(VecXd(rod1ndof));
            n1.back().setZero();
            // Zero the velocity of this particular dof of this particular vertex
            n1.back()(m_rods[rod1]->vertIdx(scriptedverts1[i], j)) = 1.0;
        }
    }
    // Determine the desired value of each constraint
    int nc1 = (int) (n1.size());
    VecXd cval1(nc1);
#ifndef NDEBUG
    cval1.setConstant(std::numeric_limits<double>::signaling_NaN());
#endif
    // First constraint, the collision constraint, is 0
    cval1(0) = 0.0;
    // The 'scripted vertex' constraints are the respective components of the velocity
    curdof = 1;
    for (size_t i = 0; i < scriptedverts1.size(); ++i)
    {
        assert(scriptedverts1[i] >= 0);
        assert((int) scriptedverts1[i] < m_rods[rod1]->nv());
        cval1(curdof++) = m_vnphalf(rod1base + 3 * scriptedverts1[i] + 0);
        cval1(curdof++) = m_vnphalf(rod1base + 3 * scriptedverts1[i] + 1);
        cval1(curdof++) = m_vnphalf(rod1base + 3 * scriptedverts1[i] + 2);
    }
    assert(curdof == nc1);
    assert(cval1.size() == nc1);
    // Ensure collision constraint adds up to actual normal
#ifndef NDEBUG
    Vec3d testn1 = n1[0].segment<3> (base_j0) + n1[0].segment<3> (base_j1);
    Vec3d actln1 = eecol.GetNormal();
    assert(approxEq(testn1, actln1, 1.0e-6));
#endif
    // Currently, all fixed vertex constraints normalized
#ifndef NDEBUG
    for (int i = 1; i < nc1; ++i)
        assert(approxEq(n1[i].norm(), 1.0, 1.0e-9));
#endif
    // Ensure the edge degrees of freedom experience no impulse
#ifndef NDEBUG
    for (int i = 0; i < nc1; ++i)
        for (ElasticRod::edge_iter eit = m_rods[rod1]->edges_begin(); eit != m_rods[rod1]->edges_end(); ++eit)
            assert(n1[i](m_rods[rod1]->edgeIdx(*eit)) == 0.0);
#endif

    //for( int i = 0; i < nc1; ++i ) TraceStream(g_log, "") << "n1:    " << n1[i] << '\n';
    //for( int i = 0; i < nc1; ++i ) TraceStream(g_log, "") << "cval1: " << cval1[i] << '\n';

    // Compute ntilde for each rod 0 constraint
    LinearSystemSolver* lss = m_solver_collection.getLinearSystemSolver(rod0ndof);
    MatrixBase* lhs0 = lss->m_lhs;
    assert(lhs0 != NULL);
    assert(lhs0->rows() == rod0ndof);
    assert(lhs0->cols() == rod0ndof);
    LinearSolverBase* solver0 = lss->m_solver;
    assert(solver0 != NULL);
    computeCompliantLHS(lhs0, rod0);
    std::vector<VecXd> ntilde0;
    for (int i = 0; i < nc0; ++i)
    {
        ntilde0.push_back(VecXd(rod0ndof));
        ntilde0.back().setZero();
        int status = solver0->solve(ntilde0.back(), n0[i]);
        if (status < 0)
            DebugStream(g_log, "") << "\033[31;1mWARNING IN IMPLICITEULER:\033[m Problem during linear solve detected. "
                    << '\n';
    }
    // Ensure the edge degrees of freedom experience no impulse
#ifndef NDEBUG
    for (int i = 0; i < nc0; ++i)
        for (ElasticRod::edge_iter eit = m_rods[rod0]->edges_begin(); eit != m_rods[rod0]->edges_end(); ++eit)
            assert(ntilde0[i](m_rods[rod0]->edgeIdx(*eit)) == 0.0);
#endif

    //for( int i = 0; i < nc0; ++i ) TraceStream(g_log, "") << "ntilde0: " << ntilde0[i] << '\n';

    // Compute ntilde for each rod 1 constraint
    lss = m_solver_collection.getLinearSystemSolver(rod1ndof);
    MatrixBase* lhs1 = lss->m_lhs;
    assert(lhs1->rows() == rod1ndof);
    assert(lhs1->cols() == rod1ndof);
    assert(lhs1 != NULL);
    LinearSolverBase* solver1 = lss->m_solver;
    assert(solver1 != NULL);
    computeCompliantLHS(lhs1, rod1);
    std::vector<VecXd> ntilde1;
    for (int i = 0; i < nc1; ++i)
    {
        ntilde1.push_back(VecXd(rod1ndof));
        ntilde1.back().setZero();
        int status = solver1->solve(ntilde1.back(), n1[i]);
        if (status < 0)
            DebugStream(g_log, "") << "\033[31;1mWARNING IN IMPLICITEULER:\033[m Problem during linear solve detected. "
                    << '\n';
    }
    // Ensure the edge degrees of freedom experience no impulse
#ifndef NDEBUG
    for (int i = 0; i < nc1; ++i)
        for (ElasticRod::edge_iter eit = m_rods[rod1]->edges_begin(); eit != m_rods[rod1]->edges_end(); ++eit)
            assert(ntilde1[i](m_rods[rod1]->edgeIdx(*eit)) == 0.0);
#endif

    //for( int i = 0; i < nc1; ++i ) TraceStream(g_log, "") << "ntilde1: " << ntilde1[i] << '\n';

    // Vectors restriced to just vertex DoFs for rod 0. TODO: do all of this in place in the vectors I already have
    std::vector<VecXd> posnn0;
    for (int i = 0; i < nc0; ++i)
        posnn0.push_back(VecXd(rod0nvdof));
#ifndef NDEBUG
    for (int i = 0; i < nc0; ++i)
        posnn0[i].setConstant(std::numeric_limits<double>::signaling_NaN());
#endif

    std::vector<VecXd> posnntilde0;
    for (int i = 0; i < nc0; ++i)
        posnntilde0.push_back(VecXd(rod0nvdof));
#ifndef NDEBUG
    for (int i = 0; i < nc0; ++i)
        posnntilde0[i].setConstant(std::numeric_limits<double>::signaling_NaN());
#endif

    for (int i = 0; i < nc0; ++i)
    {
        int theidx = 0;
        for (ElasticRod::vertex_iter vit = m_rods[rod0]->vertices_begin(); vit != m_rods[rod0]->vertices_end(); ++vit)
        {
            int vert_dof = m_rods[rod0]->vertIdx(*vit, 0);
            posnn0[i].segment<3> (3 * theidx) = n0[i].segment<3> (vert_dof);
            ++theidx;
        }
    }

    //for( int i = 0; i < nc0; ++i ) TraceStream(g_log, "") << "posnn0: " << posnn0[i] << '\n';

    for (int i = 0; i < nc0; ++i)
    {
        int theidx = 0;
        for (ElasticRod::vertex_iter vit = m_rods[rod0]->vertices_begin(); vit != m_rods[rod0]->vertices_end(); ++vit)
        {
            int vert_dof = m_rods[rod0]->vertIdx(*vit, 0);
            posnntilde0[i].segment<3> (3 * theidx) = ntilde0[i].segment<3> (vert_dof);
            ++theidx;
        }
    }

    //for( int i = 0; i < nc0; ++i ) TraceStream(g_log, "") << "posnntilde0: " << posnntilde0[i] << '\n';


    // Vectors restriced to just vertex DoFs for rod 1. TODO: do all of this in place in the vectors I already have
    std::vector<VecXd> posnn1;
    for (int i = 0; i < nc1; ++i)
        posnn1.push_back(VecXd(rod1nvdof));
#ifndef NDEBUG
    for (int i = 0; i < nc1; ++i)
        posnn1[i].setConstant(std::numeric_limits<double>::signaling_NaN());
#endif

    std::vector<VecXd> posnntilde1;
    for (int i = 0; i < nc1; ++i)
        posnntilde1.push_back(VecXd(rod1nvdof));
#ifndef NDEBUG
    for (int i = 0; i < nc1; ++i)
        posnntilde1[i].setConstant(std::numeric_limits<double>::signaling_NaN());
#endif

    for (int i = 0; i < nc1; ++i)
    {
        int theidx = 0;
        for (ElasticRod::vertex_iter vit = m_rods[rod1]->vertices_begin(); vit != m_rods[rod1]->vertices_end(); ++vit)
        {
            int vert_dof = m_rods[rod1]->vertIdx(*vit, 0);
            posnn1[i].segment<3> (3 * theidx) = n1[i].segment<3> (vert_dof);
            ++theidx;
        }
    }

    //for( int i = 0; i < nc1; ++i ) TraceStream(g_log, "") << "posnn1: " << posnn1[i] << '\n';

    for (int i = 0; i < nc1; ++i)
    {
        int theidx = 0;
        for (ElasticRod::vertex_iter vit = m_rods[rod1]->vertices_begin(); vit != m_rods[rod1]->vertices_end(); ++vit)
        {
            int vert_dof = m_rods[rod1]->vertIdx(*vit, 0);
            posnntilde1[i].segment<3> (3 * theidx) = ntilde1[i].segment<3> (vert_dof);
            ++theidx;
        }
    }

    //for( int i = 0; i < nc1; ++i ) TraceStream(g_log, "") << "posnntilde1: " << posnntilde1[i] << '\n';


    // Compute the Lagrange multipliers for all constraints
    // TODO: This matrix is symmetric, exploit that fact to avoid computations
    int numalpha = nc0 + nc1 - 1;
    MatXd lglhs(numalpha, numalpha);
    //#ifndef NDEBUG
    //  lglhs.setConstant(std::numeric_limits<double>::signaling_NaN());
    //#endif
    lglhs.setZero();

    //TraceStream(g_log, "") << "System size: " << numalpha << '\n';

    // TODO: Make sure indexing is correct here
    // Entry (0,0)
    lglhs(0, 0) = posnn0[0].dot(posnntilde0[0]) + posnn1[0].dot(posnntilde1[0]);
    // Entires (0,1::numalpha)
    int curel = 1;
    for (int i = 1; i < nc0; ++i)
        lglhs(0, curel++) = posnn0[0].dot(posnntilde0[i]);
    // Entries (0,numalpha::numalpha+numbeta)
    for (int i = 1; i < nc1; ++i)
        lglhs(0, curel++) = posnn1[0].dot(posnntilde1[i]);

    // Entires (1::numalpha,0)
    curel = 1;
    for (int i = 1; i < nc0; ++i)
        lglhs(curel++, 0) = posnntilde0[0].dot(posnn0[i]);
    // Entries (numalpha::numalpha+numbeta,0)
    for (int i = 1; i < nc1; ++i)
        lglhs(curel++, 0) = posnntilde1[0].dot(posnn1[i]);

    // Entries (1::numalpha,1::numalpha)
    for (int i = 1; i < nc0; ++i)
        for (int j = 1; j < nc0; ++j)
            lglhs(i, j) = posnn0[i].dot(posnntilde0[j]);
    // Entries (numalpha::numalpha+numbeta,numalpha::numalpha+numbeta)
    for (int i = 1; i < nc1; ++i)
        for (int j = 1; j < nc1; ++j)
            lglhs(nc0 + i - 1, nc0 + j - 1) = posnn1[i].dot(posnntilde1[j]);

    assert(approxSymmetric(lglhs, 1.0e-6));

    Eigen::VectorXd lgrhs(numalpha);
#ifndef NDEBUG
    lgrhs.setConstant(std::numeric_limits<double>::signaling_NaN());
#endif
    // Entry 0
    assert(posnn0[0].size() == m_vnphalf.segment(rod0base, rod0nvdof).size());
    assert(posnn1[0].size() == m_vnphalf.segment(rod1base, rod1nvdof).size());
    lgrhs(0) = posnn0[0].dot(m_vnphalf.segment(rod0base, rod0nvdof)) + posnn1[0].dot(m_vnphalf.segment(rod1base, rod1nvdof));
    assert(approxEq(lgrhs(0), eecol.computeRelativeVelocity(), 1.0e-6));
    assert(lgrhs(0) < 0.0);

    // Entries 1...numalpha
    for (int i = 1; i < nc0; ++i)
    {
        assert(posnn0[i].size() == m_vnphalf.segment(rod0base, rod0nvdof).size());
        lgrhs(i) = posnn0[i].dot(m_vnphalf.segment(rod0base, rod0nvdof));
    }
    // Entries numalpha...numalpha+numbeta
    for (int i = 1; i < nc1; ++i)
    {
        assert(posnn1[i].size() == m_vnphalf.segment(rod1base, rod1nvdof).size());
        lgrhs(nc0 + i - 1) = posnn1[i].dot(m_vnphalf.segment(rod1base, rod1nvdof));
    }
    //TraceStream(g_log, "") << lgrhs << '\n';
    lgrhs *= -1.0;

    // For now, model an inelastic collision
    //lglhs0(0) += 0.0;
    // Entries 1...numalpha
    lgrhs.segment(1, nc0 - 1) += cval0.segment(1, nc0 - 1);
    // Entries numalpha...numalpha+numbeta
    lgrhs.segment(nc0, nc1 - 1) += cval1.segment(1, nc1 - 1);
    //TraceStream(g_log, "") << lgrhs << '\n';

    Eigen::VectorXd alpha(numalpha);
    assert(lglhs.rows() == lglhs.cols());
    assert(lglhs.rows() == lgrhs.size());
    assert(lglhs.rows() == alpha.size());
    alpha = lglhs.lu().solve(lgrhs);

    assert(alpha(0) >= 0.0);

    if (alpha(0) < 0.0)
        DebugStream(g_log, "") << "WARNING NEGATIVE LAGRANGE MULTIPLIER alpha: " << alpha << '\n';

    // BEGIN TEMP
    //if( rod0 == 23 || rod1 == 23 ) DebugStream(g_log, "") << "Collision lagrange multiplier: " << alpha(0) << '\n';
    // END TEMP

    // Add the impulse corresponding to the collision to rod 0
    assert(m_vnphalf.segment(rod0base, rod0nvdof).size() == posnntilde0[0].size());
    m_vnphalf.segment(rod0base, rod0nvdof) += alpha(0) * posnntilde0[0];
    assert(m_vnphalf.segment(rod1base, rod1nvdof).size() == posnntilde1[0].size());
    m_vnphalf.segment(rod1base, rod1nvdof) += alpha(0) * posnntilde1[0];
    // Add the rod 0 specific impulses to rod 0
    for (int i = 1; i < nc0; ++i)
    {
        assert(m_vnphalf.segment(rod0base, rod0nvdof).size() == posnntilde0[i].size());
        m_vnphalf.segment(rod0base, rod0nvdof) += alpha(i) * posnntilde0[i];
    }
    // Add the rod 1 specific impulses to rod 1
    for (int i = 1; i < nc1; ++i)
    {
        assert(m_vnphalf.segment(rod1base, rod1nvdof).size() == posnntilde1[i].size());
        m_vnphalf.segment(rod1base, rod1nvdof) += alpha(nc0 + i - 1) * posnntilde1[i];
    }

    // Ensure that the impulse eliminated the realtive velocity
#ifndef NDEBUG
    //Vec3d postrelvel = computeRelativeVelocity( m_vnphalf, eecol.e0_v0, eecol.e0_v1, eecol.e1_v0, eecol.e1_v1, eecol.s, eecol.t );
    //double postmagrelvel = postrelvel.dot(eecol.n);
    // Ensure the inelastic impulse decreased the realtive velocity
    //assert( fabs(postmagrelvel) <= fabs(magrelvel) );
    // Should add some 'extra kick' to ensure collision gets killed, but for now just be content with small velocity
    //assert( fabs(postmagrelvel) < 1.0e-9 );
    //TraceStream(g_log, "") << magrelvel << "   " << postmagrelvel << '\n';
    //assert( postmagrelvel < 0.0 );
#endif

    // Ensure the 'scripted vertices' achieved the desired velocity for rod 0
#ifndef NDEBUG
    curdof = 1;
    for (size_t i = 0; i < scriptedverts0.size(); ++i)
    {
        assert(scriptedverts0[i] >= 0);
        assert((int) scriptedverts0[i] < m_rods[rod0]->nv());
        //TraceStream(g_log, "") << m_vnphalf(rodbase+3*scriptedverts[i]+0) << "   " << desired_values(curdof) << '\n';
        assert(approxEq(m_vnphalf(rod0base + 3 * scriptedverts0[i] + 0), cval0(curdof), 1.0e-6));
        ++curdof;
        //TraceStream(g_log, "") << m_vnphalf(rodbase+3*scriptedverts[i]+1) << "   " << desired_values(curdof) << '\n';
        assert(approxEq(m_vnphalf(rod0base + 3 * scriptedverts0[i] + 1), cval0(curdof), 1.0e-6));
        ++curdof;
        //TraceStream(g_log, "") << m_vnphalf(rodbase+3*scriptedverts[i]+2) << "   " << desired_values(curdof) << '\n';
        assert(approxEq(m_vnphalf(rod0base + 3 * scriptedverts0[i] + 2), cval0(curdof), 1.0e-6));
        ++curdof;
    }
#endif

    // Ensure the 'scripted vertices' achieved the desired velocity for rod 1
#ifndef NDEBUG
    curdof = 1;
    for (size_t i = 0; i < scriptedverts1.size(); ++i)
    {
        assert(scriptedverts1[i] >= 0);
        assert((int) scriptedverts1[i] < m_rods[rod1]->nv());
        //TraceStream(g_log, "") << m_vnphalf(rodbase+3*scriptedverts[i]+0) << "   " << desired_values(curdof) << '\n';
        assert(approxEq(m_vnphalf(rod1base + 3 * scriptedverts1[i] + 0), cval1(curdof), 1.0e-6));
        ++curdof;
        //TraceStream(g_log, "") << m_vnphalf(rodbase+3*scriptedverts[i]+1) << "   " << desired_values(curdof) << '\n';
        assert(approxEq(m_vnphalf(rod1base + 3 * scriptedverts1[i] + 1), cval1(curdof), 1.0e-6));
        ++curdof;
        //TraceStream(g_log, "") << m_vnphalf(rodbase+3*scriptedverts[i]+2) << "   " << desired_values(curdof) << '\n';
        assert(approxEq(m_vnphalf(rod1base + 3 * scriptedverts1[i] + 2), cval1(curdof), 1.0e-6));
        ++curdof;
    }
#endif
}

void BARodStepper::exertCompliantInelasticEdgeEdgeImpulseOneFixed(const EdgeEdgeCTCollision& eecol)
{
    // Must have one totally fixed and one totally free edge
    // assert(
    //         (YAEdge(eecol.e0_v0, eecol.e0_v1).IsFree(m_geodata) && YAEdge(eecol.e1_v0, eecol.e1_v1).IsFixed(m_geodata))
    //                 || (YAEdge(eecol.e1_v0, eecol.e1_v1).IsFree(m_geodata) && YAEdge(eecol.e0_v0, eecol.e0_v1).IsFixed(
    //                         m_geodata)));

    // TraceStream(g_log, "") << "BARodStepper::exertCompliantInelasticEdgeEdgeImpulseOneFixed: (x[" << eecol.e0_v0 << "]="
    //        << m_geodata.GetPoint(eecol.e0_v0) << " - x[" << eecol.e0_v1 << "]=" << m_geodata.GetPoint(eecol.e0_v1) << ",   x["
    //        << eecol.e1_v0 << "]=" << m_geodata.GetPoint(eecol.e1_v0) << " - x[" << eecol.e1_v1 << "]=" << m_geodata.GetPoint(
    //       eecol.e1_v1) << ")" << " s=" << eecol.s << " t=" << eecol.t << '\n';

    // Determine if either edge is fixed
    bool rod0fixed = YAEdge(eecol.e0_v0, eecol.e0_v1).IsFixed(m_geodata);
    bool rod1fixed = YAEdge(eecol.e1_v0, eecol.e1_v1).IsFixed(m_geodata);
    assert(rod0fixed != rod1fixed);

    double preRelativeVelocity = eecol.computeRelativeVelocity();

    // Compute the relative velocity, which must be negative for a collision to have happened
    //    Vec3d relvel = m_geodata.computeRelativeVelocity(eecol.e0_v0, eecol.e0_v1, eecol.e1_v0, eecol.e1_v1, eecol.s, eecol.t);
    //    double magrelvel = relvel.dot(eecol.n);
    //   double magrelvel = eecol.computeRelativeVelocity(m_geodata);
    assert(eecol.computeRelativeVelocity() < 0.0);

    // Extract the rod index, barycentric coordinate of the collision, and indices of the edge involved in the collision
    int rodidx = -1;
    double u = -1.0;
    int v0 = -1;
    int v1 = -1;
    if (!rod0fixed)
    {
        rodidx = getContainingRod(eecol.e0_v0);
        u = eecol.s;
        v0 = eecol.e0_v0;
        v1 = eecol.e0_v1;
    }
    if (!rod1fixed)
    {
        rodidx = getContainingRod(eecol.e1_v0);
        u = eecol.t;
        v0 = eecol.e1_v0;
        v1 = eecol.e1_v1;
    }
    assert(u >= 0.0);
    assert(u <= 1.0);
    assert(rodidx >= 0);
    assert(rodidx < (int) m_number_of_rods);

    // If the rod has not solved properly, no need to compute its collision response
    if (!m_steppers[rodidx]->HasSolved())
    {
        DebugStream(g_log, "") << "WARNING: attempt to do edge-edge collision with non-dependable rod\n";
        return;
    }

    TraceStream(g_log, "") << "CTcollision: rod " << rodidx << " vs. mesh\n";

    // Convert the vertices' global indices to rod indices
    assert(m_base_dof_indices[rodidx] % 3 == 0);
    v0 = v0 - m_base_dof_indices[rodidx] / 3;
    v1 = v1 - m_base_dof_indices[rodidx] / 3;
    assert(v0 >= 0);
    assert(v0 < m_rods[rodidx]->nv());
    assert(v1 >= 0);
    assert(v1 < m_rods[rodidx]->nv());

    // Get storage for lhs of linear system, get a solver
    LinearSystemSolver* lss = m_solver_collection.getLinearSystemSolver(m_rods[rodidx]->ndof());
    MatrixBase* lhs = lss->m_lhs;
    assert(lhs != NULL);
    LinearSolverBase* solver = lss->m_solver;
    assert(solver != NULL);

    // Compute M - h^2*dF/dx
    computeCompliantLHS(lhs, rodidx);

    int rodbase = m_base_dof_indices[rodidx];

    // Compute the 'base' index of each vertex
    int base0 = m_rods[rodidx]->vertIdx(v0, 0);
    assert(base0 % 4 == 0);
    int base1 = m_rods[rodidx]->vertIdx(v1, 0);
    assert(base1 % 4 == 0);

    int ndof = m_rods[rodidx]->ndof();
    int nvdof = 3 * m_rods[rodidx]->nv();

    // Compute the 'global' normal
    std::vector<VecXd> n;

    n.push_back(VecXd(ndof));
    n.back().setZero();
    n.back().segment<3> (base0) = (1 - u) * eecol.GetNormal();
    n.back().segment<3> (base1) = u * eecol.GetNormal();

    // If rod 0 has the free edge, need to invert sign of normal
    if (!rod0fixed)
    {
        n.back().segment<3> (base0) *= -1.0;
        n.back().segment<3> (base1) *= -1.0;
    }

    // Compute three normals per scripted vertex
    const std::vector<int>& scriptedverts = m_rods[rodidx]->getBoundaryCondition()->scriptedVertices();
    for (size_t i = 0; i < scriptedverts.size(); ++i)
    {
        for (int j = 0; j < 3; ++j)
        {
            //TraceStream(g_log, "") << m_rods[rodidx]->vertIdx(scriptedverts[i],j) << " ";
            n.push_back(VecXd(ndof));
            n.back().setZero();
            // Zero the velocity of this particular dof of this particular vertex
            n.back()(m_rods[rodidx]->vertIdx(scriptedverts[i], j)) = 1.0;
        }
    }
    //TraceStream(g_log, "") << '\n';

    int numconstraints = (int) (n.size());

    // Determine the desired values for each constraint
    VecXd desired_values(numconstraints);
#ifndef NDEBUG
    desired_values.setConstant(std::numeric_limits<double>::signaling_NaN());
#endif
    // First constraint, the collision constraint, is 0
    desired_values(0) = 0.0;
    // The 'scripted vertex' constraints are the respective components of the velocity
    int curdof = 1;
    for (size_t i = 0; i < scriptedverts.size(); ++i)
    {
        assert(scriptedverts[i] >= 0);
        assert((int) scriptedverts[i] < m_rods[rodidx]->nv());
        desired_values(curdof++) = m_vnphalf(rodbase + 3 * scriptedverts[i] + 0);
        desired_values(curdof++) = m_vnphalf(rodbase + 3 * scriptedverts[i] + 1);
        desired_values(curdof++) = m_vnphalf(rodbase + 3 * scriptedverts[i] + 2);
    }
    assert(curdof == numconstraints);
    assert(desired_values.size() == numconstraints);

    //TraceStream(g_log, "") << desired_values << '\n';

    //for( size_t i = 0; i < vertexConstraints.size(); ++i ) TraceStream(g_log, "") << vertexConstraints[i] << '\n';

    // Ensure collision constraint adds up to actual normal
#ifndef NDEBUG
    Vec3d testn = n[0].segment<3> (base0) + n[0].segment<3> (base1);
    if (!rod0fixed)
        testn *= -1.0;
    assert(approxEq(testn, eecol.GetNormal(), 1.0e-6));
#endif

    // Currently, all fixed vertex constraints normalized
#ifndef NDEBUG
    for (int i = 1; i < numconstraints; ++i)
        assert(approxEq(n[i].norm(), 1.0, 1.0e-9));
#endif

    std::vector<VecXd> ntilde;
    for (int i = 0; i < numconstraints; ++i)
    {
        ntilde.push_back(VecXd(ndof));
        ntilde.back().setZero();
        int status = solver->solve(ntilde.back(), n[i]);
        if (status < 0)
            DebugStream(g_log, "")
                    << "\033[31;1mWARNING IN IMPLICITEULER:\033[m Problem during linear solve detected. exertCompliantInelasticEdgeEdgeImpulseOneFixedThree. Time: "
                    << m_t << '\n';
    }

    //for( size_t i = 0; i < n.size(); ++i ) TraceStream(g_log, "") << ntilde[i].transpose() << '\n';


    // Ensure the edge degrees of freedom experience no impulse
#ifndef NDEBUG
    for (ElasticRod::edge_iter eit = m_rods[rodidx]->edges_begin(); eit != m_rods[rodidx]->edges_end(); ++eit)
    {
        for (int i = 0; i < numconstraints; ++i)
            assert(ntilde[i](m_rods[rodidx]->edgeIdx(*eit)) == 0.0);
    }
#endif

    // Vectors restriced to just verted DoFs. TODO: do all of this in place in the vectors I already have
    std::vector<VecXd> posnn;
    for (int i = 0; i < numconstraints; ++i)
        posnn.push_back(VecXd(nvdof));
#ifndef NDEBUG
    for (int i = 0; i < numconstraints; ++i)
        posnn[i].setConstant(std::numeric_limits<double>::signaling_NaN());
#endif

    std::vector<VecXd> posnntilde;
    for (int i = 0; i < numconstraints; ++i)
        posnntilde.push_back(VecXd(nvdof));
#ifndef NDEBUG
    for (int i = 0; i < numconstraints; ++i)
        posnntilde[i].setConstant(std::numeric_limits<double>::signaling_NaN());
#endif

    for (int i = 0; i < numconstraints; ++i)
    {
        int theidx = 0;
        for (ElasticRod::vertex_iter vit = m_rods[rodidx]->vertices_begin(); vit != m_rods[rodidx]->vertices_end(); ++vit)
        {
            int vert_dof = m_rods[rodidx]->vertIdx(*vit, 0);
            posnn[i].segment<3> (3 * theidx) = n[i].segment<3> (vert_dof);
            ++theidx;
        }
    }

    for (int i = 0; i < numconstraints; ++i)
    {
        int theidx = 0;
        for (ElasticRod::vertex_iter vit = m_rods[rodidx]->vertices_begin(); vit != m_rods[rodidx]->vertices_end(); ++vit)
        {
            int vert_dof = m_rods[rodidx]->vertIdx(*vit, 0);
            posnntilde[i].segment<3> (3 * theidx) = ntilde[i].segment<3> (vert_dof);
            ++theidx;
        }
    }

    // Compute the Lagrange multipliers
    // TODO: This matrix is symmetric, exploit that fact to avoid computations
    MatXd lglhs(numconstraints, numconstraints);
#ifndef NDEBUG
    lglhs.setConstant(std::numeric_limits<double>::signaling_NaN());
#endif
    for (int i = 0; i < numconstraints; ++i)
        for (int j = 0; j < numconstraints; ++j)
            lglhs(i, j) = posnn[i].dot(posnntilde[j]);
    assert(approxSymmetric(lglhs, 1.0e-6));

    Eigen::VectorXd lgrhs(numconstraints);
#ifndef NDEBUG
    lgrhs.setConstant(std::numeric_limits<double>::signaling_NaN());
#endif
    lgrhs(0) = eecol.computeRelativeVelocity(); // CRAZY SIGN PROBLEM HERE?
    assert(lgrhs(0) < 0.0);
    for (int i = 1; i < numconstraints; ++i)
        lgrhs(i) = posnn[i].dot(m_vnphalf.segment(rodbase, nvdof));
    lgrhs *= -1.0;

    //TraceStream(g_log, "") << lgrhs << '\n';

    lgrhs += desired_values;

    //TraceStream(g_log, "") << lgrhs << '\n';
    //TraceStream(g_log, "") << desired_values << '\n';

    Eigen::VectorXd alpha(numconstraints); // = lglhs.inverse()*lgrhs;
    alpha = lglhs.lu().solve(lgrhs);

    // Contact constraint should 'push not pull'
    assert(alpha(0) >= 0.0);

    // TEMP DEBUG STUFF
    //  if( rodidx == 20 && m_t >= 0.85 )
    //  {
    //    TraceStream(g_log, "") << "---------- ROD 20 ROD 20 ROD 20 ROD 20 ROD 20 ----------" << '\n';
    //    TraceStream(g_log, "") << "alpha: " << alpha << '\n';
    //    TraceStream(g_log, "") << "time: " << m_t << '\n';
    //
    //    if( !rod0fixed ) TraceStream(g_log, "")  << "Free: " << 0 << '\n';
    //    if( !rod1fixed ) TraceStream(g_log, "")  << "Free: " << 1 << '\n';
    //
    //    TraceStream(g_log, "") << "temprelvel: " << posnn[0].dot(m_vnphalf.segment(rodbase,nvdof)) << '\n';
    //  }
    // END TEMP DEBUG STUFF

    // TraceStream(g_log, "") << "Pre-impulse velocity = " << m_vnphalf.segment(rodbase, nvdof) << '\n' << '\n';

    // for (int i = 0; i < nvdof; i += 3)
    // {

    //    double vx = m_vnphalf.segment(rodbase, nvdof)[i];
    //    double vy = m_vnphalf.segment(rodbase, nvdof)[i + 1];
    //   double vz = m_vnphalf.segment(rodbase, nvdof)[i + 2];

    //    double px = m_xn.segment(rodbase, nvdof)[i];
    //    double py = m_xn.segment(rodbase, nvdof)[i + 1];
    //   double pz = m_xn.segment(rodbase, nvdof)[i + 2];

    // TraceStream(g_log, "") << "curve -p " << px << " " << py << " " << pz << " -p " << px + vx << " " << py + vy << " " << pz + vz
    //         << ";" << '\n';
    // }
    // TraceStream(g_log, "") << '\n';

    for (int i = 0; i < numconstraints; ++i)
    {
        //  TraceStream(g_log, "") << "alpha[" << i << "]=" << alpha(i) << " and " << " posnntilde[" << i << "] = " << posnntilde[i]
        //          << '\n' << '\n';

        m_vnphalf.segment(rodbase, nvdof) += alpha(i) * posnntilde[i];
    }

    // TraceStream(g_log, "") << "Post-impulse velocity = " << m_vnphalf.segment(rodbase, nvdof) << '\n' << '\n';

    // for (int i = 0; i < nvdof; i += 3)
    // {
    // double vx = m_vnphalf.segment(rodbase, nvdof)[i];
    // double vy = m_vnphalf.segment(rodbase, nvdof)[i + 1];
    // double vz = m_vnphalf.segment(rodbase, nvdof)[i + 2];

    // double px = m_xn.segment(rodbase, nvdof)[i];
    // double py = m_xn.segment(rodbase, nvdof)[i + 1];
    // double pz = m_xn.segment(rodbase, nvdof)[i + 2];

    // TraceStream(g_log, "") << "curve -p " << px << " " << py << " " << pz << " -p " << px + vx << " " << py + vy << " " << pz + vz
    //         << ";" << '\n';
    //}
    // TraceStream(g_log, "") << '\n';

    applyInextensibilityVelocityFilter(rodidx);

    // ////////////////////////////////////////////////////////
    // // Inextensibility Filter
    // // for my reference: m_xnp1 = m_xn + m_dt * m_vnphalf;
    // //

    // {
    //   Vec3d v0  = m_vnphalf.segment<3> (rodbase+3); // velocity of vertex 1
    //   Vec3d x0  = m_xn.segment<3>      (rodbase+3); // start-of-step position of vertex 1
    //   Vec3d x0N = x0 + m_dt * v0;                   // end-of-step position of vertex 1
    //   for (int i=2; i<m_rods[rodidx]->nv(); ++i)
    //  {
    //    Vec3d x1 = m_xn.segment<3>      (rodbase + 3*i);
    //    Vec3d v1 = m_vnphalf.segment<3> (rodbase + 3*i);

    //    // filter the velocity of vertex i

    //    Vec3d x1N = x1 + m_dt * v1;

    //    double l  = (x1 -x0 ).norm();
    //    double lN = (x1N-x0N).norm();

    //    Vec3d x1Nrevised = x0N + (x1N-x0N)*l/lN;

    //    double lNrevised = (x1Nrevised-x0N).norm();

    //    Vec3d v1revised = (x1Nrevised - x1) / m_dt;

    //    m_vnphalf.segment<3> (rodbase + 3*i) = v1revised;

    //    TraceStream(g_log, "") << "inextensibility: x0N = " << x0N << " x1Nrevised = " << x1Nrevised << " l = " << l << " lNrevised-l = " << (lNrevised-l) << " lNrevised/l = " << (lNrevised/l) << '\n';//x1Nrevised = " << x1Nrevised << " strain = " << (lN/l) << " revised: " << (lNrevised/l) << " l = " << l << " lN = " << lN << " lNrevised = " << lNrevised << " x0N = " << x0N << " x1Nrevised = " << x1Nrevised << " m_vnphalf revised = " << m_vnphalf.segment<3> (rodbase + 3*i) << '\n';

    //    // Vec3d t = (x1-x0).normalized(); // unit tangent

    //    // double relvel = t.dot(v1-v0);

    //    // Vec3d impulse = -relvel*t;

    //    // double relvelAfter = t.dot(v1 + impulse - v0);

    //    // m_vnphalf.segment<3> (rodbase + 3*i) = v1 + impulse;

    //    // TraceStream(g_log, "") << "inextensibility: i=" << i
    //    //        << " x0 = " << x0
    //    //        << " x1 = " << x1
    //    //        << " v0 = " << v0
    //    //        << " v1 = " << v1
    //    //           << " t = " << t
    //    //           << " relvel before = " << relvel
    //    //           << " after = " << relvelAfter << '\n';

    //    x0N = x1Nrevised;
    //    x0  = x1;
    //    v0  = v1;
    //  }
    // }


    ////////////////////////////////////////////////////////


    // TEMP DEBUG STUFF
    //  Vec3d postrelvel = computeRelativeVelocity( m_vnphalf, eecol.e0_v0, eecol.e0_v1, eecol.e1_v0, eecol.e1_v1, eecol.s, eecol.t );
    //  double postmagrelvel = postrelvel.dot(eecol.n);
    //  // Ensure the inelastic impulse decreased the realtive velocity
    //  if( fabs(postmagrelvel) > fabs(magrelvel) )
    //  {
    //    TraceStream(g_log, "") << "Rod: " << rodidx << '\n';
    //    TraceStream(g_log, "") << "Time is: " << m_t << '\n';
    //    TraceStream(g_log, "") << "Incoming relative velocity: " << magrelvel << "      Post-impulse relative velocity: " << postmagrelvel << '\n';
    //
    //    exit(1);
    //  }
    // END TEMP DEBUG STUFF


    // Ensure that the impulse eliminated the realtive velocity
#ifndef NDEBUG
    //Vec3d postrelvel = computeRelativeVelocity( m_vnphalf, eecol.e0_v0, eecol.e0_v1, eecol.e1_v0, eecol.e1_v1, eecol.s, eecol.t );
    //double postmagrelvel = postrelvel.dot(eecol.n);
    // Ensure the inelastic impulse decreased the realtive velocity
    // if( fabs(postmagrelvel) > fabs(magrelvel) )
    // {
    //     TraceStream(g_log, "") << "Rod: " << rodidx << '\n';
    //     TraceStream(g_log, "") << "Time is: " << m_t << '\n';
    //     TraceStream(g_log, "") << "Incoming relative velocity: " << magrelvel << "      Post-impulse relative velocity: " << postmagrelvel << '\n';
    // }
    // assert( fabs(postmagrelvel) <= fabs(magrelvel) );
    // // Should add some 'extra kick' to ensure collision gets killed, but for now just be content with small velocity
    // assert( fabs(postmagrelvel) < 1.0e-9 );
    //assert( postmagrelvel < 0.0 );
#endif

    // Ensure the 'scripted vertices' achieved the desired velocity
#ifndef NDEBUG
    curdof = 1;
    for (size_t i = 0; i < scriptedverts.size(); ++i)
    {
        assert(scriptedverts[i] >= 0);
        assert((int) scriptedverts[i] < m_rods[rodidx]->nv());
        //TraceStream(g_log, "") << m_vnphalf(rodbase+3*scriptedverts[i]+0) << "   " << desired_values(curdof) << '\n';
        assert(approxEq(m_vnphalf(rodbase + 3 * scriptedverts[i] + 0), desired_values(curdof), 1.0e-6));
        ++curdof;
        //TraceStream(g_log, "") << m_vnphalf(rodbase+3*scriptedverts[i]+1) << "   " << desired_values(curdof) << '\n';
        assert(approxEq(m_vnphalf(rodbase + 3 * scriptedverts[i] + 1), desired_values(curdof), 1.0e-6));
        ++curdof;
        //TraceStream(g_log, "") << m_vnphalf(rodbase+3*scriptedverts[i]+2) << "   " << desired_values(curdof) << '\n';
        assert(approxEq(m_vnphalf(rodbase + 3 * scriptedverts[i] + 2), desired_values(curdof), 1.0e-6));
        ++curdof;
    }
#endif

    double postRelativeVelocity = eecol.computeRelativeVelocity();

    //  TraceStream(g_log, "") << "BARodStepper::exertCompliantInelasticEdgeEdgeImpulseOneFixed: Relative velocity before = "
    //         << preRelativeVelocity << " after = " << postRelativeVelocity << '\n';
}

bool BARodStepper::checkExplosions(std::vector<bool>& exploding_rods, const std::vector<bool>& failed_collisions_rods,
        const RodSelectionType& selected_rods)
{
    bool explosions_detected = false;
    double maxRate = 0;
    double maxStart = 0;
    double minStart = std::numeric_limits<double>::max();
    int worstViolator = std::numeric_limits<int>::signaling_NaN();
    TraceStream(g_log, "") << "Checking for explosions\n";

    for (RodSelectionType::const_iterator rod = selected_rods.begin(); rod != selected_rods.end(); rod++)
    {
        if (true || m_steppers[*rod]->HasSolved() && !failed_collisions_rods[*rod])
        {
            if (hadExplosion(*rod))
            {
                explosions_detected = true;
                exploding_rods[*rod] = true;
            }
        }
    }
    if (explosions_detected)
        DebugStream(g_log, "") << "Some rods had explosions\n";

    return explosions_detected;
}

bool BARodStepper::hadExplosion(int rodIdx) const // All forces must have been computed
{
    TraceStream(g_log, "") << "Rod " << rodIdx << " Force norms: initial: " << (*(m_startForces[rodIdx])).norm()
            << " pre-dynamic: " << (*(m_preDynamicForces[rodIdx])).norm() << " pre-collisions: "
            << (*(m_preCollisionForces[rodIdx])).norm() << " post-collisions: " << (*(m_endForces[rodIdx])).norm() << "\n";

    for (int j = 0; j < m_rods[rodIdx]->ndof(); ++j)
    {
        const double s = (*(m_startForces[rodIdx]))[j];
        const double p = (*(m_preCollisionForces[rodIdx]))[j];
        const double e = (*(m_endForces[rodIdx]))[j];
        const double rate = fabs(e - p) / (fabs(p) + m_perf_param.m_explosion_damping);
        if (isnan(rate) || rate > m_perf_param.m_explosion_threshold)
        {
            TraceStream(g_log, "") << "Rod " << rodIdx << " had an explosion on DOF " << j << ": s = " << s << " p = " << p
                    << " e = " << e << " rate = " << rate << " \n";
            return true;
        }
    }
    return false;
}

bool BARodStepper::checkLengths(std::vector<bool>& stretching_rods)
{
    bool stretching_detected = false;
    TraceStream(g_log, "") << "Checking for lengths\n";

    for (RodSelectionType::iterator rod = m_simulated_rods.begin(); rod != m_simulated_rods.end(); rod++)
    {
        stretching_rods[*rod] = !checkLength(*rod);

        if (stretching_rods[*rod])
        {
            // Declare the rod collision-immune for the rest of the time step.
            for (int j = 0; j < m_rods[*rod]->nv(); ++j)
                m_collision_immune[m_base_vtx_indices[*rod] + j] = false;
        }

        stretching_detected = stretching_detected || stretching_rods[*rod];
    }

    if (stretching_detected)
        DebugStream(g_log, "") << "Some rods were stretched by a factor > " << m_perf_param.m_stretching_threshold << '\n';

    return stretching_detected;
}

bool BARodStepper::checkLength(int rodIdx) const
{
    double length = 0.0;

    for (int j = 1; j < m_rods[rodIdx]->nv(); j++)
        length += (m_rods[rodIdx]->getVertex(j) - m_rods[rodIdx]->getVertex(j - 1)).norm();

    if (length > m_initialLengths[rodIdx] * m_perf_param.m_stretching_threshold)
    {
        TraceStream(g_log, "") << "Rod number " << rodIdx << " was stretched by a factor " << length / m_initialLengths[rodIdx]
                << '\n';
        return false;
    }

    return true;
}

void BARodStepper::applyInextensibilityVelocityFilter(int rodidx)
{

    // TraceStream(g_log, "") << "WARNING! SKIPPING INEXTENSIBILITY FILTER!" << '\n';

    // return;


    if (m_level < m_perf_param.m_inextensibility_threshold)
        return;

    int rodbase = m_base_dof_indices[rodidx];

    // if (boundary->isVertexScripted(j))
    //   {
    //  //TraceStream(g_log, "") << "BridsonTimeStepper is calling RodBoundaryCondition at m_t = " << m_t << '\n';
    //  Vec3d desiredposition = boundary->getDesiredVertexPosition(j, m_t);
    //  Vec3d actualvalue = m_xnp1.segment<3> (rodbase + 3 * j);
    //  assert(approxEq(desiredposition, actualvalue, 1.0e-6));
    //   }

    // DebugStream(g_log, "") << "Velocity Filter begin: rod " << rodidx << " vertex 0 = " << m_xn.segment<3>(rodbase) << '\n';

    // m_vnphalf.segment<3> (rodbase + 0) = Vec3d(0, 0, 0);
    // m_vnphalf.segment<3> (rodbase + 3) = Vec3d(0, 0, 0);

    {
        Vec3d v0 = m_vnphalf.segment<3> (rodbase + 3); // velocity of vertex 1
        Vec3d x0 = m_xn.segment<3> (rodbase + 3); // start-of-step position of vertex 1
        Vec3d x0N = x0 + m_dt * v0; // end-of-step position of vertex 1
        for (int i = 2; i < m_rods[rodidx]->nv(); ++i)
        {
            Vec3d x1 = m_xn.segment<3> (rodbase + 3 * i);
            Vec3d v1 = m_vnphalf.segment<3> (rodbase + 3 * i);

            // filter the velocity of vertex i

            Vec3d x1N = x1 + m_dt * v1;

            double l = (x1 - x0).norm();
            double lN = (x1N - x0N).norm();

            Vec3d x1Nrevised = x0N + (x1N - x0N) * l / lN;

            double lNrevised = (x1Nrevised - x0N).norm();

            Vec3d v1revised = (x1Nrevised - x1) / m_dt;

            m_vnphalf.segment<3> (rodbase + 3 * i) = v1revised;

            //TraceStream(g_log, "") << "inextensibility: x0N = " << x0N << " x1N = " << x1N << " x1Nrevised = " << x1Nrevised << " l = " << l << " lNrevised = "
            //            << lNrevised << " l = " << l << '\n';

            //x1Nrevised = " << x1Nrevised << " strain = " << (lN/l) << " revised: " << (lNrevised/l) << " l = " << l << " lN = " << lN << " lNrevised = " << lNrevised << " x0N = " << x0N << " x1Nrevised = " << x1Nrevised << " m_vnphalf revised = " << m_vnphalf.segment<3> (rodbase + 3*i) << '\n';

            // Vec3d t = (x1-x0).normalized(); // unit tangent

            // double relvel = t.dot(v1-v0);

            // Vec3d impulse = -relvel*t;

            // double relvelAfter = t.dot(v1 + impulse - v0);

            // m_vnphalf.segment<3> (rodbase + 3*i) = v1 + impulse;

            // TraceStream(g_log, "") << "inextensibility: i=" << i
            //      << " x0 = " << x0
            //      << " x1 = " << x1
            //      << " v0 = " << v0
            //      << " v1 = " << v1
            //           << " t = " << t
            //           << " relvel before = " << relvel
            //           << " after = " << relvelAfter << '\n';

            x0N = x1Nrevised;
            x0 = x1;
            v0 = v1;
        }
    }

    // TraceStream(g_log, "") << "Edge lengths after inextensibility: ";
    // for (int i = 0; i < m_rods[rodidx]->nv() - 1; ++i)
    // {
    //     Vec3d x0 = m_xn.segment<3> (rodbase + 3 * i);
    //     Vec3d v0 = m_vnphalf.segment<3> (rodbase + 3 * i);
    //     Vec3d x1 = m_xn.segment<3> (rodbase + 3 * i + 3);
    //     Vec3d v1 = m_vnphalf.segment<3> (rodbase + 3 * i + 3);
    //  Vec3d x0N = x0 + m_dt * v0;
    //  Vec3d x1N = x1 + m_dt * v1;
    //  Vec3d eN = (x1N - x0N);
    //  TraceStream(g_log, "") << " " << eN.norm();
    // }
    // TraceStream(g_log, "") << '\n';

    // DebugStream(g_log, "") << "Velocity Filter end: rod " << rodidx << " vertex 0 = " << m_xn.segment<3>(rodbase) << '\n';


}

/**
 * Implicit penalty response
 */
void BARodStepper::setupPenaltyForces(std::list<Collision*>& collisions, const RodSelectionType& selected_rods)
{
    // Detect proximity collisions
    m_collision_detector->getCollisions(collisions, Proximity);

    for (int rod_id = 0; rod_id < m_number_of_rods; rod_id++)
        assert(m_implicit_pnlty_forces[rod_id]->cleared());

    // Store the proximity collision in the RodPenaltyForce
    for (std::list<Collision*>::const_iterator col = collisions.begin(); col != collisions.end(); col++)
    {
        VertexFaceProximityCollision* vfpcol = dynamic_cast<VertexFaceProximityCollision*> (*col);
        if (vfpcol)
        {
            assert(vfpcol->isAnalysed());
            int rod_id = getContainingRod(vfpcol->v0);
            TraceStream(g_log, "") << "Creating penalty force for rod " << rod_id << " address " << &m_rods[rod_id] << '\n';
            int v_id = vfpcol->v0 - m_base_dof_indices[rod_id] / 3;
            m_implicit_pnlty_forces[rod_id]->registerProximityCollision(v_id, vfpcol);
        }
        // TODO: delete vfpcol once used. Not here though, as vfpcol is used each time RodPenaltyForce::computeForce is called
    }

}

void BARodStepper::setImplicitPenaltyExtraThickness(const double& h)
{
    m_perf_param.m_implicit_thickness = h;
}

void BARodStepper::setVertexFacePenalty(const double& k)
{
    assert(k >= 0.0);
    m_perf_param.m_implicit_stiffness = k;
}

}
