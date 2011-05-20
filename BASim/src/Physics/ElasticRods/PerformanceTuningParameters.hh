/*
 * PerformanceTuningParameters.hh
 *
 *  Created on: 19/05/2011
 *      Author: jaubry
 */

#ifndef PERFORMANCETUNINGPARAMETERS_HH_
#define PERFORMANCETUNINGPARAMETERS_HH_

namespace BASim
{

struct PerformanceTuningParameters
{
    /**
     * General parameters
     */

    /**
     * Implicit rod/mesh penalty
     */
    // Whether we want to apply implicit penalty response for rod/mesh collisions
    bool m_enable_penalty_response;
    // Extra thickness in rod/mesh implicit penalty response
    double m_implicit_thickness;
    // Penalty stiffness in rod/mesh implicit penalty response
    double m_implicit_stiffness;

    /**
     * Rod-rod collisions
     */
    // Whether we should ignore rod-rod collisions
    bool m_skipRodRodCollisions;

    /**
     * Inextensibility filter
     */
    // Number of times the original step has to be halved before the inextensibility filter is applied
    int m_inextensibility_threshold;

    /**
     * Failure detection
     */

    /**
     * Solver
     */
    // Maximum number of iterations allowed in the solver
    int m_maximum_number_of_solver_iterations;

    /**
     * Collisions
     */
    // Maximum number of iterations allowed in collision response. Set to 0 to disable collision response, to 1 to disable iterations.
    int m_maximum_number_of_collisions_iterations;

    /**
     * Explosions
     */
    // Explosion detection
    bool m_enable_explosion_detection;
    // Damping parameter in explosion detection
    double m_explosion_damping;
    // Threshold in explosion detection
    double m_explosion_threshold;

    /**
     * Response to failure
     */
    enum ResponseSeverity
    {
        IgnoreError, KillTheRod, HaltSimulation
    };

    ResponseSeverity m_in_case_of_solver_failure;
    ResponseSeverity m_in_case_of_collision_failure;
    ResponseSeverity m_in_case_of_explosion_failure;

    /**
     * Substepping: if in KillTheRod mode, how deep are we going before shedding.
     */
    // Maximum number of times the original step will be halved by the adaptive substepping for solver-failing rods
    int m_max_number_of_substeps_for_solver;
    // Maximum number of times the original step will be halved by the adaptive substepping for collision-failing rods
    int m_max_number_of_substeps_for_collision;
    // Maximum number of times the original step will be halved by the adaptive substepping for explosion-failing rods
    int m_max_number_of_substeps_for_explosion;

    PerformanceTuningParameters();

};

}

#endif /* PERFORMANCETUNINGPARAMETERS_HH_ */
