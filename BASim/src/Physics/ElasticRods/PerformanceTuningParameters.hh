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

struct FailureMode
{
    enum ResponseSeverity
    {
        IgnoreError, KillTheRod, HaltSimulation
    };

    int m_max_iterations; // For solver and collision only I suppose
    int m_max_substeps;
    ResponseSeverity m_in_case_of;

    FailureMode(int max_iterations, int max_substeps, ResponseSeverity in_case_of) :
        m_max_iterations(max_iterations), m_max_substeps(max_substeps), m_in_case_of(in_case_of)
    {
    }
};

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
    // Explosion detection
    bool m_enable_explosion_detection;
    // Damping parameter in explosion detection
    double m_explosion_damping;
    // Threshold in explosion detection
    double m_explosion_threshold;
    // Stretching detection
    double m_stretching_threshold;

    /**
     * Failure modes
     */
    FailureMode m_solver;
    FailureMode m_collision;
    FailureMode m_explosion;
    FailureMode m_stretching;

    PerformanceTuningParameters();

};

}

#endif /* PERFORMANCETUNINGPARAMETERS_HH_ */
