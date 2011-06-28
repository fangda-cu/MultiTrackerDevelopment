/*
 * PerformanceTuningParameters.hh
 *
 *  Created on: 19/05/2011
 *      Author: jaubry
 */

#ifndef PERFORMANCETUNINGPARAMETERS_HH_
#define PERFORMANCETUNINGPARAMETERS_HH_

#include <string>
#include <boost/lexical_cast.hpp>

namespace BASim
{

struct FailureMode
{
    enum ResponseSeverity
    {
        IgnoreError, KillTheRod, HaltSimulation
    };

    // Name of the failure
    std::string m_name;
    // Maximum number of iterations when applicable (solver and collision)
    int m_max_iterations;
    // Maximum level of binary substepping that this failure mode will trigger
    int m_max_substeps;
    // What to do in case of failure past the max level of substepping
    ResponseSeverity m_in_case_of;
    // Counter of rods killed because of this failure mode during the time step
    int m_num_killed;
    // Cumulative counter of rods killed because of this failure mode
    int m_total_killed;

    FailureMode(std::string name, int max_iterations, int max_substeps, ResponseSeverity in_case_of) :
        m_name(name), m_max_iterations(max_iterations), m_max_substeps(max_substeps), m_in_case_of(in_case_of),
                m_num_killed(0), m_total_killed(0)
    {
    }

    void resetNum()
    {
        m_num_killed = 0;
    }

    std::string sumMessage()
    {
        return "Rods killed because of " + m_name + " failure: " + boost::lexical_cast<std::string>(m_num_killed)
                + " (this step), " + boost::lexical_cast<std::string>(m_total_killed) + " (total)";
    }

    FailureMode& operator++() // Prefix operator
    {
        ++m_num_killed;
        ++m_total_killed;

        return *this;
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
    // Subsampling rate for level set detection / response in RodLevelSetForce
    int m_levelset_subsampling;

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
