/*
 * PerformanceTuningParameters.cc
 *
 *  Created on: 19/05/2011
 *      Author: jaubry
 */

#include "PerformanceTuningParameters.hh"

namespace BASim
{

PerformanceTuningParameters::PerformanceTuningParameters() :
    m_enable_penalty_response(true), //
            m_implicit_thickness(1.0), //
            m_implicit_stiffness(200.0), //
            m_skipRodRodCollisions(true), //
            m_inextensibility_threshold(3), //
            m_maximum_number_of_solver_iterations(50), //
            m_maximum_number_of_collisions_iterations(10), //
            m_enable_explosion_detection(true), //
            m_explosion_damping(100.0), //
            m_explosion_threshold(100.0), //
            m_in_case_of_solver_failure(KillTheRod), //
            m_in_case_of_collision_failure(KillTheRod), //
            m_in_case_of_explosion_failure(KillTheRod), //
            m_max_number_of_substeps_for_solver(0), //
            m_max_number_of_substeps_for_collision(0), //
            m_max_number_of_substeps_for_explosion(0)
{
}

}
