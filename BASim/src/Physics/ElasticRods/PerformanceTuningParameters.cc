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
            m_levelset_subsampling(0), //
            m_skipRodRodCollisions(true), //
            m_inextensibility_threshold(3), //
            m_enable_explosion_detection(true), //
            m_explosion_damping(100.0), //
            m_explosion_threshold(100.0), //
            m_stretching_threshold(2.0), //
            m_solver("solver", 50, 0, FailureMode::KillTheRod), //
            m_collision("collision", 10, 0, FailureMode::KillTheRod), //
            m_explosion("explosion", 0, 0, FailureMode::KillTheRod), //
            m_stretching("stretching", 0, 0, FailureMode::KillTheRod)
{
}

}
