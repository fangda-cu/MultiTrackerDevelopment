//
//  PositionDofsModel.cc
//  BASim
//
//  Created by Fang Da (fang@cs.columbia.edu) on 5/11/12.
//  Copyright (c) 2012 Columbia. All rights reserved.
//

#include "BASim/src/Physics/DeformableObjects/PositionDofsModel.hh"

namespace BASim
{
  void PositionDofsModel::startStep(Scalar time, Scalar timestep)
  {
    m_damping_undeformed_positions = m_positions;
  }

  void PositionDofsModel::endStep(Scalar time, Scalar timestep)
  {
    
  }
  
}
