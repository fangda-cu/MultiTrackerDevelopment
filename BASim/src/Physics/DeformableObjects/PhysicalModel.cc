#include "BASim/src/Physics/DeformableObjects/PhysicalModel.hh"
#include "BASim/src/Physics/DeformableObjects/DeformableObject.hh"
namespace BASim {


PhysicalModel::PhysicalModel(BASim::DeformableObject &obj):
  m_obj(obj), m_vertexDofIdxs(&obj), m_edgeDofIdxs(&obj), m_faceDofIdxs(&obj), m_tetDofIdxs(&obj) {}


}
