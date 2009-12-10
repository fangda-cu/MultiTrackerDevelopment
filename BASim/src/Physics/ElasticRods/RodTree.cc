#include "RodTree.hh"
#include <BASim/Render>
#include <Eigen/Geometry>

namespace BASim {

RodTree::RodTree(ElasticRod* parentRod, ElasticRod *rod, Vec3d rootPoint)
{
  m_parentRod = parentRod;
  if (m_parentRod == NULL)
    m_rootPoint = rootPoint;
  m_rod = rod;
  m_childNodes.clear();

  getRodMaterialFrame( m_rod, 0, m_material1, m_material2, m_tangent );

  Vec3d pMaterial1, pMaterial2, pTangent;
  if ( m_parentRod != NULL )
    getRodMaterialFrame( m_parentRod, m_parentRod->ne()-1, pMaterial1, pMaterial2, pTangent );
  else
    getRodMaterialFrame( NULL, 0, pMaterial1, pMaterial2, pTangent );

  m_material1.normalize();
  m_material2.normalize();
  pMaterial1.normalize();
  pMaterial2.normalize();

  Quaternion mat1Rot;
  //mat1Rot.setFromTwoVectors( pMaterial1, m_material1 );
  double angle = acos(pMaterial1.dot(m_material1));
  Vec3d axis = pMaterial1.cross(m_material1);
  axis.normalize();
  mat1Rot = Quaternion(Eigen::AngleAxis<Scalar>(angle, axis));
  Vec3d rotatedMaterial2 = mat1Rot * pMaterial2;

  Quaternion mat2Rot;
  //mat1Rot.setFromTwoVectors( rotatedMaterial2, m_material2 );
  angle = acos(rotatedMaterial2.dot(m_material2));
  axis = m_material1;
  axis.normalize();
  mat2Rot = Quaternion(Eigen::AngleAxis<Scalar>(angle, axis));
  

  /*angle = acos(pMaterial2.dot(m_material2));
  axis = pMaterial2.cross(m_material2);;
  axis.normalize();
  mat2Rot = Quaternion(Eigen::AngleAxis<Scalar>(angle, axis));
  */

  m_rotationFromParentFrame = mat2Rot * mat1Rot;
  //m_rotationFromParentFrame = mat2Rot;

  // Find the Quaternion that will rotate the parent's end material frame 
  // to align with the start edge material frame.
  /*double angle = acos(pMaterial1.dot(m_material1));
  Vec3d axis = pMaterial1.cross(m_material1);
  axis.normalize();
  m_material1Rot = Quaternion(Eigen::AngleAxis<Scalar>(angle, axis));

  angle = acos(pMaterial2.dot(m_material2));
  axis = pMaterial2.cross(m_material2);
  axis.normalize();
  m_material2Rot = Quaternion(Eigen::AngleAxis<Scalar>(angle, axis));

  m_rotationFromParentFrame = m_material1Rot * m_material2Rot;*/
}

RodTree::~RodTree() 
{
}

void RodTree::addChildNode(RodTree* child)
{
  m_childNodes.push_back(child);
}
  
vector<RodTree*>& RodTree::getChildNodes()
{
  return m_childNodes;
}

void RodTree::enforceParentConnections()
{
  // First move all the children to sit at the end of the parent
  Vec3d parentEndPosition;
  if (m_parentRod == NULL )
  {
    // We are the root so make sure we are stuck where we should be
    parentEndPosition = m_rootPoint;
  }
  else 
  {
    parentEndPosition = m_parentRod->getVertex(m_parentRod->nv()-1);  
  }

  Vec3d delta = parentEndPosition - m_rod->getVertex(0);
  for (int v=0; v<m_rod->nv(); ++v)
  {
    Vec3d oldPosition = m_rod->getVertex(v);
    m_rod->setVertex(v, oldPosition+delta);
  }

  for (size_t c=0; c<m_childNodes.size(); ++c)
  {
    m_childNodes[c]->enforceParentConnections();
  }
}

void RodTree::getRodMaterialFrame( ElasticRod* rod, const int edge, Vec3d& m1, Vec3d& m2, Vec3d&t  )
{
  if ( rod != NULL )
  {
    m1 = rod->getMaterial1( edge );
    m2 = rod->getMaterial2( edge );
    t = rod->getTangent( edge );
  }
  else
  {
    // They passed in no rod so they must be looking for a parent frame
    // but there is no parent.
    t =  Vec3d(1,0,0);
    m1 = Vec3d(0,1,0);
    m2 = Vec3d(0,0,1);
  }
}

void RodTree::draw()
{
  Vec3d material1;
  Vec3d material2;
  Vec3d tangent;

  // Start of rod frame
  getRodMaterialFrame( m_rod, m_rod->ne()-1, material1, material2, tangent );

  glLineWidth( 3.0 );
  glBegin( GL_LINES );
  glColor3f(0,1,0);
  Vec3d vertStart = ( m_rod->getVertex( m_rod->nv()-2 ) + m_rod->getVertex( m_rod->nv()-1 ) ) / 2.0;
  glVertex3f( vertStart[0], vertStart[1], vertStart[2] );
  Vec3d vert = vertStart+material1;
  glVertex3f( vert[0], vert[1], vert[2] );
  glColor3f(0,0,1);
  glVertex3f( vertStart[0], vertStart[1], vertStart[2] );
  vert = vertStart+material2;
  glVertex3f( vert[0], vert[1], vert[2] );
  glEnd();

  // End of rod frame
  getRodMaterialFrame( m_rod, 0, material1, material2, tangent );
  
  glBegin( GL_LINES );
  glColor3f(0,1,0);
  vertStart = ( m_rod->getVertex( 0 ) + m_rod->getVertex( 1 ) ) / 2.0;
  glVertex3f( vertStart[0], vertStart[1], vertStart[2] );
  vert = vertStart+material1;
  glVertex3f( vert[0], vert[1], vert[2] );
  glColor3f(0,0,1);
  glVertex3f( vertStart[0], vertStart[1], vertStart[2] );
  vert = vertStart+material2;
  glVertex3f( vert[0], vert[1], vert[2] );
  glEnd();

  // Start of rod rotated parent frame
  Vec3d pMaterial1, pMaterial2, pTangent;
  
  if ( m_parentRod != NULL )
  {
    getRodMaterialFrame( m_parentRod, m_parentRod->ne()-1, pMaterial1, pMaterial2, pTangent );
    
    material1 = m_rotationFromParentFrame * pMaterial1;
    material2 = m_rotationFromParentFrame * pMaterial2;
   // tangent = rotationFromParentFrame * pMaterial1;
    
    glBegin( GL_LINES );
    glColor3f(1,0,0);
    vertStart = ( m_rod->getVertex( 1 ) + m_rod->getVertex( 0 ) ) / 2.0;
    glVertex3f( vertStart[0], vertStart[1], vertStart[2] );
    vert = vertStart+material1 * 2;
    glVertex3f( vert[0], vert[1], vert[2] );
    glVertex3f( vertStart[0], vertStart[1], vertStart[2] );
    vert = vertStart+material2 * 2;
    glVertex3f( vert[0], vert[1], vert[2] );
    glEnd();
  }

  glLineWidth( 1.0 );

  for (size_t c=0; c<m_childNodes.size(); ++c)
  {
    m_childNodes[c]->draw();
  }  
}

}
