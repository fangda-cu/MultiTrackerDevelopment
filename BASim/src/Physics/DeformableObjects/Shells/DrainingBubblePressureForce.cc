#include "BASim/src/Physics/DeformableObjects/Shells/DrainingBubblePressureForce.hh"
#include "BASim/src/Physics/DeformableObjects/DeformableObject.hh"
#include "BASim/src/Core/TopologicalObject/TopObjIterators.hh"
#include "BASim/src/Math/Math.hh"

namespace BASim {

  DrainingBubblePressureForce::DrainingBubblePressureForce( ElasticShell& shell, const std::string& name, std::vector<EdgeHandle> holeEdges, std::vector<EdgeHandle> baseEdges, Scalar gas_density, Scalar time_step ) : 
    ElasticShellForce(shell, name), m_hole_edges(holeEdges), m_base_edges(baseEdges), m_gas_density(gas_density), m_timestep(time_step)
{
  update();
}

Scalar DrainingBubblePressureForce::globalEnergy() const
{
  return 0;
}

void DrainingBubblePressureForce::globalForce( VecXd& force ) const
{
  //Compute the normals
  VertexProperty<Vec3d> normals(&m_shell.getDefoObj());
  m_shell.getVertexNormals(normals);

  //Add normal pressure force to all nodes
  DeformableObject& obj = m_shell.getDefoObj();
  for(VertexIterator vit = obj.vertices_begin(); vit != obj.vertices_end(); ++vit) {
    VertexHandle& vh = *vit;
    int dofIdx = m_shell.getDefoObj().getPositionDofBase(vh);//getVertexDofBase(vh);
    Vec3d curNormal = normals[vh];
    
    force[dofIdx]   += m_current_pressure * curNormal[0];
    force[dofIdx+1] += m_current_pressure * curNormal[1];
    force[dofIdx+2] += m_current_pressure * curNormal[2];
  }
  
}

void DrainingBubblePressureForce::globalJacobian( Scalar scale, MatrixBase& Jacobian ) const
{
  //Assume explicit integration of the pressure force, therefore no Jacobian
  return;
}

void DrainingBubblePressureForce::update() {
  

  //Compute current volume of the bubble
  Scalar currentVolume = 0; 

  //We'll assume the top and bottom are closed with triangle fans, and use the usual method to integrate the volume of a polyhedron
  DeformableObject& obj = m_shell.getDefoObj();
  
  //Choose a reference position arbitrarily (first vertex in the mesh)
  Vec3d refPoint = m_shell.getVertexPosition(*(obj.vertices_begin()));
  for(FaceIterator fit = m_shell.getDefoObj().faces_begin(); fit != m_shell.getDefoObj().faces_end();++fit) {
    FaceHandle fh = *fit;
    FaceVertexIterator fvit = obj.fv_iter(fh);
    Vec3d v0 = m_shell.getVertexPosition(*fvit); ++fvit;
    Vec3d v1 = m_shell.getVertexPosition(*fvit); ++fvit;
    Vec3d v2 = m_shell.getVertexPosition(*fvit);
    currentVolume += (v0-refPoint).dot((v1-refPoint).cross(v2-refPoint))/6.0; //add the tet volume
  }
  
  //Approximate the area of the hole with a simple triangle fan
  Scalar holeArea = 0; 
  if(m_hole_edges.size() < 2) return; //not a hole

  //Get the base vertex for the fan
  EdgeHandle edge0 = m_hole_edges[0];
  VertexHandle vh0 = m_shell.getDefoObj().fromVertex(edge0);
  Vec3d v0 = m_shell.getVertexPosition(vh0);

  //Iterate over all the edges, summing triangle areas, and the top loop volume
  for(unsigned int i = 0; i < m_hole_edges.size(); ++i) {
    EdgeHandle curEdge = m_hole_edges[i];
    VertexHandle vh1 = m_shell.getDefoObj().fromVertex(curEdge);
    VertexHandle vh2 = m_shell.getDefoObj().toVertex(curEdge);
    
    if(vh1 == vh0 || vh2 == vh0) continue; //skip edges bordering the base vertex, since they don't form a triangle

    Vec3d v1 = m_shell.getVertexPosition(vh1);
    Vec3d v2 = m_shell.getVertexPosition(vh2);
    holeArea += 0.5f*((v2-v0).cross(v1-v0)).norm();

    //Add "volume" of top hole loop (is this oriented properly)
    currentVolume += (v0-refPoint).dot((v1-refPoint).cross(v2-refPoint))/6.0; //add the tet volume
  }

  //Now do the base circle of the bubble
  if(m_base_edges.size() < 2) return; //not a hole

  //Get the base vertex for the fan
  edge0 = m_base_edges[0];
  vh0 = m_shell.getDefoObj().fromVertex(edge0);
  v0 = m_shell.getVertexPosition(vh0);

  //Iterate over all the edges, summing triangle areas, and the top loop volume
  for(unsigned int i = 0; i < m_base_edges.size(); ++i) {
    EdgeHandle curEdge = m_base_edges[i];
    VertexHandle vh1 = m_shell.getDefoObj().fromVertex(curEdge);
    VertexHandle vh2 = m_shell.getDefoObj().toVertex(curEdge);

    if(vh1 == vh0 || vh2 == vh0) continue; //skip edges bordering the base vertex, since they don't form a triangle

    Vec3d v1 = m_shell.getVertexPosition(vh1);
    Vec3d v2 = m_shell.getVertexPosition(vh2);
    
    //Add "volume" of top hole loop (is this oriented properly)
    currentVolume += (v0-refPoint).dot((v1-refPoint).cross(v2-refPoint))/6.0; //add the tet volume
  }

  //Approximate the gas velocity as (rho*v^2)/2 (Bernoulli)
  Scalar volumeChange = currentVolume - m_old_volume; 
  Scalar velocity = volumeChange / m_timestep / holeArea;
  
  //Store the pressure to add as a force
  m_current_pressure = 0.5* m_gas_density * square(velocity);
  
  //Update the volume for next time
  m_old_volume = currentVolume;
}

} //namespace BASim