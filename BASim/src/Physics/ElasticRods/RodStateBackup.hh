/**
 * \file RodStateBackup.hh
 *
 * \author smith@cs.columbia.edu
 * \date 06/28/2010
 */


#ifndef RODSTATEBACKUP_HH
#define RODSTATEBACKUP_HH

#include <fstream>

#ifdef WETA
#include "ElasticRod.hh"
#include "../../Core/Property.hh"
#include "RodStretchingForce.hh"
#include "../../IO/SerializationUtils.hh"
#else
#include "BASim/src/Physics/ElasticRods/ElasticRod.hh"
#include "BASim/src/Core/Property.hh"
#include "BASim/src/Physics/ElasticRods/RodStretchingForce.hh"
#include "BASim/src/IO/SerializationUtils.hh"
#endif

// TODO: Add some more specific tests when comparing rod state... 
//       Properties might not be 'attatched'

namespace BASim 
{

class RodState
{
public:

  RodState();

  ~RodState();

  void copyState( ElasticRod& rod );

  // Restores state, assuming the EXACT same properties are present.
  void restoreState( ElasticRod& rod );
  
  // Clears the state of the rod, sets entirely new state
  void setState( ElasticRod& rod );

  // A debug check to ensure no properties of the same name appear.
  void ensureNoRepeatProperties();
  
  void clear();
  
  void print( ElasticRod& rod );

  bool compareProperties( const ElasticRod& other, bool printdiffs );

private:

  const PropertyContainer& getVertexPropertyContainer( const ElasticRod& rod );

  PropertyContainer& getVertexPropertyContainer( ElasticRod& rod );

  const PropertyContainer& getEdgePropertyContainer( const ElasticRod& rod );

  PropertyContainer& getEdgePropertyContainer( ElasticRod& rod );

  const PropertyContainer& getFacePropertyContainer( const ElasticRod& rod );
  
  PropertyContainer& getFacePropertyContainer( ElasticRod& rod );
  
  const PropertyContainer& getObjectPropertyContainer( const ElasticRod& rod );
  
  PropertyContainer& getObjectPropertyContainer( ElasticRod& rod );

  ///////////////////////////////////////////////////////////////////////////
  // Generic implementation for properties that do not need special treatment
  //template<class T>
  //void printProperty( Property<T>* prp )
  //{
  //  for( int i  = 0; i < (int) prp->size(); ++i ) std::cout << (*prp)[i] << " ";
  //}

  //////////////////////////
  // INT STUFF
  void printProperty( Property<int>* prp )
  {
    for( int i  = 0; i < (int) prp->size(); ++i )
    {
      std::cout << (*prp)[i] << " ";
    }
  }  
  
  bool propertiesEqual( Property<int>* prpA, Property<int>* prpB )
  {
    if( prpA->size() != prpB->size() ) return false;
    for( int i  = 0; i < (int) prpA->size(); ++i ) if( (*prpA)[i] != (*prpB)[i] ) return false;
    return true;
  }
  
  
  //////////////////////////
  // BOOL STUFF
  void printProperty( Property<bool>* prp )
  {
    for( int i  = 0; i < (int) prp->size(); ++i )
    {
      std::cout << (*prp)[i] << " ";
    }
  }
  
  bool propertiesEqual( Property<bool>* prpA, Property<bool>* prpB )
  {
    if( prpA->size() != prpB->size() ) return false;
    for( int i  = 0; i < (int) prpA->size(); ++i ) if( (*prpA)[i] != (*prpB)[i] ) return false;
    return true;
  }


  //////////////////////////
  // SCALAR STUFF
//  void printProperty( Property<Scalar>* prp )
//  {
//    for( int i  = 0; i < (int) prp->size(); ++i )
//    {
//      std::cout << (*prp)[i] << " ";
//    }
//  }
//  
//  bool propertiesEqual( Property<Scalar>* prpA, Property<Scalar>* prpB, double eps = 1.0e-6 )
//  {
//    if( prpA->size() != prpB->size() ) return false;
//    for( int i  = 0; i < (int) prpA->size(); ++i ) if( !approxEq((*prpA)[i],(*prpB)[i],eps) ) return false;
//    return true;
//  }
  
  //////////////////////////
  // REFRAME STUFF
  void printProperty( Property<ElasticRod::RefFrameType>* prp )
  {
    for( int i  = 0; i < (int) prp->size(); ++i )
    {
      switch( (*prp)[i] ) 
      {
        case ElasticRod::SpaceParallel:
        {
          std::cout << "SpaceParallel" << " ";
          break;
        }
        case ElasticRod::TimeParallel:
        {
          std::cout << "TimeParallel" << " ";
          break;
        }
        default:
        {
          std::cout << "wow, unknown frame type. i'm just a gonna exit here" << std::endl;
          exit(1);
          break;
        }
      }
    }
  }
  
  bool propertiesEqual( Property<ElasticRod::RefFrameType>* prpA, Property<ElasticRod::RefFrameType>* prpB )
  {
    if( prpA->size() != prpB->size() ) return false;
    for( int i  = 0; i < (int) prpA->size(); ++i ) if( (*prpA)[i] != (*prpB)[i] ) return false;
    return true;
  }
  
  ////////////////////////////
  // RODBOUNDARYCONDITIONSTUFF
  void printProperty( Property<RodBoundaryCondition::BCList>* prp )
  {
    for( int i = 0; i < (int) prp->size(); ++i )
    {
      for( int j = 0; j < (int) (*prp)[i].size(); ++j )
      {
        std::cout << (*prp)[i][j] << " ";
      }
    }
  }

  bool propertiesEqual( Property<RodBoundaryCondition::BCList>* prpA, Property<RodBoundaryCondition::BCList>* prpB )
  {
    if( prpA->size() != prpB->size() ) return false;
    for( int i = 0; i < (int) prpA->size(); ++i ) if( (*prpA)[i].size() != (*prpB)[i].size() ) return false;
    for( int i = 0; i < (int) prpA->size(); ++i ) for( int j = 0; j < (int) (*prpA)[i].size(); ++j ) if( (*prpA)[i][j] != (*prpB)[i][j] ) return false;
    return true;
  }

  ////////////////////////////
  // FORCE STUFF
  void printProperty( Property<ElasticRod::RodForces>* prp )
  {
    for( int i  = 0; i < (int) prp->size(); ++i )
    {
      for( int j = 0; j < (int) (*prp)[i].size(); ++j )
      {
        std::cout << (*prp)[i][j]->getName() << " ";
      }
    }
  }

  bool propertiesEqual( Property<ElasticRod::RodForces>* prpA, Property<ElasticRod::RodForces>* prpB )
  {
    if( prpA->size() != prpB->size() ) return false;
    for( int i = 0; i < (int) prpA->size(); ++i ) for( int j = 0; j < (int) (*prpA)[i].size(); ++j ) if( (*prpA)[i][j]->getName() != (*prpB)[i][j]->getName() ) return false;
    return true;
  }  

  ////////////////////////////
  // DOFMAP STUFF
  // Used by: void printProperty( Property<DOFMap>* prp )
  void printDofHandle( const DofHandle& dh )
  {
    std::cout << "<";
    
    switch(dh.getType()) 
    {
      case DofHandle::VERTEX_DOF:
      {
        std::cout << "VERTEX_DOF";
        break;
      }
      case DofHandle::EDGE_DOF:
      {
        std::cout << "EDGE_DOF";
        break;
      }
      default:
      {
        std::cerr << "blablabla...exiting cause i dont recognize dofhandle" << std::endl;
        exit(1);
        break;
      }
    }
    
    std::cout << "," << dh.getNum();
    std::cout << "," << dh.idx();
    std::cout << "," << dh.getHandle().idx();
    
    std::cout << ">";
  }
  
  void printProperty( Property<DOFMap>* prp )
  {
    for( int i = 0; i < (int) prp->size(); ++i )
    {
      std::cout << "DOF_TO_HANDLE: ";
      const std::map<DofHandle,int>& dofToIndex = (*prp)[i].getDofToIndexMap();
      for( std::map<DofHandle,int>::const_iterator itr = dofToIndex.begin(); itr != dofToIndex.end(); ++itr )
      {
        printDofHandle( itr->first );
        std::cout << "->" << itr->second;
        std::cout << "   ";
      }
      std::cout << "        ";
      std::cout << "HANDLE_TO_DOF: ";
      const std::map<int,DofHandle>& indexToDof = (*prp)[i].getIndexToDofMap();
      for( std::map<int,DofHandle>::const_iterator itr = indexToDof.begin(); itr != indexToDof.end(); ++itr )
      {
        std::cout << itr->first << "->";
        printDofHandle( itr->second );
        std::cout << "   ";
      }      
    }
  }

  bool propertiesEqual( Property<DOFMap>* prpA, Property<DOFMap>* prpB )
  {    
    if( prpA->size() != prpB->size() ) return false;
    for( int i = 0; i < (int) prpA->size(); ++i )
    {
      const std::map<DofHandle,int>& dofToIndexA = (*prpA)[i].getDofToIndexMap();
      const std::map<DofHandle,int>& dofToIndexB = (*prpB)[i].getDofToIndexMap();
      if( dofToIndexA.size() != dofToIndexB.size() ) return false;
      for( std::map<DofHandle,int>::const_iterator itrA = dofToIndexA.begin(), itrB = dofToIndexB.begin(); itrA != dofToIndexA.end(); ++itrA, ++itrB )
      {
        if( itrA->first != itrB->first )   return false;
        if( itrA->second != itrB->second ) return false;
      }
      
      const std::map<int,DofHandle>& indexToDofA = (*prpA)[i].getIndexToDofMap();
      const std::map<int,DofHandle>& indexToDofB = (*prpB)[i].getIndexToDofMap();
      if( indexToDofA.size() != indexToDofB.size() ) return false;
      for( std::map<int,DofHandle>::const_iterator itrA = indexToDofA.begin(), itrB = indexToDofB.begin(); itrA != indexToDofA.end(); ++itrA, ++itrB )
      {
        if( itrA->first != itrB->first )   return false;
        if( itrA->second != itrB->second ) return false;
      }      
    }
    return true;
  }  
  
  ///////////////////////
  // VEC2D STUFF
  void printProperty( Property<Vec2d>* prp )
  {
    for( int i  = 0; i < (int) prp->size(); ++i )
    {
      std::cout << (*prp)[i].transpose() << " ";
    }
  }
  
  bool propertiesEqual( Property<Vec2d>* prpA, Property<Vec2d>* prpB, double eps = 1.0e-6 )
  {
    if( prpA->size() != prpB->size() ) return false;
    for( int i = 0; i < (int) prpA->size(); ++i ) if( !approxEq((*prpA)[i],(*prpB)[i],eps) ) return false;
    return true;    
  }

  ///////////////////////
  // VEC3D STUFF
//  void printProperty( Property<Vec3d>* prp )
//  {
//    for( int i  = 0; i < (int) prp->size(); ++i )
//    {
//      std::cout << (*prp)[i].transpose() << " ";
//    }
//  }
//
//  bool propertiesEqual( Property<Vec3d>* prpA, Property<Vec3d>* prpB, double eps = 1.0e-6 )
//  {
//    if( prpA->size() != prpB->size() ) return false;
//    for( int i = 0; i < (int) prpA->size(); ++i ) if( !approxEq((*prpA)[i],(*prpB)[i],eps) ) return false;
//    return true;    
//  }

  ///////////////////////
  // VECXD STUFF
  void printProperty( Property<VecXd>* prp )
  {
    for( int i  = 0; i < (int) prp->size(); ++i )
    {
      std::cout << (*prp)[i].size() << "-" << (*prp)[i].transpose() << " ";
    }
  }

  bool propertiesEqual( Property<VecXd>* prpA, Property<VecXd>* prpB, double eps = 1.0e-6 )
  {
    if( prpA->size() != prpB->size() ) return false;
    for( int i = 0; i < (int) prpA->size(); ++i ) if( !approxEq((*prpA)[i],(*prpB)[i],eps) ) return false;
    return true;    
  }

  ///////////////////////
  // MAT2D STUFF
  void printProperty( Property<Mat2d>* prp )
  {
    for( int i  = 0; i < (int) prp->size(); ++i )
    {
      std::cout << (*prp)[i] << " ";
    }
  }
  
  bool propertiesEqual( Property<Mat2d>* prpA, Property<Mat2d>* prpB, double eps = 1.0e-6 )
  {
    if( prpA->size() != prpB->size() ) return false;
    for( int i = 0; i < (int) prpA->size(); ++i ) if( !approxEq((*prpA)[i],(*prpB)[i],eps) ) return false;
    return true;    
  }
  

  ///////////////////////
  // MATXD STUFF
  void printProperty( Property<MatXd>* prp )
  {
    for( int i  = 0; i < (int) prp->size(); ++i )
    {
      std::cout << "(" << (*prp)[i].rows() << "," << (*prp)[i].cols() << ")-" << (*prp)[i] << " ";
    }
  }

  bool propertiesEqual( Property<MatXd>* prpA, Property<MatXd>* prpB, double eps = 1.0e-6 )
  {
    if( prpA->size() != prpB->size() ) return false;
    for( int i = 0; i < (int) prpA->size(); ++i ) if( !approxEq((*prpA)[i],(*prpB)[i],eps) ) return false;
    return true;    
  }

  ////////////////////////
  // std::vector< std::pair<Scalar,Scalar> >
  void printProperty( Property<std::vector< std::pair<Scalar,Scalar> > >* prp )
  {
    for( int i  = 0; i < (int) prp->size(); ++i )
    {
      std::cout << (*prp)[i].size() << "-[";
      for( int j = 0; j < (int) (*prp)[i].size(); ++j )
      {
        std::cout << "<" << (*prp)[i][j].first << "," << (*prp)[i][j].second << ">";
        if( j != (int) ((*prp)[i].size()-1) ) std::cout << " ";
      }
      std::cout << "] ";
    }
  }

  bool propertiesEqual( Property<std::vector< std::pair<Scalar,Scalar> > >* prpA, Property<std::vector< std::pair<Scalar,Scalar> > >* prpB, double eps = 1.0e-6 )
  {
    if( prpA->size() != prpB->size() ) return false;
    for( int i  = 0; i < (int) prpA->size(); ++i ) if( (*prpA)[i].size() != (*prpB)[i].size() ) return false;
    
    for( int i  = 0; i < (int) prpA->size(); ++i )
    {
      for( int j = 0; j < (int) (*prpA)[i].size(); ++j )
      {
        if( !approxEq((*prpA)[i][j].first,  (*prpB)[i][j].first, eps) ) return false;
        if( !approxEq((*prpA)[i][j].second, (*prpB)[i][j].second, eps) ) return false;
      }
    }
    return true;
  }
  
  //////////////////////////////////////
  // std::vector<Vec3d>
  void printProperty( Property<std::vector<Vec3d> >* prp )
  {
    for( int i  = 0; i < (int) prp->size(); ++i )
    {
      std::cout << (*prp)[i].size() << "-[";
      for( int j = 0; j < (int) (*prp)[i].size(); ++j )
      {
        std::cout << (*prp)[i][j].transpose();
        if( j != (int) ((*prp)[i].size()-1) ) std::cout << " ";
      }
      std::cout << "] ";
    }
  }
  
  bool propertiesEqual( Property<std::vector<Vec3d> >* prpA, Property<std::vector<Vec3d> >* prpB, double eps = 1.0e-6 )
  {
    if( prpA->size() != prpB->size() ) return false;
    for( int i  = 0; i < (int) prpA->size(); ++i ) if( (*prpA)[i].size() != (*prpB)[i].size() ) return false;
    
    for( int i  = 0; i < (int) prpA->size(); ++i )
    {
      for( int j = 0; j < (int) (*prpA)[i].size(); ++j )
      {
        if( !approxEq((*prpA)[i][j],(*prpB)[i][j], eps) ) return false;
      }
    }
    return true;
  }
  
  /////////////////////////////////////
  // std::pair<MatXd, MatXd>
  void printProperty( Property<std::pair<MatXd,MatXd> >* prp )
  {
    for( int i  = 0; i < (int) prp->size(); ++i )
    {
      std::cout << "<" << "(" << (*prp)[i].first.rows() << "," << (*prp)[i].first.cols() << ")-" << (*prp)[i].first << "," <<  "(" << (*prp)[i].first.rows() << "," << (*prp)[i].first.cols() << ")-" << (*prp)[i].second << "> ";
    }
  }  

  bool propertiesEqual( Property<std::pair<MatXd,MatXd> >* prpA, Property<std::pair<MatXd,MatXd> >* prpB, double eps = 1.0e-6 )
  {
    if( prpA->size() != prpB->size() ) return false;
    
    for( int i  = 0; i < (int) prpA->size(); ++i )
    {
      if( !approxEq((*prpA)[i].first,  (*prpB)[i].first,  eps) ) return false;
      if( !approxEq((*prpA)[i].second, (*prpB)[i].second, eps) ) return false;
    }
    return true;
  }

  //////////////////////////////////////
  // VertexTopology<TopologicalObject>
  void printProperty( Property<VertexTopology<TopologicalObject> >* prp )
  {
    for( int i  = 0; i < (int) prp->size(); ++i )
    {
      std::cout << "[";
      for( int j = 0; j < (int) (*prp)[i].size(); ++j )
      {
        VertexTopology<TopologicalObject>::edge_handle ehndl;
        ehndl = (*prp)[i][j];
        std::cout << ehndl.idx();
        if( j != (int)((*prp)[i].size()-1) ) std::cout << " ";
      }
      std::cout << "] ";
    }
  }

  bool propertiesEqual( Property<VertexTopology<TopologicalObject> >* prpA, Property<VertexTopology<TopologicalObject> >* prpB )
  {
    if( prpA->size() != prpB->size() ) return false;
    for( int i  = 0; i < (int) prpA->size(); ++i ) if( (*prpA)[i].size() != (*prpB)[i].size() ) return false;    
    for( int i  = 0; i < (int) prpA->size(); ++i )
    {
      for( int j = 0; j < (int) (*prpA)[i].size(); ++j )
      {
        VertexTopology<TopologicalObject>::edge_handle ehndlA;
        ehndlA = (*prpA)[i][j];
        VertexTopology<TopologicalObject>::edge_handle ehndlB;
        ehndlB = (*prpB)[i][j];
        if( ehndlA.idx() != ehndlB.idx() ) return false;
      }      
    }
    return true;
  }
  
  ////////////////////////////////////////////////////////////
  // Util::pair<Scalar,Scalar>
  void printProperty( Property<Util::pair<Scalar,Scalar> >* prp )
  {
    for( int i  = 0; i < (int) prp->size(); ++i )
    {
      std::cout << "<" << (*prp)[i].first << "," << (*prp)[i].second << "> ";
    }
  }

  bool propertiesEqual( Property<Util::pair<Scalar,Scalar> >* prpA, Property<Util::pair<Scalar,Scalar> >* prpB, double eps = 1.0e-6 )
  {
    if( prpA->size() != prpB->size() ) return false;
    for( int i = 0; i < (int) prpA->size(); ++i )
    {
      if( !approxEq((*prpA)[i].first,  (*prpB)[i].first,  eps) ) return false;
      if( !approxEq((*prpA)[i].second, (*prpB)[i].second, eps) ) return false;
    }
    return true;
  }

  //////////////////////////////////////////////////////////////
  // Util::pair<Vec3d,Vec3d>
  void printProperty( Property<Util::pair<Vec3d,Vec3d> >* prp )
  {
    for( int i  = 0; i < (int) prp->size(); ++i )
    {
      std::cout << "<" << (*prp)[i].first.transpose() << "," << (*prp)[i].second.transpose() << "> ";
    }
  }

  bool propertiesEqual( Property<Util::pair<Vec3d,Vec3d> >* prpA, Property<Util::pair<Vec3d,Vec3d> >* prpB, double eps = 1.0e-6 )
  {
    if( prpA->size() != prpB->size() ) return false;
    for( int i = 0; i < (int) prpA->size(); ++i )
    {
      if( !approxEq((*prpA)[i].first,  (*prpB)[i].first,  eps) ) return false;
      if( !approxEq((*prpA)[i].second, (*prpB)[i].second, eps) ) return false;
    }
    return true;
  }
  
  //////////////////////////////////////////////////////////////
  // Util::pair<Vec3d,Vec3d>
  void printProperty( Property<EdgeTopology<TopologicalObject> >* prp )
  {
    for( int i  = 0; i < (int) prp->size(); ++i )
    {
      std::cout << "[";
      for( int j = 0; j < (int) (*prp)[i].size(); ++j )
      {
        EdgeTopology<TopologicalObject>::vertex_handle vhndl;
        vhndl = (*prp)[i][j];
        std::cout << vhndl.idx();
        if( j != (int)((*prp)[i].size()-1) ) std::cout << " ";
      }
      std::cout << "] ";
    }
  }
  
  bool propertiesEqual( Property<EdgeTopology<TopologicalObject> >* prpA, Property<EdgeTopology<TopologicalObject> >* prpB )
  {
    if( prpA->size() != prpB->size() ) return false;
    for( int i = 0; i < (int) prpA->size(); ++i )
    {
      for( int j = 0; j < (int) (*prpA)[i].size(); ++j )
      {
        EdgeTopology<TopologicalObject>::vertex_handle vhndlA;
        vhndlA = (*prpA)[i][j];
        EdgeTopology<TopologicalObject>::vertex_handle vhndlB;
        vhndlB = (*prpB)[i][j];
        if( vhndlA.idx() != vhndlB.idx() ) return false;
      }
    }
    return true;
  }
  
  ///////////////////////////////////////////////////////
  // Face Topology Stuff
  void printProperty( Property<FaceTopology<TopologicalObject> >* prp )
  {
    for( int i  = 0; i < (int) prp->size(); ++i )
    {
      std::cout << "[";
      for( int j = 0; j < (int) (*prp)[i].size(); ++j )
      {
        FaceTopology<TopologicalObject>::vertex_handle vhndl;
        vhndl = (*prp)[i][j];
        std::cout << vhndl.idx();
        if( j != (int)((*prp)[i].size()-1) ) std::cout << " ";
      }
      std::cout << "] ";
    }
  }
  
  bool propertiesEqual( Property<FaceTopology<TopologicalObject> >* prpA, Property<FaceTopology<TopologicalObject> >* prpB )
  {
    if( prpA->size() != prpB->size() ) return false;
    for( int i = 0; i < (int) prpA->size(); ++i ) if( (*prpA)[i].size() != (*prpB)[i].size() ) return false;
    for( int i = 0; i < (int) prpA->size(); ++i )
    {
      for( int j = 0; j < (int) (*prpA)[i].size(); ++j )
      {
        FaceTopology<TopologicalObject>::vertex_handle vhndlA;
        vhndlA = (*prpA)[i][j];
        FaceTopology<TopologicalObject>::vertex_handle vhndlB;
        vhndlB = (*prpB)[i][j];
        if( vhndlA.idx() != vhndlB.idx() ) return false;
      }
    }
    return true;
  }


  PropertyContainer::Properties m_v_props;
  PropertyContainer::Properties m_e_props;
  PropertyContainer::Properties m_f_props;
  PropertyContainer::Properties m_o_props;
};

} // namespace BASim

#endif // RODSTATEBACKUP_HH

