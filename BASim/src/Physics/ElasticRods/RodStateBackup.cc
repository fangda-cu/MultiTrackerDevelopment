#include "RodStateBackup.hh"

namespace BASim 
{

RodState::RodState()
{
}

RodState::~RodState()
{
  clear();
}

void RodState::copyState( ElasticRod& rod )
{
  assert( m_v_props.size() == 0 );
  assert( m_e_props.size() == 0 );
  assert( m_f_props.size() == 0 );
  assert( m_o_props.size() == 0 );
  
  const PropertyContainer& vprops = getVertexPropertyContainer(rod);
  for( int i = 0; i < (int) vprops.getProperties().size(); ++i ) m_v_props.push_back(vprops.getProperties()[i]->clone());
  
  const PropertyContainer& eprops = getEdgePropertyContainer(rod);
  for( int i = 0; i < (int) eprops.getProperties().size(); ++i ) m_e_props.push_back(eprops.getProperties()[i]->clone());
  
  const PropertyContainer& fprops = getFacePropertyContainer(rod);
  for( int i = 0; i < (int) fprops.getProperties().size(); ++i ) m_f_props.push_back(fprops.getProperties()[i]->clone());
  
  const PropertyContainer& oprops = getObjectPropertyContainer(rod);
  for( int i = 0; i < (int) oprops.getProperties().size(); ++i ) m_o_props.push_back(oprops.getProperties()[i]->clone());
  
  #ifndef NDEBUG
    ensureNoRepeatProperties();
  #endif
}

void RodState::restoreState( ElasticRod& rod )
{
  assert( getVertexPropertyContainer(rod).getProperties().size() == m_v_props.size() );
  assert( getEdgePropertyContainer(rod).getProperties().size() == m_e_props.size() );
  assert( getFacePropertyContainer(rod).getProperties().size() == m_f_props.size() );
  assert( getObjectPropertyContainer(rod).getProperties().size() == m_o_props.size() );
  
  //std::cout << "Restoring rod: " << m_v_props.size() << " " << m_e_props.size() << " " << m_f_props.size() << std::endl;
  
  PropertyContainer& vprops = getVertexPropertyContainer(rod);
  for( int i = 0; i < (int) vprops.getProperties().size(); ++i ) vprops.getProperties()[i]->copy(m_v_props[i]);
  
  PropertyContainer& eprops = getEdgePropertyContainer(rod);
  for( int i = 0; i < (int) eprops.getProperties().size(); ++i ) eprops.getProperties()[i]->copy(m_e_props[i]);
  
  PropertyContainer& fprops = getFacePropertyContainer(rod);
  for( int i = 0; i < (int) fprops.getProperties().size(); ++i ) fprops.getProperties()[i]->copy(m_f_props[i]);
  
  PropertyContainer& oprops = getObjectPropertyContainer(rod);
  for( int i = 0; i < (int) oprops.getProperties().size(); ++i ) oprops.getProperties()[i]->copy(m_o_props[i]);
}

void RodState::setState( ElasticRod& rod )
{
  PropertyContainer::Properties& vprops = getVertexPropertyContainer(rod).getProperties();
  for( int i = 0; i < (int) vprops.size(); ++i ) { assert( vprops[i] != NULL ); delete vprops[i]; }
  vprops.resize(m_v_props.size());
  for( int i = 0; i < (int) vprops.size(); ++i ) vprops[i] = m_v_props[i]->clone();
  
  PropertyContainer::Properties& eprops = getEdgePropertyContainer(rod).getProperties();
  for( int i = 0; i < (int) eprops.size(); ++i ) { assert( eprops[i] != NULL ); delete eprops[i]; }
  eprops.resize(m_e_props.size());
  for( int i = 0; i < (int) eprops.size(); ++i ) eprops[i] = m_e_props[i]->clone();
  
  PropertyContainer::Properties& fprops = getFacePropertyContainer(rod).getProperties();
  for( int i = 0; i < (int) fprops.size(); ++i ) { assert( fprops[i] != NULL ); delete fprops[i]; }
  fprops.resize(m_f_props.size());
  for( int i = 0; i < (int) fprops.size(); ++i ) fprops[i] = m_f_props[i]->clone();
  
  PropertyContainer::Properties& oprops = getObjectPropertyContainer(rod).getProperties();
  for( int i = 0; i < (int) oprops.size(); ++i ) { assert( oprops[i] != NULL ); delete oprops[i]; }
  oprops.resize(m_o_props.size());
  for( int i = 0; i < (int) oprops.size(); ++i ) oprops[i] = m_o_props[i]->clone();
}

void RodState::ensureNoRepeatProperties()
{
  // Check vertex properties
  for( int i = 0; i < (int) m_v_props.size(); ++i )
  {
    std::string name = m_v_props[i]->name();
    for( int j = i+1; j < (int) m_v_props.size(); ++j )
    {
      if( name == m_v_props[j]->name() ) 
      {
        std::cerr << "\033[31;1mERROR IN RODSTATEBACKUP:\033[m Repeat vertex property encountered: " << name << ". Exiting." << std::endl;
        exit(1);
      }
    }
  }
  
  // Check edge properties
  for( int i = 0; i < (int) m_e_props.size(); ++i )
  {
    std::string name = m_e_props[i]->name();
    for( int j = i+1; j < (int) m_e_props.size(); ++j )
    {
      if( name == m_e_props[j]->name() ) 
      {
        std::cerr << "\033[31;1mERROR IN RODSTATEBACKUP:\033[m Repeat edge property encountered: " << name << ". Exiting." << std::endl;
        exit(1);
      }
    }
  }
  
  // Check face properties
  for( int i = 0; i < (int) m_f_props.size(); ++i )
  {
    std::string name = m_f_props[i]->name();
    for( int j = i+1; j < (int) m_f_props.size(); ++j )
    {
      if( name == m_f_props[j]->name() ) 
      {
        std::cerr << "\033[31;1mERROR IN RODSTATEBACKUP:\033[m Repeat face property encountered: " << name << ". Exiting." << std::endl;
        exit(1);
      }
    }
  }
  
  // Check object properties
  for( int i = 0; i < (int) m_o_props.size(); ++i )
  {
    std::string name = m_o_props[i]->name();
    for( int j = i+1; j < (int) m_o_props.size(); ++j )
    {
      if( name == m_o_props[j]->name() ) 
      {
        std::cerr << "\033[31;1mERROR IN RODSTATEBACKUP:\033[m Repeat object property encountered: " << name << ". Exiting." << std::endl;
        exit(1);
      }
    }
  }
}

void RodState::clear()
{
  for( int i = 0; i < (int) m_v_props.size(); ++i ) { assert( m_v_props[i] != NULL ); delete m_v_props[i]; }
  m_v_props.clear();
  for( int i = 0; i < (int) m_e_props.size(); ++i ) { assert( m_e_props[i] != NULL ); delete m_e_props[i]; }
  m_e_props.clear();
  for( int i = 0; i < (int) m_f_props.size(); ++i ) { assert( m_f_props[i] != NULL ); delete m_f_props[i]; }
  m_f_props.clear();
  for( int i = 0; i < (int) m_o_props.size(); ++i ) { assert( m_o_props[i] != NULL ); delete m_o_props[i]; }
  m_o_props.clear();
}

void RodState::print( ElasticRod& rod )
{
  PropertyContainer::Properties v_props = getVertexPropertyContainer(rod).getProperties();
  PropertyContainer::Properties e_props = getEdgePropertyContainer(rod).getProperties();
  PropertyContainer::Properties f_props = getFacePropertyContainer(rod).getProperties();
  PropertyContainer::Properties o_props = getObjectPropertyContainer(rod).getProperties();
  
  //std::cout << "Backed rods up: " << m_v_props.size() << " " << m_e_props.size() << " " << m_f_props.size() << std::endl;
  std::cout << "Vertex properties (" << v_props.size() << ")" << std::endl;
  for( int i = 0; i < (int) v_props.size(); ++i )
  {
    assert( v_props[i] != NULL );
    
    std::cout << "\033[34;1m" << v_props[i]->name() << "\033[m"; std::cout.flush();
    if( dynamic_cast<Property<Scalar>*>(v_props[i]) )
    {
      std::cout << "(" << dynamic_cast<Property<Scalar>*>(v_props[i]) << ")";
      std::cout << " : " << "Scalar";
      std::cout << "   "; BASim::printProperty(dynamic_cast<Property<Scalar>*>(v_props[i]));
    }
    else if( dynamic_cast<Property<int>*>(v_props[i]) )
    {
      std::cout << "(" << dynamic_cast<Property<int>*>(v_props[i]) << ")";
      std::cout << " : " << "int";
      std::cout << "   "; printProperty(dynamic_cast<Property<int>*>(v_props[i]));
    }
    else if( dynamic_cast<Property<bool>*>(v_props[i]) )
    {
      std::cout << "(" << dynamic_cast<Property<bool>*>(v_props[i]) << ")";
      std::cout << " : " << "bool";
      std::cout << "   "; printProperty(dynamic_cast<Property<bool>*>(v_props[i]));
    }
    else if( dynamic_cast<Property<Vec2d>*>(v_props[i]) )
    {
      std::cout << "(" << dynamic_cast<Property<Vec2d>*>(v_props[i]) << ")";
      std::cout << " : " << "Vec2d";
      std::cout << "   "; printProperty(dynamic_cast<Property<Vec2d>*>(v_props[i]));
    }
    else if( dynamic_cast<Property<Vec3d>*>(v_props[i]) )
    {
      std::cout << "(" << dynamic_cast<Property<Vec3d>*>(v_props[i]) << ")";
      std::cout << " : " << "Vec3d";
      std::cout << "   "; BASim::printProperty(dynamic_cast<Property<Vec3d>*>(v_props[i]));
    }
    else if( dynamic_cast<Property<VecXd>*>(v_props[i]) )
    {
      std::cout << "(" << dynamic_cast<Property<VecXd>*>(v_props[i]) << ")";
      std::cout << " : " << "VecXd";
      std::cout << "   "; printProperty(dynamic_cast<Property<VecXd>*>(v_props[i]));
    }
    else if( dynamic_cast<Property<Mat2d>*>(v_props[i]) )
    {
      std::cout << "(" << dynamic_cast<Property<Mat2d>*>(v_props[i]) << ")";
      std::cout << " : " << "Mat2d";
      std::cout << "   "; printProperty(dynamic_cast<Property<Mat2d>*>(v_props[i]));
    }
    else if( dynamic_cast<Property<MatXd>*>(v_props[i]) )
    {
      std::cout << "(" << dynamic_cast<Property<MatXd>*>(v_props[i]) << ")";
      std::cout << " : " << "MatXd";
      std::cout << "   "; printProperty(dynamic_cast<Property<MatXd>*>(v_props[i]));
    }
    else if( dynamic_cast<Property<std::vector< std::pair<Scalar,Scalar> > >*>(v_props[i]) )
    {
      std::cout << "(" << dynamic_cast<Property<std::vector< std::pair<Scalar,Scalar> > >*>(v_props[i]) << ")";
      std::cout << " : " << "std::vector< std::pair<Scalar,Scalar> >";
      std::cout << "   "; printProperty(dynamic_cast<Property<std::vector< std::pair<Scalar,Scalar> > >*>(v_props[i]));
    }
    else if( dynamic_cast<Property< std::vector<Vec3d> >*>(v_props[i]) )
    {
      std::cout << "(" << dynamic_cast<Property<std::vector<Vec3d> >*>(v_props[i]) << ")";
      std::cout << " : " << "std::vector<Vec3d>";
      std::cout << "   "; printProperty(dynamic_cast<Property<std::vector<Vec3d> >*>(v_props[i]));
    }
    else if( dynamic_cast<Property< std::pair<MatXd, MatXd> >*>(v_props[i]) )
    {
      std::cout << "(" << dynamic_cast<Property<std::pair<MatXd,MatXd> >*>(v_props[i]) << ")";
      std::cout << " : " << "std::pair<MatXd,MatXd>";
      std::cout << "   "; printProperty(dynamic_cast<Property<std::pair<MatXd,MatXd> >*>(v_props[i]));
    }
    else if( dynamic_cast<Property<VertexTopology<TopologicalObject> >*>(v_props[i]) )
    {
      std::cout << "(" << dynamic_cast<Property<VertexTopology<TopologicalObject> >*>(v_props[i]) << ")";
      std::cout << " : " << "VertexTopology<TopologicalObject>";
      std::cout << "   "; printProperty(dynamic_cast<Property<VertexTopology<TopologicalObject> >*>(v_props[i]));
    }
    else 
    {
      std::cerr << "\033[31;1mERROR IN RODSTATEBACKUP:\033[m Unrecognized veretx property type encountered for property: " << v_props[i]->name() << ". Exiting." << std::endl;
      exit(1);
    }
    
    std::cout << std::endl;
  }
  
  std::cout << "Edge properties (" << e_props.size() << ")" << std::endl;
  for( int i = 0; i < (int) e_props.size(); ++i )
  {
    std::cout << "\033[35;1m" << e_props[i]->name()  << "\033[m";
    
    if( dynamic_cast<Property<Scalar>*>(e_props[i]) )
    {
      std::cout << "(" << dynamic_cast<Property<Scalar>*>(e_props[i]) << ")";
      std::cout << " : " << "Scalar";
      std::cout << "   "; BASim::printProperty(dynamic_cast<Property<Scalar>*>(e_props[i]));
    }
    else if( dynamic_cast<Property<int>*>(e_props[i]) )
    {
      std::cout << "(" << dynamic_cast<Property<int>*>(e_props[i]) << ")";
      std::cout << " : " << "int";
      std::cout << "   "; printProperty(dynamic_cast<Property<int>*>(e_props[i]));
    }
    else if( dynamic_cast<Property<bool>*>(e_props[i]) )
    {
      std::cout << "(" << dynamic_cast<Property<bool>*>(e_props[i]) << ")";
      std::cout << " : " << "bool";
      std::cout << "   "; printProperty(dynamic_cast<Property<bool>*>(e_props[i]));
    }
    else if( dynamic_cast<Property<Vec3d>*>(e_props[i]) )
    {
      std::cout << "(" << dynamic_cast<Property<Vec3d>*>(e_props[i]) << ")";
      std::cout << " : " << "Vec3d";
      std::cout << "   "; BASim::printProperty(dynamic_cast<Property<Vec3d>*>(e_props[i]));
    }
    else if( dynamic_cast<Property<Util::pair<Scalar, Scalar> >*>(e_props[i]) )
    {
      std::cout << "(" << dynamic_cast<Property<Util::pair<Scalar, Scalar> >*>(e_props[i]) << ")";
      std::cout << " : " << "Util::pair<Scalar, Scalar>";
      std::cout << "   "; printProperty(dynamic_cast<Property<Util::pair<Scalar, Scalar> >*>(e_props[i]));
    }
    else if( dynamic_cast<Property< Util::pair<Vec3d,Vec3d> >* >(e_props[i]) )
    {
      std::cout << "(" << dynamic_cast<Property< Util::pair<Vec3d,Vec3d> >* >(e_props[i]) << ")";
      std::cout << " : " << "Util::pair<Vec3d,Vec3d>";
      std::cout << "   "; printProperty(dynamic_cast<Property< Util::pair<Vec3d,Vec3d> >* >(e_props[i]));
    }
    else if( dynamic_cast<Property<EdgeTopology<TopologicalObject> >*>(e_props[i]) )
    {
      std::cout << "(" << dynamic_cast<Property<EdgeTopology<TopologicalObject> >*>(e_props[i]) << ")";
      std::cout << " : " << "EdgeTopology<TopologicalObject>";
      std::cout << "   "; printProperty(dynamic_cast<Property<EdgeTopology<TopologicalObject> >*>(e_props[i]));
    }
    else 
    {
      std::cerr << "\033[31;1mERROR IN RODSTATEBACKUP:\033[m Unrecognized edge property type encountered for property: " << e_props[i]->name() << ". Exiting." << std::endl;
      exit(1);
    }
    std::cout << std::endl;
  }
  
  std::cout << "Face properties (" << f_props.size() << ")" << std::endl;
  for( int i = 0; i < (int) f_props.size(); ++i )
  {
    std::cout << "\033[36;1m" << f_props[i]->name()  << "\033[m";
    
    if( dynamic_cast<Property<FaceTopology<TopologicalObject> >*>(f_props[i]) )
    {
      std::cout << "(" << dynamic_cast<Property<FaceTopology<TopologicalObject> >*>(f_props[i]) << ")";
      std::cout << " : " << "FaceTopology<TopologicalObject>";
      std::cout << "   "; printProperty(dynamic_cast<Property<FaceTopology<TopologicalObject> >*>(f_props[i]));
    }
    else 
    {
      std::cerr << "\033[31;1mERROR IN RODSTATEBACKUP:\033[m Unrecognized face property type encountered for property: " << f_props[i]->name() << ". Exiting." << std::endl;
      exit(1);
    }
    
    std::cout << std::endl;
  }
  
  std::cout << "Object properties (" << o_props.size() << ")" << std::endl;
  for( int i = 0; i < (int) o_props.size(); ++i )
  {
    std::cout << "\033[33;1m" << o_props[i]->name()  << "\033[m";
    
    if( dynamic_cast<Property<int>*>(o_props[i]) )
    {
      std::cout << "(" << dynamic_cast<Property<int>*>(o_props[i]) << ")";
      std::cout << " : " << "int";
      std::cout << "   "; printProperty(dynamic_cast<Property<int>*>(o_props[i]));
    }
    else if( dynamic_cast<Property<Scalar>*>(o_props[i]) )
    {
      std::cout << "(" << dynamic_cast<Property<Scalar>*>(o_props[i]) << ")";
      std::cout << " : " << "Scalar";
      std::cout << "   "; BASim::printProperty(dynamic_cast<Property<Scalar>*>(o_props[i]));
    }
    else if( dynamic_cast<Property<bool>*>(o_props[i]) )
    {
      std::cout << "(" << dynamic_cast<Property<bool>*>(o_props[i]) << ")";
      std::cout << " : " << "bool";
      std::cout << "   "; printProperty(dynamic_cast<Property<bool>*>(o_props[i]));
    }
    else if( dynamic_cast<Property<DOFMap>*>(o_props[i]) )
    {
      std::cout << "(" << dynamic_cast<Property<DOFMap>*>(o_props[i]) << ")";
      std::cout << " : " << "DOFMap";
      std::cout << "   "; printProperty(dynamic_cast<Property<DOFMap>*>(o_props[i]));
    }
    else if( dynamic_cast<Property<ElasticRod::RodForces>*>(o_props[i]) )
    {
      std::cout << "(" << dynamic_cast<Property<ElasticRod::RodForces>*>(o_props[i]) << ")";
      std::cout << " : " << "ElasticRod::RodForces";
      std::cout << "   "; printProperty(dynamic_cast<Property<ElasticRod::RodForces>*>(o_props[i]));
    }
    else if( dynamic_cast<Property<ElasticRod::RefFrameType>*>(o_props[i]) )
    {
      std::cout << "(" << dynamic_cast<Property<ElasticRod::RefFrameType>*>(o_props[i]) << ")";
      std::cout << " : " << "ElasticRod::RefFrameType";
      std::cout << "   "; printProperty(dynamic_cast<Property<ElasticRod::RefFrameType>*>(o_props[i]));
    }
    else if( dynamic_cast<Property<RodBoundaryCondition::BCList>*>(o_props[i]) )
    {
      std::cout << "(" << dynamic_cast<Property<RodBoundaryCondition::BCList>*>(o_props[i]) << ")";
      std::cout << " : " << "RodBoundaryCondition::BCList";
      std::cout << "   "; printProperty(dynamic_cast<Property<RodBoundaryCondition::BCList>*>(o_props[i]));
    }
    else 
    {
      std::cerr << "\033[31;1mERROR IN RODSTATEBACKUP:\033[m Unrecognized object property type encountered for property: " << o_props[i]->name() << ". Exiting." << std::endl;
      //exit(1);
    }
    
    std::cout << std::endl;
  }
}

bool RodState::compareProperties( const ElasticRod& other, bool printdiffs )
{
  PropertyContainer::Properties o_o_props = getObjectPropertyContainer(other).getProperties();
  PropertyContainer::Properties o_v_props = getVertexPropertyContainer(other).getProperties();
  PropertyContainer::Properties o_e_props = getEdgePropertyContainer(other).getProperties();
  PropertyContainer::Properties o_f_props = getFacePropertyContainer(other).getProperties();
  
  // Ensure object property containers are the same size
  if( m_o_props.size() != o_o_props.size() ) 
  {
    if( printdiffs ) std::cout << "\033[31;1mWARNING OBJECT PROPERTY CONTAINERS OF DIFFERENT SIZE.\033[m" << std::endl;
    return false;
  }
  // Ensure object property containers are the same size
  if( m_v_props.size() != o_v_props.size() ) 
  {
    if( printdiffs ) std::cout << "\033[31;1mWARNING VERTEX PROPERTY CONTAINERS OF DIFFERENT SIZE.\033[m" << std::endl;
    return false;
  }
  // Ensure object property containers are the same size
  if( m_e_props.size() != o_e_props.size() ) 
  {
    if( printdiffs ) std::cout << "\033[31;1mWARNING EDGE PROPERTY CONTAINERS OF DIFFERENT SIZE.\033[m" << std::endl;
    return false;
  }
  // Ensure object property containers are the same size
  if( m_f_props.size() != o_f_props.size() ) 
  {
    if( printdiffs ) std::cout << "\033[31;1mWARNING FACE PROPERTY CONTAINERS OF DIFFERENT SIZE.\033[m" << std::endl;
    return false;
  }
  
  if( printdiffs ) std::cout << "\033[33;1mCOMPARING OBJECT PROPERTIES\033[m" << std::endl;
  for( int i = 0; i < (int) o_o_props.size(); ++i )
  {
    assert( o_o_props[i] != NULL );
    assert( m_o_props[i] != NULL );
    
    if( o_o_props[i]->name() != m_o_props[i]->name() ) 
    {
      if( printdiffs ) std::cout << "\033[31;1mWARNING OBJECT PROPERTYS IN DIFFERENT ORDER.\033[m" << std::endl;
      return false;
    }
    
    if( printdiffs ) std::cout << "\033[33;1m  COMPARING: " << o_o_props[i]->name() << "\033[m: ";
    
    if( dynamic_cast<Property<int>*>(o_o_props[i]) )
    {
      if( !dynamic_cast<Property<int>*>(m_o_props[i]) )
      {
        if( printdiffs ) std::cout << "\033[31;1mWARNING OBJECT PROPERTYS OF SAME NAME NOT THE SAME TYPE.\033[m" << std::endl;
        return false;
      }
      if( !propertiesEqual(dynamic_cast<Property<int>*>(o_o_props[i]),dynamic_cast<Property<int>*>(m_o_props[i])) )
      {
        if( printdiffs ) std::cout << "\033[31;1mWARNING OBJECT PROPERTYS NOT EQUAL.\033[m" << std::endl;
        return false;
      }
      if( printdiffs ) std::cout << "EQUAL";
    }
    else if( dynamic_cast<Property<bool>*>(o_o_props[i]) )
    {
      if( !dynamic_cast<Property<bool>*>(m_o_props[i]) )
      {
        if( printdiffs ) std::cout << "\033[31;1mWARNING OBJECT PROPERTYS OF SAME NAME NOT THE SAME TYPE.\033[m" << std::endl;
        return false;
      }
      if( !propertiesEqual(dynamic_cast<Property<bool>*>(o_o_props[i]),dynamic_cast<Property<bool>*>(m_o_props[i])) )
      {
        if( printdiffs ) std::cout << "\033[31;1mWARNING OBJECT PROPERTYS NOT EQUAL.\033[m" << std::endl;
        return false;
      }
      if( printdiffs ) std::cout << "EQUAL";
    }
    else if( dynamic_cast<Property<Scalar>*>(o_o_props[i]) )
    {
      if( !dynamic_cast<Property<Scalar>*>(m_o_props[i]) )
      {
        if( printdiffs ) std::cout << "\033[31;1mWARNING OBJECT PROPERTYS OF SAME NAME NOT THE SAME TYPE.\033[m" << std::endl;
        return false;
      }
      if( !BASim::propertiesEqual(dynamic_cast<Property<Scalar>*>(o_o_props[i]),dynamic_cast<Property<Scalar>*>(m_o_props[i]),1.0e-9) )
      {
        if( printdiffs ) std::cout << "\033[31;1mWARNING OBJECT PROPERTYS NOT EQUAL.\033[m" << std::endl;
        return false;
      }
      if( printdiffs ) std::cout << "EQUAL";
    }
    else if( dynamic_cast<Property<ElasticRod::RefFrameType>*>(o_o_props[i]) )
    {
      if( !dynamic_cast<Property<ElasticRod::RefFrameType>*>(m_o_props[i]) )
      {
        if( printdiffs ) std::cout << "\033[31;1mWARNING OBJECT PROPERTYS OF SAME NAME NOT THE SAME TYPE.\033[m" << std::endl;
        return false;
      }
      if( !propertiesEqual(dynamic_cast<Property<ElasticRod::RefFrameType>*>(o_o_props[i]),dynamic_cast<Property<ElasticRod::RefFrameType>*>(m_o_props[i])) )
      {
        if( printdiffs ) std::cout << "\033[31;1mWARNING OBJECT PROPERTYS NOT EQUAL.\033[m" << std::endl;
        return false;
      }
      if( printdiffs ) std::cout << "EQUAL";
    }      
    else if( dynamic_cast<Property<RodBoundaryCondition::BCList>*>(o_o_props[i]) )
    {
      if( !dynamic_cast<Property<RodBoundaryCondition::BCList>*>(m_o_props[i]) )
      {
        if( printdiffs ) std::cout << "\033[31;1mWARNING OBJECT PROPERTYS OF SAME NAME NOT THE SAME TYPE.\033[m" << std::endl;
        return false;
      }
      if( !propertiesEqual(dynamic_cast<Property<RodBoundaryCondition::BCList>*>(o_o_props[i]),dynamic_cast<Property<RodBoundaryCondition::BCList>*>(m_o_props[i])) )
      {
        if( printdiffs ) std::cout << "\033[31;1mWARNING OBJECT PROPERTYS NOT EQUAL.\033[m" << std::endl;
        return false;
      }
      if( printdiffs ) std::cout << "EQUAL";        
    }
    else if( dynamic_cast<Property<ElasticRod::RodForces>*>(o_o_props[i]) )
    {
      if( !dynamic_cast<Property<ElasticRod::RodForces>*>(m_o_props[i]) )
      {
        if( printdiffs ) std::cout << "\033[31;1mWARNING OBJECT PROPERTYS OF SAME NAME NOT THE SAME TYPE.\033[m" << std::endl;
        return false;
      }
      if( !propertiesEqual(dynamic_cast<Property<ElasticRod::RodForces>*>(o_o_props[i]),dynamic_cast<Property<ElasticRod::RodForces>*>(m_o_props[i])) )
      {
        if( printdiffs ) std::cout << "\033[31;1mWARNING OBJECT PROPERTYS NOT EQUAL.\033[m" << std::endl;
        return false;
      }
      if( printdiffs ) std::cout << "EQUAL";        
    }
    else if( dynamic_cast<Property<DOFMap>*>(o_o_props[i]) )
    {
      if( !dynamic_cast<Property<DOFMap>*>(m_o_props[i]) )
      {
        if( printdiffs ) std::cout << "\033[31;1mWARNING OBJECT PROPERTYS OF SAME NAME NOT THE SAME TYPE.\033[m" << std::endl;
        return false;
      }
      if( !propertiesEqual(dynamic_cast<Property<DOFMap>*>(o_o_props[i]),dynamic_cast<Property<DOFMap>*>(m_o_props[i])) )
      {
        if( printdiffs ) std::cout << "\033[31;1mWARNING OBJECT PROPERTYS NOT EQUAL.\033[m" << std::endl;
        return false;
      }
      if( printdiffs ) std::cout << "EQUAL";        
    }
    else 
    {
      if( printdiffs ) std::cerr << "\033[31;1mERROR IN RODSTATEBACKUP:\033[m Unrecognized object property type encountered for property: " << o_o_props[i]->name() << ". Exiting." << std::endl;
      exit(1);
    }
    if( printdiffs ) std::cout << std::endl;      
  }
  
  if( printdiffs ) std::cout << "\033[34;1mCOMPARING VERTEX PROPERTIES\033[m" << std::endl;
  for( int i = 0; i < (int) o_v_props.size(); ++i )
  {
    assert( o_v_props[i] != NULL );
    assert( m_v_props[i] != NULL );
    
    if( o_v_props[i]->name() != m_v_props[i]->name() ) 
    {
      if( printdiffs ) std::cout << "\033[31;1mWARNING VERTEX PROPERTIES IN DIFFERENT ORDER.\033[m" << std::endl;
      return false;
    }
    
    if( printdiffs ) std::cout << "\033[34;1m  COMPARING: " << o_v_props[i]->name() << "\033[m: ";
    
    if( dynamic_cast<Property<int>*>(o_v_props[i]) )
    {
      if( !dynamic_cast<Property<int>*>(m_v_props[i]) )
      {
        if( printdiffs ) std::cout << "\033[31;1mWARNING VERTEX PROPERTYS OF SAME NAME NOT THE SAME TYPE.\033[m" << std::endl;
        return false;
      }
      if( !propertiesEqual(dynamic_cast<Property<int>*>(o_v_props[i]),dynamic_cast<Property<int>*>(m_v_props[i])) )
      {
        if( printdiffs ) std::cout << "\033[31;1mWARNING VERTEX PROPERTYS NOT EQUAL.\033[m" << std::endl;
        return false;
      }
      if( printdiffs ) std::cout << "EQUAL";
    }      
    else if( dynamic_cast<Property<bool>*>(o_v_props[i]) )
    {
      if( !dynamic_cast<Property<bool>*>(m_v_props[i]) )
      {
        if( printdiffs ) std::cout << "\033[31;1mWARNING VERTEX PROPERTYS OF SAME NAME NOT THE SAME TYPE.\033[m" << std::endl;
        return false;
      }
      if( !propertiesEqual(dynamic_cast<Property<bool>*>(o_v_props[i]),dynamic_cast<Property<bool>*>(m_v_props[i])) )
      {
        if( printdiffs ) std::cout << "\033[31;1mWARNING VERTEX PROPERTYS NOT EQUAL.\033[m" << std::endl;
        return false;
      }
      if( printdiffs ) std::cout << "EQUAL";
    }      
    else if( dynamic_cast<Property<Scalar>*>(o_v_props[i]) )
    {
      if( !dynamic_cast<Property<Scalar>*>(m_v_props[i]) )
      {
        if( printdiffs ) std::cout << "\033[31;1mWARNING VERTEX PROPERTYS OF SAME NAME NOT THE SAME TYPE.\033[m" << std::endl;
        return false;
      }
      if( !BASim::propertiesEqual(dynamic_cast<Property<Scalar>*>(o_v_props[i]),dynamic_cast<Property<Scalar>*>(m_v_props[i]),1.0e-9) )
      {
        if( printdiffs ) std::cout << "\033[31;1mWARNING VERTEX PROPERTYS NOT EQUAL.\033[m" << std::endl;
        return false;
      }
      if( printdiffs ) std::cout << "EQUAL";
    }
    else if( dynamic_cast<Property<Vec2d>*>(o_v_props[i]) )
    {
      if( !dynamic_cast<Property<Vec2d>*>(m_v_props[i]) )
      {
        if( printdiffs ) std::cout << "\033[31;1mWARNING VERTEX PROPERTYS OF SAME NAME NOT THE SAME TYPE.\033[m" << std::endl;
        return false;
      }
      if( !propertiesEqual(dynamic_cast<Property<Vec2d>*>(o_v_props[i]),dynamic_cast<Property<Vec2d>*>(m_v_props[i]),1.0e-9) )
      {
        if( printdiffs ) std::cout << "\033[31;1mWARNING VERTEX PROPERTYS NOT EQUAL.\033[m" << std::endl;
        return false;
      }
      if( printdiffs ) std::cout << "EQUAL";
    }
    else if( dynamic_cast<Property<Vec3d>*>(o_v_props[i]) )
    {
      if( !dynamic_cast<Property<Vec3d>*>(m_v_props[i]) )
      {
        if( printdiffs ) std::cout << "\033[31;1mWARNING VERTEX PROPERTYS OF SAME NAME NOT THE SAME TYPE.\033[m" << std::endl;
        return false;
      }
      if( !BASim::propertiesEqual(dynamic_cast<Property<Vec3d>*>(o_v_props[i]),dynamic_cast<Property<Vec3d>*>(m_v_props[i]),1.0e-9) )
      {
        if( printdiffs ) std::cout << "\033[31;1mWARNING VERTEX PROPERTYS NOT EQUAL.\033[m" << std::endl;
        return false;
      }
      if( printdiffs ) std::cout << "EQUAL";
    }
    else if( dynamic_cast<Property<VecXd>*>(o_v_props[i]) )
    {
      if( !dynamic_cast<Property<VecXd>*>(m_v_props[i]) )
      {
        if( printdiffs ) std::cout << "\033[31;1mWARNING VERTEX PROPERTYS OF SAME NAME NOT THE SAME TYPE.\033[m" << std::endl;
        return false;
      }
      if( !propertiesEqual(dynamic_cast<Property<VecXd>*>(o_v_props[i]),dynamic_cast<Property<VecXd>*>(m_v_props[i]),1.0e-9) )
      {
        if( printdiffs ) std::cout << "\033[31;1mWARNING VERTEX PROPERTYS NOT EQUAL.\033[m" << std::endl;
        return false;
      }
      if( printdiffs ) std::cout << "EQUAL";
    }
    else if( dynamic_cast<Property<Mat2d>*>(o_v_props[i]) )
    {
      if( !dynamic_cast<Property<Mat2d>*>(m_v_props[i]) )
      {
        if( printdiffs ) std::cout << "\033[31;1mWARNING VERTEX PROPERTYS OF SAME NAME NOT THE SAME TYPE.\033[m" << std::endl;
        return false;
      }
      if( !propertiesEqual(dynamic_cast<Property<Mat2d>*>(o_v_props[i]),dynamic_cast<Property<Mat2d>*>(m_v_props[i]),1.0e-9) )
      {
        if( printdiffs ) std::cout << "\033[31;1mWARNING VERTEX PROPERTYS NOT EQUAL.\033[m" << std::endl;
        return false;
      }
      if( printdiffs ) std::cout << "EQUAL";
    }
    else if( dynamic_cast<Property<MatXd>*>(o_v_props[i]) )
    {
      if( !dynamic_cast<Property<MatXd>*>(m_v_props[i]) )
      {
        if( printdiffs ) std::cout << "\033[31;1mWARNING VERTEX PROPERTYS OF SAME NAME NOT THE SAME TYPE.\033[m" << std::endl;
        return false;
      }
      if( !propertiesEqual(dynamic_cast<Property<MatXd>*>(o_v_props[i]),dynamic_cast<Property<MatXd>*>(m_v_props[i]),1.0e-9) )
      {
        if( printdiffs ) std::cout << "\033[31;1mWARNING VERTEX PROPERTYS NOT EQUAL.\033[m" << std::endl;
        return false;
      }
      if( printdiffs ) std::cout << "EQUAL";
    }
    else if( dynamic_cast<Property<std::vector< std::pair<Scalar,Scalar> > >*>(o_v_props[i]) )
    {
      if( !dynamic_cast<Property<std::vector< std::pair<Scalar,Scalar> > >*>(m_v_props[i]) )
      {
        if( printdiffs ) std::cout << "\033[31;1mWARNING VERTEX PROPERTYS OF SAME NAME NOT THE SAME TYPE.\033[m" << std::endl;
        return false;
      }
      if( !propertiesEqual(dynamic_cast<Property<std::vector< std::pair<Scalar,Scalar> > >*>(o_v_props[i]),dynamic_cast<Property<std::vector< std::pair<Scalar,Scalar> > >*>(m_v_props[i]),1.0e-9) )
      {
        if( printdiffs ) std::cout << "\033[31;1mWARNING VERTEX PROPERTYS NOT EQUAL.\033[m" << std::endl;
        return false;
      }
      if( printdiffs ) std::cout << "EQUAL";
    }
    else if( dynamic_cast<Property< std::vector<Vec3d> >*>(o_v_props[i]) )
    {
      if( !dynamic_cast<Property< std::vector<Vec3d> >*>(m_v_props[i]) )
      {
        if( printdiffs ) std::cout << "\033[31;1mWARNING VERTEX PROPERTYS OF SAME NAME NOT THE SAME TYPE.\033[m" << std::endl;
        return false;
      }
      if( !propertiesEqual(dynamic_cast<Property< std::vector<Vec3d> >*>(o_v_props[i]),dynamic_cast<Property< std::vector<Vec3d> >*>(m_v_props[i]),1.0e-9) )
      {
        if( printdiffs ) std::cout << "\033[31;1mWARNING VERTEX PROPERTYS NOT EQUAL.\033[m" << std::endl;
        return false;
      }
      if( printdiffs ) std::cout << "EQUAL";
    }
    else if( dynamic_cast<Property< std::pair<MatXd, MatXd> >*>(o_v_props[i]) )
    {
      if( !dynamic_cast<Property< std::pair<MatXd, MatXd> >*>(m_v_props[i]) )
      {
        if( printdiffs ) std::cout << "\033[31;1mWARNING VERTEX PROPERTYS OF SAME NAME NOT THE SAME TYPE.\033[m" << std::endl;
        return false;
      }
      if( !propertiesEqual(dynamic_cast<Property< std::pair<MatXd, MatXd> >*>(o_v_props[i]),dynamic_cast<Property< std::pair<MatXd, MatXd> >*>(m_v_props[i]),1.0e-9) )
      {
        if( printdiffs ) std::cout << "\033[31;1mWARNING VERTEX PROPERTYS NOT EQUAL.\033[m" << std::endl;
        return false;
      }
      if( printdiffs ) std::cout << "EQUAL";
    }
    else if( dynamic_cast<Property<VertexTopology<TopologicalObject> >*>(o_v_props[i]) )
    {
      if( !dynamic_cast<Property<VertexTopology<TopologicalObject> >*>(m_v_props[i]) )
      {
        if( printdiffs ) std::cout << "\033[31;1mWARNING VERTEX PROPERTYS OF SAME NAME NOT THE SAME TYPE.\033[m" << std::endl;
        return false;
      }
      if( !propertiesEqual(dynamic_cast<Property<VertexTopology<TopologicalObject> >*>(o_v_props[i]),dynamic_cast<Property<VertexTopology<TopologicalObject> >*>(m_v_props[i])) )
      {
        if( printdiffs ) std::cout << "\033[31;1mWARNING VERTEX PROPERTYS NOT EQUAL.\033[m" << std::endl;
        return false;
      }
      if( printdiffs ) std::cout << "EQUAL";
    }
    else 
    {
      if( printdiffs ) std::cerr << "\033[31;1mERROR IN RODSTATEBACKUP:\033[m Unrecognized veretx property type encountered for property: " << o_v_props[i]->name() << ". Exiting." << std::endl;
      exit(1);
    }
    if( printdiffs ) std::cout << std::endl;
  }
  
  if( printdiffs ) std::cout << "\033[35;1mCOMPARING EDGE PROPERTIES\033[m" << std::endl;
  for( int i = 0; i < (int) o_e_props.size(); ++i )
  {
    assert( o_e_props[i] != NULL );
    assert( m_e_props[i] != NULL );
    
    if( o_e_props[i]->name() != m_e_props[i]->name() ) 
    {
      if( printdiffs ) std::cout << "\033[35;1mWARNING EDGE PROPERTIES IN DIFFERENT ORDER.\033[m" << std::endl;
      return false;
    }
    
    if( printdiffs ) std::cout << "\033[35;1m  COMPARING: " << o_e_props[i]->name() << "\033[m: ";
    
    if( dynamic_cast<Property<int>*>(o_e_props[i]) )
    {
      if( !dynamic_cast<Property<int>*>(m_e_props[i]) )
      {
        if( printdiffs ) std::cout << "\033[31;1mWARNING EDGE PROPERTIES OF SAME NAME NOT THE SAME TYPE.\033[m" << std::endl;
        return false;
      }
      if( !propertiesEqual(dynamic_cast<Property<int>*>(o_e_props[i]),dynamic_cast<Property<int>*>(m_e_props[i])) )
      {
        if( printdiffs ) std::cout << "\033[31;1mWARNING EDGE PROPERTIES NOT EQUAL.\033[m" << std::endl;
        return false;
      }
      if( printdiffs ) std::cout << "EQUAL";
    }
    else if( dynamic_cast<Property<bool>*>(o_e_props[i]) )
    {
      if( !dynamic_cast<Property<bool>*>(m_e_props[i]) )
      {
        if( printdiffs ) std::cout << "\033[31;1mWARNING EDGE PROPERTIES OF SAME NAME NOT THE SAME TYPE.\033[m" << std::endl;
        return false;
      }
      if( !propertiesEqual(dynamic_cast<Property<bool>*>(o_e_props[i]),dynamic_cast<Property<bool>*>(m_e_props[i])) )
      {
        if( printdiffs ) std::cout << "\033[31;1mWARNING EDGE PROPERTIES NOT EQUAL.\033[m" << std::endl;
        return false;
      }
      if( printdiffs ) std::cout << "EQUAL";
    }
    else if( dynamic_cast<Property<Scalar>*>(o_e_props[i]) )
    {
      if( !dynamic_cast<Property<Scalar>*>(m_e_props[i]) )
      {
        if( printdiffs ) std::cout << "\033[31;1mWARNING EDGE PROPERTIES OF SAME NAME NOT THE SAME TYPE.\033[m" << std::endl;
        return false;
      }
      if( !BASim::propertiesEqual(dynamic_cast<Property<Scalar>*>(o_e_props[i]),dynamic_cast<Property<Scalar>*>(m_e_props[i]),1.0e-9) )
      {
        if( printdiffs ) std::cout << "\033[31;1mWARNING EDGE PROPERTIES NOT EQUAL.\033[m" << std::endl;
        return false;
      }
      if( printdiffs ) std::cout << "EQUAL";
    }      
    else if( dynamic_cast<Property<Vec3d>*>(o_e_props[i]) )
    {
      if( !dynamic_cast<Property<Vec3d>*>(m_e_props[i]) )
      {
        if( printdiffs ) std::cout << "\033[31;1mWARNING EDGE PROPERTIES OF SAME NAME NOT THE SAME TYPE.\033[m" << std::endl;
        return false;
      }
      if( !BASim::propertiesEqual(dynamic_cast<Property<Vec3d>*>(o_e_props[i]),dynamic_cast<Property<Vec3d>*>(m_e_props[i]),1.0e-9) )
      {
        if( printdiffs ) std::cout << "\033[31;1mWARNING EDGE PROPERTIES NOT EQUAL.\033[m" << std::endl;
        return false;
      }
      if( printdiffs ) std::cout << "EQUAL";
    }
    else if( dynamic_cast<Property<Util::pair<Scalar,Scalar> >*>(o_e_props[i]) )
    {
      if( !dynamic_cast<Property<Util::pair<Scalar, Scalar> >*>(m_e_props[i]) )
      {
        if( printdiffs ) std::cout << "\033[31;1mWARNING EDGE PROPERTIES OF SAME NAME NOT THE SAME TYPE.\033[m" << std::endl;
        return false;
      }
      if( !propertiesEqual(dynamic_cast<Property<Util::pair<Scalar, Scalar> >*>(o_e_props[i]),dynamic_cast<Property<Util::pair<Scalar, Scalar> >*>(m_e_props[i]),1.0e-9) )
      {
        if( printdiffs ) std::cout << "\033[31;1mWARNING EDGE PROPERTIES NOT EQUAL.\033[m" << std::endl;
        return false;
      }
      if( printdiffs ) std::cout << "EQUAL";
    }
    else if( dynamic_cast<Property< Util::pair<Vec3d,Vec3d> >* >(o_e_props[i]) )
    {
      if( !dynamic_cast<Property< Util::pair<Vec3d,Vec3d> >* >(m_e_props[i]) )
      {
        if( printdiffs ) std::cout << "\033[31;1mWARNING EDGE PROPERTIES OF SAME NAME NOT THE SAME TYPE.\033[m" << std::endl;
        return false;
      }
      if( !propertiesEqual(dynamic_cast<Property< Util::pair<Vec3d,Vec3d> >* >(o_e_props[i]),dynamic_cast<Property< Util::pair<Vec3d,Vec3d> >* >(m_e_props[i]),1.0e-9) )
      {
        if( printdiffs ) std::cout << "\033[31;1mWARNING EDGE PROPERTIES NOT EQUAL.\033[m" << std::endl;
        return false;
      }
      if( printdiffs ) std::cout << "EQUAL";
    }
    else if( dynamic_cast<Property<EdgeTopology<TopologicalObject> >*>(o_e_props[i]) )
    {
      if( !dynamic_cast<Property<EdgeTopology<TopologicalObject> >*>(m_e_props[i]) )
      {
        if( printdiffs ) std::cout << "\033[31;1mWARNING EDGE PROPERTIES OF SAME NAME NOT THE SAME TYPE.\033[m" << std::endl;
        return false;
      }
      if( !propertiesEqual(dynamic_cast<Property<EdgeTopology<TopologicalObject> >*>(o_e_props[i]),dynamic_cast<Property<EdgeTopology<TopologicalObject> >*>(m_e_props[i])) )
      {
        if( printdiffs ) std::cout << "\033[31;1mWARNING EDGE PROPERTIES NOT EQUAL.\033[m" << std::endl;
        return false;
      }
      if( printdiffs ) std::cout << "EQUAL";
    }
    else 
    {
      if( printdiffs ) std::cerr << "\033[31;1mERROR IN RODSTATEBACKUP:\033[m Unrecognized edge property type encountered for property: " << o_e_props[i]->name() << ". Exiting." << std::endl;
      exit(1);
    }
    if( printdiffs ) std::cout << std::endl;
  }
  
  if( printdiffs ) std::cout << "\033[36;1mCOMPARING FACE PROPERTIES\033[m" << std::endl;
  for( int i = 0; i < (int) o_f_props.size(); ++i )
  {
    assert( o_f_props[i] != NULL );
    assert( m_f_props[i] != NULL );
    
    if( o_f_props[i]->name() != m_f_props[i]->name() ) 
    {
      if( printdiffs ) std::cout << "\033[36;1mWARNING FACE PROPERTIES IN DIFFERENT ORDER.\033[m" << std::endl;
      return false;
    }
    
    if( printdiffs ) std::cout << "\033[36;1m  COMPARING: " << o_f_props[i]->name() << "\033[m: ";
    
    if( dynamic_cast<Property<FaceTopology<TopologicalObject> >*>(o_f_props[i]) )
    {
      if( !dynamic_cast<Property<FaceTopology<TopologicalObject> >*>(m_f_props[i]) )
      {
        if( printdiffs ) std::cout << "\033[31;1mWARNING FACE PROPERTIES OF SAME NAME NOT THE SAME TYPE.\033[m" << std::endl;
        return false;
      }
      if( !propertiesEqual(dynamic_cast<Property<FaceTopology<TopologicalObject> >*>(o_f_props[i]),dynamic_cast<Property<FaceTopology<TopologicalObject> >*>(m_f_props[i])) )
      {
        if( printdiffs ) std::cout << "\033[31;1mWARNING FACE PROPERTIES NOT EQUAL.\033[m" << std::endl;
        return false;
      }
      if( printdiffs ) std::cout << "EQUAL";
    }
    else 
    {
      if( printdiffs ) std::cerr << "\033[31;1mERROR IN RODSTATEBACKUP:\033[m Unrecognized face property type encountered for property: " << o_f_props[i]->name() << ". Exiting." << std::endl;
      exit(1);
    }
    if( printdiffs ) std::cout << std::endl;
  }
  
  return true;
}

const PropertyContainer& RodState::getVertexPropertyContainer( const ElasticRod& rod )
{
  VPropHandle<int> vh;
  return rod.container(vh);
}

PropertyContainer& RodState::getVertexPropertyContainer( ElasticRod& rod )
{
  VPropHandle<int> vh;
  return rod.container(vh);
}  

const PropertyContainer& RodState::getEdgePropertyContainer( const ElasticRod& rod )
{
  EPropHandle<int> eh;
  return rod.container(eh);
}

PropertyContainer& RodState::getEdgePropertyContainer( ElasticRod& rod )
{
  EPropHandle<int> eh;
  return rod.container(eh);
}  

const PropertyContainer& RodState::getFacePropertyContainer( const ElasticRod& rod )
{
  FPropHandle<int> fh;
  return rod.container(fh);
}

PropertyContainer& RodState::getFacePropertyContainer( ElasticRod& rod )
{
  FPropHandle<int> fh;
  return rod.container(fh);
}  

const PropertyContainer& RodState::getObjectPropertyContainer( const ElasticRod& rod )
{
  ObjPropHandle<int> oh;
  return rod.container(oh);
}

PropertyContainer& RodState::getObjectPropertyContainer( ElasticRod& rod )
{
  ObjPropHandle<int> oh;
  return rod.container(oh);
}






}
