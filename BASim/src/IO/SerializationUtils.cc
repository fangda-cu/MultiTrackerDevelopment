#include "SerializationUtils.hh"

namespace BASim
{

///////////////////////
// Integer Functions

void serializeVal( std::ofstream& of, const int& val )
{
  assert( of.is_open() );
  of.write((char*)&val,sizeof(int));
}

void loadVal( std::ifstream& ifs, int& val )
{
  assert( ifs.is_open() );
  ifs.read((char*)&val,sizeof(int));
}

void serializeProperty( std::ofstream& of, Property<int>* prop )
{
  for( size_t i = 0; i < prop->size(); ++i )
  {
    serializeVal(of, (*prop)[i]);
  }
}

void loadProperty( std::ifstream& ifs, Property<int>* prop )
{
  assert( ifs.is_open() );
  
  for( size_t i = 0; i < prop->size(); ++i )
  {
    int newint;
    loadVal(ifs,newint);
    (*prop)[i] = newint;
  }
}
  
///////////////////////
// Vec3d Functions
void printProperty( Property<Vec3d>* prp )
{
  for( int i  = 0; i < (int) prp->size(); ++i ) std::cout << (*prp)[i].transpose() << " ";
}

bool propertiesEqual( Property<Vec3d>* prpA, Property<Vec3d>* prpB, Scalar eps )
{
  if( prpA->size() != prpB->size() ) return false;
  for( int i = 0; i < (int) prpA->size(); ++i ) if( !approxEq((*prpA)[i],(*prpB)[i],eps) ) return false;
  return true;    
}

void serializeVec3d( std::ofstream& of, const Vec3d& val )
{
  assert( of.is_open() );
  for( int i = 0; i < val.size(); ++i )
  {
    Scalar x = val[i];
    of.write((char*)&x,sizeof(Scalar));
  }
}

void loadVec3d( std::ifstream& ifs, Vec3d& val )
{
  assert( ifs.is_open() );
  for( int i = 0; i < val.size(); ++i ) ifs.read((char*)&val[i],sizeof(Scalar));
}  

void serializeProperty( std::ofstream& of, Property<Vec3d>* prop )
{
  for( size_t i = 0; i < prop->size(); ++i ) serializeVec3d(of, (*prop)[i]);
}

void loadProperty( std::ifstream& ifs, Property<Vec3d>* prop )
{
  assert( ifs.is_open() );
  for( size_t i = 0; i < prop->size(); ++i )
  {
    Vec3d newval;
    loadVec3d(ifs,newval);
    (*prop)[i] = newval;
  }
}


//////////////////////////
// Scalar Functions
void printProperty( Property<Scalar>* prp )
{
  for( int i  = 0; i < (int) prp->size(); ++i ) std::cout << (*prp)[i] << " ";
}

bool propertiesEqual( Property<Scalar>* prpA, Property<Scalar>* prpB, Scalar eps )
{
  if( prpA->size() != prpB->size() ) return false;
  for( int i  = 0; i < (int) prpA->size(); ++i ) if( !approxEq((*prpA)[i],(*prpB)[i],eps) ) return false;
  return true;
}

void serializeScalar( std::ofstream& of, const Scalar& val )
{
  assert( of.is_open() );
  of.write((char*)&val,sizeof(Scalar));
}

void loadScalar( std::ifstream& ifs, Scalar& val )
{
  assert( ifs.is_open() );
  ifs.read((char*)&val,sizeof(Scalar));
}  
  
void serializeProperty( std::ofstream& of, Property<Scalar>* prop )
{
  for( size_t i = 0; i < prop->size(); ++i ) serializeScalar(of, (*prop)[i]);
}

void loadProperty( std::ifstream& ifs, Property<Scalar>* prop )
{
  assert( ifs.is_open() );

  for( size_t i = 0; i < prop->size(); ++i )
  {
    Scalar newval = std::numeric_limits<Scalar>::signaling_NaN();
    loadScalar(ifs,newval);
    (*prop)[i] = newval;
  }
}

//////////////////////////
// String Functions

void serializeString( std::ofstream& of, const std::string& name )
{
  assert( of.is_open() );
  serializeVal(of,(int)name.size());
  of.write((char*)name.c_str(),name.length()*sizeof(char));
}

void loadString( std::ifstream& ifs, std::string& name )
{
  assert( ifs.is_open() );
  int len;
  ifs.read((char*)&len,sizeof(int));
  assert( len >= 0 );    
  char* charname = new char[len];
  ifs.read(charname,len*sizeof(char));
  
  name.resize(len);
  for( int i = 0; i < len; ++i ) name[i] = charname[i];
  delete[] charname; charname = NULL;
}

//////////////////////////
// PropertyID Functions

void serializePropertyID( std::ofstream& of, const PropertyID& id )
{
  assert( of.is_open() );
  of.write((char*)&id,sizeof(PropertyID));
}  

//////////////////////////
// Force Functions
BASim::ForceID getForceID( RodForce* rforce )
{
  if( dynamic_cast<RodStretchingForce*>(rforce) )
  {
    if( !rforce->viscous() ) return RODSTRETCHINGFORCE;
    else return RODSTRETCHINGFORCEVISCOUS;
  }
  else if( dynamic_cast<RodTwistingForceSym*>(rforce) )
  {
    if( !rforce->viscous() ) return RODTWISTINGFORCESYM;
    else return RODTWISTINGFORCESYMVISCOUS;
  }
  else if( dynamic_cast<RodBendingForceSym*>(rforce) )
  {
    if( !rforce->viscous() ) return RODBENDINGFORCESYM;
    else return RODBENDINGFORCESYMVISCOUS;
  }
  
  std::cerr << "INVALID FORCE TYPE" << std::endl;
  exit(1);
  return RODSTRETCHINGFORCE;
}

RodForce* getRodForce( const ForceID& fid, ElasticRod* erod )
{
  switch (fid) 
  {
    case RODSTRETCHINGFORCE:
    {
      return new RodStretchingForce(*erod,false,false);
      break;
    }
    case RODSTRETCHINGFORCEVISCOUS:
    {
      return new RodStretchingForce(*erod,true,false);
      break;
    }
    case RODTWISTINGFORCESYM:
    {
      return new RodTwistingForceSym(*erod,false,false);
      break;
    }
    case RODTWISTINGFORCESYMVISCOUS:
    {
      return new RodTwistingForceSym(*erod,true,false);
      break;
    }
    case RODBENDINGFORCESYM:
    {
      return new RodBendingForceSym(*erod,false,false);
      break;
    }
    case RODBENDINGFORCESYMVISCOUS:
    {
      return new RodBendingForceSym(*erod,true,false);
      break;
    }
    default:
    {
      std::cerr << "INVALID FORCE TYPE" << std::endl;
      exit(1);
      break;
    }
  }
  
  std::cerr << "INVALID FORCE TYPE" << std::endl;
  exit(1);
  return NULL;
}

void serializeElasticRodRodForces( std::ofstream& of, const ElasticRod::RodForces& val )
{
  assert( of.is_open() );
  
  int sze = val.size();
  serializeVal(of,sze);
  
  for( size_t i = 0; i < val.size(); ++i )
  {
    ForceID fid = getForceID(val[i]);
    of.write((char*)&fid,sizeof(ForceID));
  }
}

void loadElasticRodRodForce( std::ifstream& ifs,  ElasticRod::RodForces& val, ElasticRod* erod )
{
  assert( ifs.is_open() );
  assert( erod != NULL );
  
  int sze = -1;
  loadVal(ifs,sze);
  assert( sze >= 0 );
  
  val.resize(sze);
  for( int i = 0; i < sze; ++i )
  {
    ForceID fid;
    ifs.read((char*)&fid,sizeof(ForceID));
    RodForce* newforce = getRodForce(fid,erod);
    assert( newforce != NULL );
    val[i] = newforce;
  }
  assert( (int) val.size() == sze );
}

void serializeProperty( std::ofstream& of, Property<ElasticRod::RodForces>* prop )
{
  for( size_t i = 0; i < prop->size(); ++i ) serializeElasticRodRodForces(of, (*prop)[i]);
}

void loadProperty( std::ifstream& ifs, Property<ElasticRod::RodForces>* prop, ElasticRod* erod )
{
  assert( ifs.is_open() );
  
  for( size_t i = 0; i < prop->size(); ++i )
  {
    ElasticRod::RodForces newval;
    loadElasticRodRodForce(ifs,newval,erod);
    #ifdef DEBUG
      for( size_t j = 0; j < newval.size(); ++j ) assert( newval[j] != NULL );
    #endif
    (*prop)[i] = newval;
  }
}
  
//////////////////////////
// Bool Functions

void serializeBool( std::ofstream& of, const bool& val )
{
  assert( of.is_open() );
  of.write((char*)&val,sizeof(bool));
}

void loadBool( std::ifstream& ifs, bool& val )
{
  assert( ifs.is_open() );
  ifs.read((char*)&val,sizeof(bool));
}

void serializeProperty( std::ofstream& of, Property<bool>* prop )
{
  for( size_t i = 0; i < prop->size(); ++i ) serializeBool(of, (*prop)[i]);
}

void loadProperty( std::ifstream& ifs, Property<bool>* prop )
{
  assert( ifs.is_open() );
  
  for( size_t i = 0; i < prop->size(); ++i )
  {
    bool newval;
    loadBool(ifs,newval);
    (*prop)[i] = newval;
  }
}

//////////////////////////
// ElasticRodRefFrame Functions

void serializeElasticRodRefFrameType( std::ofstream& of, const ElasticRod::RefFrameType& val )
{
  assert( of.is_open() );
  of.write((char*)&val,sizeof(ElasticRod::RefFrameType));
}

void loadElasticRodRefFrameType( std::ifstream& ifs, ElasticRod::RefFrameType& val )
{
  assert( ifs.is_open() );
  ifs.read((char*)&val,sizeof(ElasticRod::RefFrameType));
}  
  
void serializeProperty( std::ofstream& of, Property<ElasticRod::RefFrameType>* prop )
{
  for( size_t i = 0; i < prop->size(); ++i ) serializeElasticRodRefFrameType(of, (*prop)[i]);
}

void loadProperty( std::ifstream& ifs, Property<ElasticRod::RefFrameType>* prop )
{
  assert( ifs.is_open() );
  for( size_t i = 0; i < prop->size(); ++i )
  {
    ElasticRod::RefFrameType newval;
    loadElasticRodRefFrameType(ifs,newval);
    (*prop)[i] = newval;
  }
}

//////////////////////////
// RodBoundaryConditionBCList Functions
void serializeRodBoundaryConditionBCList( std::ofstream& of, const RodBoundaryCondition::BCList& val )
{
  assert( of.is_open() );
  
  int sze = val.size();
  serializeVal(of,sze);
  
  for( size_t i = 0; i < val.size(); ++i ) serializeVal(of,val[i]);
}

void loadRodBoundaryConditionBCList( std::ifstream& ifs, RodBoundaryCondition::BCList& val )
{
  assert( ifs.is_open() );

  int sze = -1;
  loadVal(ifs,sze);
  assert( sze >= 0 );

  val.resize(sze);
  for( int i = 0; i < sze; ++i ) loadVal(ifs,val[i]);
}
  
void serializeProperty( std::ofstream& of, Property<RodBoundaryCondition::BCList>* prop )
{
  for( size_t i = 0; i < prop->size(); ++i ) serializeRodBoundaryConditionBCList(of, (*prop)[i]);
}

void loadProperty( std::ifstream& ifs, Property<RodBoundaryCondition::BCList>* prop )
{
  assert( ifs.is_open() );
  
  for( size_t i = 0; i < prop->size(); ++i )
  {
    RodBoundaryCondition::BCList newval;
    loadRodBoundaryConditionBCList(ifs,newval);
    (*prop)[i] = newval;
  }
}  

//////////////////////////
// DofHandleType Functions
void serializeDofHandleType( std::ofstream& of, const DofHandle::Type& val )
{
  assert( of.is_open() );
  of.write((char*)&val,sizeof(DofHandle::Type));
}

void loadDofHandleType( std::ifstream& ifs, DofHandle::Type& val )
{
  assert( ifs.is_open() );
  ifs.read((char*)&val,sizeof(DofHandle::Type));
}
  
//////////////////////////
// DofHandle Functions
void serializeDofHandle( std::ofstream& of, const DofHandle& val )
{
  assert( of.is_open() );
  
  DofHandle::Type dhtype = val.getType();
  int memberhndlidx = val.getHandle().idx();
  int parenthndlidx = val.idx();
  int dofnum = val.getNum();
  
  serializeDofHandleType(of,dhtype);
  serializeVal(of,memberhndlidx);
  serializeVal(of,parenthndlidx);
  serializeVal(of,dofnum);
}

void loadDofHandle( std::ifstream& ifs, DofHandle& val )
{
  assert( ifs.is_open() );
  
  DofHandle::Type dhtype;
  int memberhndlidx;
  int parenthndlidx;
  int dofnum;
  
  loadDofHandleType(ifs,dhtype);
  loadVal(ifs,memberhndlidx);
  loadVal(ifs,parenthndlidx);
  loadVal(ifs,dofnum);
  
  val.setType(dhtype);
  val.setHandle(HandleBase(memberhndlidx));
  val.idx() = parenthndlidx;
  val.setNum(dofnum);
}
  
//////////////////////////
// DOFMap Functions
void serializeDOFMap( std::ofstream& of, DOFMap& val )
{
  assert( of.is_open() );
  
  int sze = (int) val.size();
  serializeVal(of,sze);
  
  for( std::map<DofHandle, int>::const_iterator itr = val.getDofToIndexMap().begin(); itr != val.getDofToIndexMap().end(); ++itr )
  {
    serializeDofHandle(of,itr->first);
    serializeVal(of,itr->second);
  }
}

void loadDOFMap( std::ifstream& ifs, DOFMap& val )
{
  assert( ifs.is_open() );
  
  int sze = -1;
  loadVal(ifs,sze);
  
  val.clearMappings();
  assert( val.getDofToIndexMap().size() == 0 );
  assert( val.getIndexToDofMap().size() == 0 );
  
  for( int i = 0; i < sze; ++i )
  {
    assert( val.getDofToIndexMap().size() == val.getIndexToDofMap().size() );
    
    DofHandle dh;
    loadDofHandle(ifs,dh);
    int idx;
    loadVal(ifs,idx);
    
    val.addMapping(dh,idx);
  }
  
  assert( (int) val.size() == sze );
}
  
void serializeProperty( std::ofstream& of, Property<DOFMap>* prop )
{
  for( size_t i = 0; i < prop->size(); ++i ) serializeDOFMap(of, (*prop)[i]);
}

void loadProperty( std::ifstream& ifs, Property<DOFMap>* prop )
{
  assert( ifs.is_open() );
  
  for( size_t i = 0; i < prop->size(); ++i )
  {
    DOFMap newval;
    loadDOFMap(ifs,newval);
    (*prop)[i] = newval;
  }
}
  
//////////////////////////
// Vec2d Functions
void serializeVec2d( std::ofstream& of, const Vec2d& val )
{
  assert( of.is_open() );
  for( int i = 0; i < val.size(); ++i )
  {
    Scalar x = val[i];
    of.write((char*)&x,sizeof(Scalar));
  }
}

void loadVec2d( std::ifstream& ifs, Vec2d& val )
{
  assert( ifs.is_open() );
  for( int i = 0; i < val.size(); ++i ) ifs.read((char*)&val[i],sizeof(Scalar));
}

void serializeProperty( std::ofstream& of, Property<Vec2d>* prop )
{
  for( size_t i = 0; i < prop->size(); ++i ) serializeVec2d(of, (*prop)[i]);
}

void loadProperty( std::ifstream& ifs, Property<Vec2d>* prop )
{
  assert( ifs.is_open() );
  
  for( size_t i = 0; i < prop->size(); ++i )
  {
    Vec2d newval;
    loadVec2d(ifs,newval);
    (*prop)[i] = newval;
  }
}

//////////////////////////
// VecXd Functions
void serializeVecXd( std::ofstream& of, const VecXd& val )
{
  assert( of.is_open() );
  int len = val.size();
  of.write((char*)&len,sizeof(int));
  for( int i = 0; i < val.size(); ++i )
  {
    Scalar x = val[i];
    of.write((char*)&x,sizeof(Scalar));
  }
}

void loadVecXd( std::ifstream& ifs, VecXd& val )
{
  assert( ifs.is_open() );
  int len = -1;
  ifs.read((char*)&len,sizeof(int));
  assert( len >= 0 );
  val.resize(len);
  for( int i = 0; i < val.size(); ++i ) ifs.read((char*)&val[i],sizeof(Scalar));
}  
  
void serializeProperty( std::ofstream& of, Property<VecXd>* prop )
{
  for( size_t i = 0; i < prop->size(); ++i ) serializeVecXd(of, (*prop)[i]);
}

void loadProperty( std::ifstream& ifs, Property<VecXd>* prop )
{
  assert( ifs.is_open() );
  
  for( size_t i = 0; i < prop->size(); ++i )
  {
    VecXd newval;
    loadVecXd(ifs,newval);
    (*prop)[i] = newval;
  }
}
  
//////////////////////////
// Mat2d Functions
void serializeMat2d( std::ofstream& of, const Mat2d& val )
{
  assert( of.is_open() );
  for( int i = 0; i < val.rows(); ++i )
  {
    for( int j = 0; j < val.cols(); ++j )
    {
      Scalar x = val(i,j);
      of.write((char*)&x,sizeof(Scalar));
    }
  }
}

void loadMat2d( std::ifstream& ifs, Mat2d& val )
{
  assert( ifs.is_open() );
  for( int i = 0; i < val.rows(); ++i ) for( int j = 0; j < val.cols(); ++j ) ifs.read((char*)&val(i,j),sizeof(Scalar));
}

void serializeProperty( std::ofstream& of, Property<Mat2d>* prop )
{
  for( size_t i = 0; i < prop->size(); ++i ) serializeMat2d(of, (*prop)[i]);
}

void loadProperty( std::ifstream& ifs, Property<Mat2d>* prop )
{
  assert( ifs.is_open() );
  
  for( size_t i = 0; i < prop->size(); ++i )
  {
    Mat2d newval;
    loadMat2d(ifs,newval);
    (*prop)[i] = newval;
  }
}
  
//////////////////////////
// MatXd Functions

void serializeMatXd( std::ofstream& of, const MatXd& val )
{
  assert( of.is_open() );
  int rows = val.rows();
  int cols = val.cols();
  of.write((char*)&rows,sizeof(int));
  of.write((char*)&cols,sizeof(int));
  for( int i = 0; i < val.rows(); ++i )
  {
    for( int j = 0; j < val.cols(); ++j )
    {
      Scalar x = val(i,j);
      of.write((char*)&x,sizeof(Scalar));
    }
  }
}

void loadMatXd( std::ifstream& ifs, MatXd& val )
{
  assert( ifs.is_open() );
  int rows = -1;
  int cols = -1;
  ifs.read((char*)&rows,sizeof(int));
  ifs.read((char*)&cols,sizeof(int));
  assert( rows >= 0 );
  assert( cols >= 0 );
  val.resize(rows,cols);
  for( int i = 0; i < val.rows(); ++i )
  {
    for( int j = 0; j < val.cols(); ++j )
    {
      ifs.read((char*)&val(i,j),sizeof(Scalar));
    }
  }
}  

void serializeProperty( std::ofstream& of, Property<MatXd>* prop )
{
  for( size_t i = 0; i < prop->size(); ++i )
  {
    serializeMatXd(of, (*prop)[i]);
  }
}

void loadProperty( std::ifstream& ifs, Property<MatXd>* prop )
{
  assert( ifs.is_open() );
  
  for( size_t i = 0; i < prop->size(); ++i )
  {
    MatXd newval;
    loadMatXd(ifs,newval);
    (*prop)[i] = newval;
  }
}

//////////////////////////
// std::pair<Scalar,Scalar> Functions
void serializePairScalarScalar( std::ofstream& of, const std::pair<Scalar,Scalar>& val )
{
  assert( of.is_open() );
  serializeScalar(of,val.first);
  serializeScalar(of,val.second);
}

void loadPairScalarScalar( std::ifstream& ifs, std::pair<Scalar,Scalar>& val )
{
  assert( ifs.is_open() );
  loadScalar(ifs,val.first);
  loadScalar(ifs,val.second);
}

//////////////////////////
// std::vector< std::pair<Scalar,Scalar> > Functions
void serializeVectorPairScalar( std::ofstream& of, const std::vector< std::pair<Scalar,Scalar> >& val )
{
  assert( of.is_open() );
  int len = val.size();
  serializeVal(of,len);
  for( size_t i = 0; i < val.size(); ++i )
  {
    serializePairScalarScalar(of,val[i]);
  }
}

void loadVectorPairScalar( std::ifstream& ifs, std::vector< std::pair<Scalar,Scalar> >& val )
{
  assert( ifs.is_open() );
  int len;
  loadVal(ifs,len);
  val.resize(len);
  for( size_t i = 0; i < val.size(); ++i )
  {
    loadPairScalarScalar(ifs,val[i]);
  }
}

void serializeProperty( std::ofstream& of, Property<std::vector< std::pair<Scalar,Scalar> > >* prop )
{
  for( size_t i = 0; i < prop->size(); ++i )
  {
    serializeVectorPairScalar(of, (*prop)[i]);
  }
}

void loadProperty( std::ifstream& ifs, Property<std::vector< std::pair<Scalar,Scalar> > >* prop )
{
  assert( ifs.is_open() );
  
  for( size_t i = 0; i < prop->size(); ++i )
  {
    std::vector< std::pair<Scalar,Scalar> > newval;
    loadVectorPairScalar(ifs,newval);
    (*prop)[i] = newval;
  }
}

//////////////////////////
// std::vector<Vec3d> Functions
void serializeVectorVec3d( std::ofstream& of, const std::vector<Vec3d>& val )
{
  assert( of.is_open() );
  int len = val.size();
  serializeVal(of,len);
  for( size_t i = 0; i < val.size(); ++i )
  {
    serializeVec3d(of,val[i]);
  }
}

void loadVectorVec3d( std::ifstream& ifs, std::vector<Vec3d>& val )
{
  assert( ifs.is_open() );
  int len;
  loadVal(ifs,len);
  val.resize(len);
  for( size_t i = 0; i < val.size(); ++i )
  {
    loadVec3d(ifs,val[i]);
  }
}

void serializeProperty( std::ofstream& of, Property<std::vector<Vec3d> >* prop )
{
  for( size_t i = 0; i < prop->size(); ++i )
  {
    serializeVectorVec3d(of, (*prop)[i]);
  }
}

void loadProperty( std::ifstream& ifs, Property<std::vector<Vec3d> >* prop )
{
  assert( ifs.is_open() );
  
  for( size_t i = 0; i < prop->size(); ++i )
  {
    std::vector<Vec3d> newval;
    loadVectorVec3d(ifs,newval);
    (*prop)[i] = newval;
  }
}
  
//////////////////////////
// std::pair<MatXd,MatXd> Functions
void serializePairMatXdMatXd( std::ofstream& of, const std::pair<MatXd,MatXd>& val )
{
  assert( of.is_open() );
  serializeMatXd(of,val.first);
  serializeMatXd(of,val.second);
}

void loadPairMatXdMatXd( std::ifstream& ifs, std::pair<MatXd,MatXd>& val )
{
  assert( ifs.is_open() );
  loadMatXd(ifs,val.first);
  loadMatXd(ifs,val.second);
}

void serializeProperty( std::ofstream& of, Property<std::pair<MatXd,MatXd> >* prop )
{
  for( size_t i = 0; i < prop->size(); ++i )
  {
    serializePairMatXdMatXd(of, (*prop)[i]);
  }
}

void loadProperty( std::ifstream& ifs, Property<std::pair<MatXd,MatXd> >* prop )
{
  assert( ifs.is_open() );
  
  for( size_t i = 0; i < prop->size(); ++i )
  {
    std::pair<MatXd,MatXd> newval;
    loadPairMatXdMatXd(ifs,newval);
    (*prop)[i] = newval;
  }
}
  
//////////////////////////
// Util::pair<Scalar,Scalar> Functions
void serializePairScalarScalar( std::ofstream& of, const Util::pair<Scalar,Scalar>& val )
{
  assert( of.is_open() );
  serializeScalar(of,val.first);
  serializeScalar(of,val.second);
}

void loadPairScalarScalar( std::ifstream& ifs, Util::pair<Scalar,Scalar>& val )
{
  assert( ifs.is_open() );
  loadScalar(ifs,val.first);
  loadScalar(ifs,val.second);
}

void serializeProperty( std::ofstream& of, Property<Util::pair<Scalar,Scalar> >* prop )
{
  for( size_t i = 0; i < prop->size(); ++i )
  {
    serializePairScalarScalar(of, (*prop)[i]);
  }
}

void loadProperty( std::ifstream& ifs, Property<Util::pair<Scalar,Scalar> >* prop )
{
  assert( ifs.is_open() );
  
  for( size_t i = 0; i < prop->size(); ++i )
  {
    Util::pair<Scalar,Scalar> newval;
    loadPairScalarScalar(ifs,newval);
    (*prop)[i] = newval;
  }
}

//////////////////////////
// Util::pair<Vec3d,Vec3d> Functions
void serializePairVec3dVec3d( std::ofstream& of, const Util::pair<Vec3d,Vec3d>& val )
{
  assert( of.is_open() );
  serializeVec3d(of,val.first);
  serializeVec3d(of,val.second);
}

void loadPairVec3dVec3d( std::ifstream& ifs, Util::pair<Vec3d,Vec3d>& val )
{
  assert( ifs.is_open() );
  loadVec3d(ifs,val.first);
  loadVec3d(ifs,val.second);
}

void serializeProperty( std::ofstream& of, Property<Util::pair<Vec3d,Vec3d> >* prop )
{
  for( size_t i = 0; i < prop->size(); ++i )
  {
    serializePairVec3dVec3d(of, (*prop)[i]);
  }
}

void loadProperty( std::ifstream& ifs, Property<Util::pair<Vec3d,Vec3d> >* prop )
{
  assert( ifs.is_open() );
  
  for( size_t i = 0; i < prop->size(); ++i )
  {
    Util::pair<Vec3d,Vec3d> newval;
    loadPairVec3dVec3d(ifs,newval);
    (*prop)[i] = newval;
  }
}

//////////////////////////
// VertexTopology<TopologicalObject> Functions
void serializeVeretexTopology( std::ofstream& of, const VertexTopology<TopologicalObject>& val )
{
  assert( of.is_open() );
  int len = val.size();
  serializeVal(of,len);
  for( size_t i = 0; i < val.size(); ++i )
  {
    serializeVal(of,val[i].idx());
  }
}  

void loadVeretexTopology( std::ifstream& ifs, VertexTopology<TopologicalObject>& val )
{
  assert( ifs.is_open() );
  int len;
  loadVal(ifs,len);
  val.resize(len);
  for( size_t i = 0; i < val.size(); ++i )
  {
    loadVal(ifs,val[i].idx());
  }
}
  
void serializeProperty( std::ofstream& of, Property<VertexTopology<TopologicalObject> >* prop )
{
  for( size_t i = 0; i < prop->size(); ++i )
  {
    serializeVeretexTopology(of, (*prop)[i]);
  }
}

void loadProperty( std::ifstream& ifs, Property<VertexTopology<TopologicalObject> >* prop )
{
  assert( ifs.is_open() );
  
  for( size_t i = 0; i < prop->size(); ++i )
  {
    VertexTopology<TopologicalObject> newval;
    loadVeretexTopology(ifs,newval);
    (*prop)[i] = newval;
  }
}
  
//////////////////////////
// EdgeTopology<TopologicalObject>
void serializeEdgeTopology( std::ofstream& of, const EdgeTopology<TopologicalObject>& val )
{
  assert( of.is_open() );
  assert( val.size() == 2 );
  //int len = val.size();
  //serializeVal(of,len);
  for( size_t i = 0; i < val.size(); ++i )
  {
    serializeVal(of,val[i].idx());
  }
}  

void loadEdgeTopology( std::ifstream& ifs, EdgeTopology<TopologicalObject>& val )
{
  assert( ifs.is_open() );
  assert( val.size() == 2 );
  //int len;
  //loadVal(ifs,len);
  //val.resize(len);
  for( size_t i = 0; i < val.size(); ++i )
  {
    loadVal(ifs,val[i].idx());
  }
}

void serializeProperty( std::ofstream& of, Property<EdgeTopology<TopologicalObject> >* prop )
{
  for( size_t i = 0; i < prop->size(); ++i )
  {
    serializeEdgeTopology(of, (*prop)[i]);
  }
}

void loadProperty( std::ifstream& ifs, Property<EdgeTopology<TopologicalObject> >* prop )
{
  assert( ifs.is_open() );
  
  for( size_t i = 0; i < prop->size(); ++i )
  {
    EdgeTopology<TopologicalObject> newval;
    loadEdgeTopology(ifs,newval);
    (*prop)[i] = newval;
  }
}

//////////////////////////
// FaceTopology<TopologicalObject> Functions
void serializeFaceTopology( std::ofstream& of, const FaceTopology<TopologicalObject>& val )
{
  assert( of.is_open() );
  int len = val.size();
  serializeVal(of,len);
  for( size_t i = 0; i < val.size(); ++i )
  {
    serializeVal(of,val[i].idx());
  }
}  

void loadFaceTopology( std::ifstream& ifs, FaceTopology<TopologicalObject>& val )
{
  assert( ifs.is_open() );
  int len;
  loadVal(ifs,len);
  val.resize(len);
  for( size_t i = 0; i < val.size(); ++i )
  {
    loadVal(ifs,val[i].idx());
  }
}

void serializeProperty( std::ofstream& of, Property<FaceTopology<TopologicalObject> >* prop )
{
  assert( of.is_open() );
  for( size_t i = 0; i < prop->size(); ++i ) serializeFaceTopology(of, (*prop)[i]);
}

void loadProperty( std::ifstream& ifs, Property<FaceTopology<TopologicalObject> >* prop )
{
  assert( ifs.is_open() );
  
  for( size_t i = 0; i < prop->size(); ++i )
  {
    FaceTopology<TopologicalObject> newval;
    loadFaceTopology(ifs,newval);
    (*prop)[i] = newval;
  }
}

//////////////////////////
// std::queue<double> Functions
void serializeDoubleQueue( std::ofstream& of, const std::queue<double>& val )
{
  assert( of.is_open() );

  // Kind of a hack since I can't iterate over a queue's elements
  std::queue<double> queuecopy = val;

  // Save the number of elements in the queue
  serializeVal(of,queuecopy.size());
  // Save the values in the queue
  while( !queuecopy.empty() )
  {
    serializeScalar(of,queuecopy.front());
    queuecopy.pop();
  }
}

void loadDoubleQueue( std::ifstream& ifs, std::queue<double>& val )
{
  assert( ifs.is_open() );
  
  // Load the number of elements in the queue
  int size = -1;
  loadVal(ifs,size);
  assert( size >= 0 );
  // Load the values in the queue
  for( int i = 0; i < size; ++i )
  {
    double queueval = std::numeric_limits<double>::signaling_NaN();
    loadScalar(ifs,queueval);
    val.push(queueval);
  }
  assert( (int) val.size() == size );
}

//////////////////////////
// Option Functions

void serializeOptionType( std::ofstream& of, const Option::Type& val )
{
  assert( of.is_open() );
  of.write((char*)&val,sizeof(Option::Type));
}

void loadOptionType( std::ifstream& ifs, Option::Type& val )
{
  assert( ifs.is_open() );
  ifs.read((char*)&val,sizeof(Option::Type));
}
  
void serializeOption( std::ofstream& of, const Option& val )
{
  assert( of.is_open() );

  // std::string name
  serializeString(of,val.name);
  // Type type
  serializeOptionType(of,val.type);
  // bool b
  serializeBool(of,val.b);
  // int i
  serializeVal(of,val.i);
  // Scalar r
  serializeScalar(of,val.r);
  // Vec3d v
  serializeVec3d(of,val.v);
  // std::string s
  serializeString(of,val.s);
  // std::string label
  serializeString(of,val.label);
}

void loadOption( std::ifstream& ifs, Option& val )
{
  assert( ifs.is_open() );

  // std::string name
  loadString(ifs,val.name);
  // Type type
  loadOptionType(ifs,val.type);
  // bool b
  loadBool(ifs,val.b);
  // int i
  loadVal(ifs,val.i);
  // Scalar r
  loadScalar(ifs,val.r);
  // Vec3d v
  loadVec3d(ifs,val.v);
  // std::string s
  loadString(ifs,val.s);
  // std::string label
  loadString(ifs,val.label);  
}

void serializeMapStringOption( std::ofstream& of, const std::map<std::string,Option>& val )
{
  assert( of.is_open() );
  
  // Serialize the number of options
  int size = val.size();
  serializeVal(of,size);
  //std::cout << "Options saved size: " << size << std::endl;
  // Serialize each string/option combo
  for( std::map<std::string,Option>::const_iterator itr = val.begin(); itr != val.end(); ++itr )
  {
    serializeString(of,itr->first);
    serializeOption(of,itr->second);
  }
}

void loadMapStringOption( std::ifstream& ifs, std::map<std::string,Option>& val )
{
  assert( ifs.is_open() );
  
  // Load the number of options
  int size = (int)std::numeric_limits<double>::signaling_NaN();
  loadVal(ifs,size);
  assert( size >= 0 );
  //std::cout << "Size: " << size << std::endl;
  // Load each string/option combo
  for( int i = 0; i < size; ++i )
  {
    std::string optionname;
    loadString(ifs,optionname);
    Option optn;
    loadOption(ifs,optn);
    val.insert(std::pair<std::string,Option>(optionname,optn));
  }  
  assert( size == (int) val.size() );
}





  
  
  
  
void serializeProperty( std::ofstream& of, PropertyBase* prop )
{
  assert( of.is_open() );
  
  // Save the name of the property
  serializeString(of,prop->name());
  
  // Save the property itself
  if( dynamic_cast<Property<ElasticRod::RodForces>*>(prop) )
  {
    serializePropertyID(of,ELASTICRODRODFORCES);
    Property<ElasticRod::RodForces>* chldprp = dynamic_cast<Property<ElasticRod::RodForces>*>(prop);
    serializeProperty( of, chldprp ); 
  }    
  else if( dynamic_cast<Property<int>*>(prop) )
  {
    serializePropertyID(of,INT);
    Property<int>* chldprp = dynamic_cast<Property<int>*>(prop);
    serializeProperty( of, chldprp ); 
  }
  else if( dynamic_cast<Property<bool>*>(prop) )
  {
    serializePropertyID(of,BOOL);
    Property<bool>* chldprp = dynamic_cast<Property<bool>*>(prop);
    serializeProperty( of, chldprp ); 
  }
  else if( dynamic_cast<Property<Scalar>*>(prop) )
  {
    serializePropertyID(of,SCALAR);
    Property<Scalar>* chldprp = dynamic_cast<Property<Scalar>*>(prop);
    serializeProperty( of, chldprp );
  }
  else if( dynamic_cast<Property<ElasticRod::RefFrameType>*>(prop) )
  {
    serializePropertyID(of,ELASTICRODREFFRAMETYPE);
    Property<ElasticRod::RefFrameType>* chldprp = dynamic_cast<Property<ElasticRod::RefFrameType>*>(prop);
    serializeProperty( of, chldprp );
  }      
  else if( dynamic_cast<Property<RodBoundaryCondition::BCList>*>(prop) )
  {
    serializePropertyID(of,RODBOUNDARYCONDITIONBCLIST);
    Property<RodBoundaryCondition::BCList>* chldprp = dynamic_cast<Property<RodBoundaryCondition::BCList>*>(prop);
    serializeProperty( of, chldprp );
  }
  else if( dynamic_cast<Property<DOFMap>*>(prop) )
  {
    serializePropertyID(of,DOFMAP);
    Property<DOFMap>* chldprp = dynamic_cast<Property<DOFMap>*>(prop);
    serializeProperty( of, chldprp );
  }
  else if( dynamic_cast<Property<Vec2d>*>(prop) )
  {
    serializePropertyID(of,VEC2D);
    Property<Vec2d>* chldprp = dynamic_cast<Property<Vec2d>*>(prop);
    serializeProperty( of, chldprp );
  }
  else if( dynamic_cast<Property<Vec3d>*>(prop) )
  {
    serializePropertyID(of,VEC3D);
    Property<Vec3d>* chldprp = dynamic_cast<Property<Vec3d>*>(prop);
    serializeProperty( of, chldprp );
  }
  else if( dynamic_cast<Property<VecXd>*>(prop) )
  {
    serializePropertyID(of,VECXD);
    Property<VecXd>* chldprp = dynamic_cast<Property<VecXd>*>(prop);
    serializeProperty( of, chldprp );
  }
  else if( dynamic_cast<Property<Mat2d>*>(prop) )
  {
    serializePropertyID(of,MAT2D);
    Property<Mat2d>* chldprp = dynamic_cast<Property<Mat2d>*>(prop);
    serializeProperty( of, chldprp );
  }
  else if( dynamic_cast<Property<MatXd>*>(prop) )
  {
    serializePropertyID(of,MATXD);
    Property<MatXd>* chldprp = dynamic_cast<Property<MatXd>*>(prop);
    serializeProperty( of, chldprp );
  }
  else if( dynamic_cast<Property<std::vector< std::pair<Scalar,Scalar> > >*>(prop) )
  {
    serializePropertyID(of,VECTORPAIRSCALARSCALAR);
    Property<std::vector< std::pair<Scalar,Scalar> > >* chldprp = dynamic_cast<Property<std::vector< std::pair<Scalar,Scalar> > >*>(prop);
    serializeProperty( of, chldprp );
  }
  else if( dynamic_cast<Property< std::vector<Vec3d> >*>(prop) )
  {
    serializePropertyID(of,VECTORVEC3D);
    Property<std::vector<Vec3d> >* chldprp = dynamic_cast<Property<std::vector<Vec3d> >*>(prop);
    serializeProperty( of, chldprp );
  }
  else if( dynamic_cast<Property< std::pair<MatXd,MatXd> >*>(prop) )
  {
    serializePropertyID(of,PAIRMATXDMATXD);
    Property<std::pair<MatXd,MatXd> >* chldprp = dynamic_cast<Property<std::pair<MatXd,MatXd> >*>(prop);
    serializeProperty( of, chldprp );
  }
  else if( dynamic_cast<Property<Util::pair<Scalar,Scalar> >*>(prop) )
  {
    serializePropertyID(of,PAIRSCALARSCALAR);
    Property<Util::pair<Scalar,Scalar> >* chldprp = dynamic_cast<Property<Util::pair<Scalar,Scalar> >*>(prop);
    serializeProperty( of, chldprp );
  }
  else if( dynamic_cast<Property< Util::pair<Vec3d,Vec3d> >* >(prop) )
  {
    serializePropertyID(of,PAIRVEC3DVEC3D);
    Property<Util::pair<Vec3d,Vec3d> >* chldprp = dynamic_cast<Property<Util::pair<Vec3d,Vec3d> >*>(prop);
    serializeProperty( of, chldprp );
  }
  else if( dynamic_cast<Property<VertexTopology<TopologicalObject> >*>(prop) )
  {
    serializePropertyID(of,VERTEXTOPOLOGY);
    Property<VertexTopology<TopologicalObject> >* chldprp = dynamic_cast<Property<VertexTopology<TopologicalObject> >*>(prop);
    serializeProperty( of, chldprp );
  }
  else if( dynamic_cast<Property<EdgeTopology<TopologicalObject> >*>(prop) )
  {
    serializePropertyID(of,EDGETOPOLOGY);
    Property<EdgeTopology<TopologicalObject> >* chldprp = dynamic_cast<Property<EdgeTopology<TopologicalObject> >*>(prop);
    serializeProperty( of, chldprp );
  }
  else if( dynamic_cast<Property<FaceTopology<TopologicalObject> >*>(prop) )
  {
    serializePropertyID(of,FACETOPOLOGY);
    Property<FaceTopology<TopologicalObject> >* chldprp = dynamic_cast<Property<FaceTopology<TopologicalObject> >*>(prop);
    serializeProperty( of, chldprp );
  }
  else 
  {
    std::cerr << "\033[31;1mERROR IN SERIALIZATIONUTILS:\033[m Attempt to serialize invalid property type. Exiting." << std::endl;
    exit(1);
  }
}  

void loadProperty( std::ifstream& ifs, PropertyBase** prop, ObjectBase* topobj, int propsize )
{
  assert( ifs.is_open() );
  
  // Load the property name
  std::string propertyname;
  loadString( ifs, propertyname );
  
  // Load the unique property id
  PropertyID pid;
  ifs.read((char*)&pid,sizeof(PropertyID));
  
  bool print_loaded = false;
  
  switch (pid) 
  {
    case ELASTICRODRODFORCES:
    {
      ElasticRod* erod = dynamic_cast<ElasticRod*>(topobj);
      if( !erod ) 
      {
        std::cerr << "\033[31;1mERROR IN TOPOLOGICALOBJECTSERIALIZER:\033[m Attempt to load non-rod object into a rod. Exiting." << std::endl;
        exit(1);
      }
      
      Property<ElasticRod::RodForces>* newprop = new Property<ElasticRod::RodForces>(propertyname);
      newprop->resize(propsize);
      BASim::loadProperty(ifs, newprop, erod);
      *prop = newprop;
      if( print_loaded )
      {
        std::cout << "   ELASTICRODRODFORCES" << " : " << propsize << " : " << propertyname << std::endl;
        for( size_t i = 0; i < newprop->size(); ++i )
        {
          std::cout << "      ";
          for( size_t j = 0; j < (*newprop)[i].size(); ++j )
          {
            assert( (*newprop)[i][j] != NULL );
            std::cout << (*newprop)[i][j]->getName() << " ";
          }
          std::cout << std::endl;
        }
      }
      break;
    }        
    case INT:
    {
      Property<int>* newprop = new Property<int>(propertyname);
      newprop->resize(propsize);
      BASim::loadProperty(ifs, newprop);
      *prop = newprop;
      if( print_loaded )
      {
        std::cout << "   INT" << " : " << propsize << " : " << propertyname << std::endl;
        for( size_t i = 0; i < newprop->size(); ++i ) std::cout << "     " << (*newprop)[i] << " "; std::cout << std::endl;
      }
      break;
    }
    case BOOL:
    {
      Property<bool>* newprop = new Property<bool>(propertyname);
      newprop->resize(propsize);
      BASim::loadProperty(ifs, newprop);
      *prop = newprop;
      if( print_loaded )
      {
        std::cout << "   BOOL"  << " : " << propsize << " : " << propertyname << std::endl;
        for( size_t i = 0; i < newprop->size(); ++i ) std::cout << "     " << (*newprop)[i] << " "; std::cout << std::endl;
      }
      break;
    }
    case SCALAR:
    {
      Property<Scalar>* newprop = new Property<Scalar>(propertyname);
      newprop->resize(propsize);
      BASim::loadProperty(ifs, newprop);
      *prop = newprop;
      if( print_loaded )
      {
        std::cout << "   SCALAR" << " : " << propsize << " : " << propertyname << std::endl;
        for( size_t i = 0; i < newprop->size(); ++i ) std::cout << "     " << (*newprop)[i] << " "; std::cout << std::endl;
      }
      break;
    }
    case ELASTICRODREFFRAMETYPE:
    {
      Property<ElasticRod::RefFrameType>* newprop = new Property<ElasticRod::RefFrameType>(propertyname);
      newprop->resize(propsize);
      BASim::loadProperty(ifs, newprop);
      *prop = newprop;        
      if( print_loaded )
      {
        std::cout << "   ELASTICRODREFFRAMETYPE" << " : " << propsize << " : " << propertyname << std::endl;
        for( size_t i = 0; i < newprop->size(); ++i ) std::cout << "     " << (*newprop)[i] << " "; std::cout << std::endl;
      }
      break;
    }
    case RODBOUNDARYCONDITIONBCLIST:
    {
      Property<RodBoundaryCondition::BCList>* newprop = new Property<RodBoundaryCondition::BCList>(propertyname);
      newprop->resize(propsize);
      BASim::loadProperty(ifs, newprop);
      *prop = newprop;
      if( print_loaded )
      {
        std::cout << "   RODBOUNDARYCONDITIONBCLIST" << " : " << propsize << " : " << propertyname << std::endl;
        for( size_t i = 0; i < newprop->size(); ++i ) for( int j = 0; j < (int) (*newprop)[i].size(); ++j ) std::cout << "     " << (*newprop)[i][j] << " "; std::cout << std::endl;
      }
      break;
    }
    case DOFMAP:
    {
      Property<DOFMap>* newprop = new Property<DOFMap>(propertyname);
      newprop->resize(propsize);
      BASim::loadProperty(ifs, newprop);
      *prop = newprop;
      if( print_loaded )
      {
        std::cout << "   DOFMAP" << " : " << propsize << " : " << propertyname << std::endl;
      }
      break;
    }
    case VEC2D:
    {
      Property<Vec2d>* newprop = new Property<Vec2d>(propertyname);
      newprop->resize(propsize);
      BASim::loadProperty(ifs, newprop);
      *prop = newprop;        
      if( print_loaded )
      {
        std::cout << "   VEC2D" << " : " << propsize << " : " << propertyname << std::endl;
        for( size_t i = 0; i < newprop->size(); ++i ) std::cout << "     " << (*newprop)[i] << " "; std::cout << std::endl;
      }
      break;
    }
    case VEC3D:
    {
      Property<Vec3d>* newprop = new Property<Vec3d>(propertyname);
      newprop->resize(propsize);
      BASim::loadProperty(ifs, newprop);
      *prop = newprop;        
      if( print_loaded )
      {
        std::cout << "   VEC3D" << " : " << propsize << " : " << propertyname << std::endl;
        for( size_t i = 0; i < newprop->size(); ++i ) std::cout << "     " << (*newprop)[i] << " "; std::cout << std::endl;
      }
      break;
    }
    case VECXD:
    {
      Property<VecXd>* newprop = new Property<VecXd>(propertyname);
      newprop->resize(propsize);
      BASim::loadProperty(ifs, newprop);
      *prop = newprop;    
      if( print_loaded )
      {    
        std::cout << "   VECXD" << " : " << propsize << " : " << propertyname << std::endl;
        for( size_t i = 0; i < newprop->size(); ++i ) std::cout << "     " << (*newprop)[i] << " "; std::cout << std::endl;
      }
      break;
    }
    case MAT2D:
    {
      Property<Mat2d>* newprop = new Property<Mat2d>(propertyname);
      newprop->resize(propsize);
      BASim::loadProperty(ifs, newprop);
      *prop = newprop;
      if( print_loaded )
      {
        std::cout << "   MAT2D" << " : " << propsize << " : " << propertyname << std::endl;
        for( size_t i = 0; i < newprop->size(); ++i ) std::cout << "     " << (*newprop)[i] << " "; std::cout << std::endl;
      }
      break;
    }
    case MATXD:
    {
      Property<MatXd>* newprop = new Property<MatXd>(propertyname);
      newprop->resize(propsize);
      BASim::loadProperty(ifs, newprop);
      *prop = newprop;        
      if( print_loaded )
      {
        std::cout << "   MATXD" << " : " << propsize << " : " << propertyname << std::endl;
        for( size_t i = 0; i < newprop->size(); ++i ) std::cout << "     " << (*newprop)[i] << " "; std::cout << std::endl;
      }
      break;
    }
    case VECTORPAIRSCALARSCALAR:
    {
      Property<std::vector< std::pair<Scalar,Scalar> > >* newprop = new Property<std::vector< std::pair<Scalar,Scalar> > >(propertyname);
      newprop->resize(propsize);
      BASim::loadProperty(ifs, newprop);
      *prop = newprop;        
      if( print_loaded )
      {
        std::cout << "   VECTORPAIRSCALARSCALAR" << " : " << propsize << " : " << propertyname << std::endl;
        for( size_t i = 0; i < newprop->size(); ++i )
        {
          for( size_t j = 0; j < (*newprop)[i].size(); ++j ) std::cout << "     " << (*newprop)[i][j].first << "<->" << (*newprop)[i][j].second << " ";
          std::cout << std::endl;
        }
      }
      break;
    }
    case VECTORVEC3D:
    {
      Property<std::vector<Vec3d> >* newprop = new Property<std::vector<Vec3d> >(propertyname);
      newprop->resize(propsize);
      BASim::loadProperty(ifs, newprop);
      *prop = newprop;
      if( print_loaded )
      {
        std::cout << "   VECTORVEC3D" << " : " << propsize << " : " << propertyname << std::endl;
        for( size_t i = 0; i < newprop->size(); ++i )
        {
          for( size_t j = 0; j < (*newprop)[i].size(); ++j ) std::cout << "     " << (*newprop)[i][j] << " ";
          std::cout << std::endl;
        }
      }
      break;
    }
    case PAIRMATXDMATXD:
    {
      Property<std::pair<MatXd,MatXd> >* newprop = new Property<std::pair<MatXd,MatXd> >(propertyname);
      newprop->resize(propsize);
      BASim::loadProperty(ifs, newprop);
      *prop = newprop;        
      if( print_loaded )
      {
        std::cout << "   PAIRMATXDMATXD" << " : " << propsize << " : " << propertyname << std::endl;
        for( size_t i = 0; i < newprop->size(); ++i ) std::cout << "     " << (*newprop)[i].first << "<->" << (*newprop)[i].second << " "; std::cout << std::endl;
      }
      break;
    }
    case PAIRSCALARSCALAR:
    {
      Property<Util::pair<Scalar,Scalar> >* newprop = new Property<Util::pair<Scalar,Scalar> >(propertyname);
      newprop->resize(propsize);
      BASim::loadProperty(ifs, newprop);
      *prop = newprop;
      if( print_loaded )
      {
        std::cout << "   PAIRSCALARSCALAR" << " : " << propsize << " : " << propertyname << std::endl;
        for( size_t i = 0; i < newprop->size(); ++i ) std::cout << "     " << (*newprop)[i].first << "<->" << (*newprop)[i].second << " "; std::cout << std::endl;
      }
      break;
    }
    case PAIRVEC3DVEC3D:
    {
      Property<Util::pair<Vec3d,Vec3d> >* newprop = new Property<Util::pair<Vec3d,Vec3d> >(propertyname);
      newprop->resize(propsize);
      BASim::loadProperty(ifs, newprop);
      *prop = newprop;
      if( print_loaded )
      {
        std::cout << "   PAIRVEC3DVEC3D" << " : " << propsize << " : " << propertyname << std::endl;
        for( size_t i = 0; i < newprop->size(); ++i ) std::cout << "     " << (*newprop)[i].first << "<->" << (*newprop)[i].second << " "; std::cout << std::endl;
      }
      break;
    }
    case VERTEXTOPOLOGY:
    {
      Property<VertexTopology<TopologicalObject> >* newprop = new Property<VertexTopology<TopologicalObject> >(propertyname);
      newprop->resize(propsize);
      BASim::loadProperty(ifs, newprop);
      *prop = newprop;
      if( print_loaded )
      {
        std::cout << "   VERTEXTOPOLOGY" << " : " << propsize << " : " << propertyname << std::endl;
        for( size_t i = 0; i < newprop->size(); ++i )
        {
          for( size_t j = 0; j < (*newprop)[i].size(); ++j ) std::cout << "     " << (*newprop)[i][j].idx() << " ";
          std::cout << std::endl;
        }
      }
      break;
    }
    case EDGETOPOLOGY:
    {
      Property<EdgeTopology<TopologicalObject> >* newprop = new Property<EdgeTopology<TopologicalObject> >(propertyname);
      newprop->resize(propsize);
      BASim::loadProperty(ifs, newprop);
      *prop = newprop;
      if( print_loaded )
      {
        std::cout << "   EDGETOPOLOGY" << " : " << propsize << " : " << propertyname << std::endl;
        for( size_t i = 0; i < newprop->size(); ++i )
        {
          for( size_t j = 0; j < (*newprop)[i].size(); ++j ) std::cout << "     " << (*newprop)[i][j].idx() << " ";
          std::cout << std::endl;
        }
      }
      break;
    }
    case FACETOPOLOGY:
    {
      Property<FaceTopology<TopologicalObject> >* newprop = new Property<FaceTopology<TopologicalObject> >(propertyname);
      newprop->resize(propsize);
      BASim::loadProperty(ifs, newprop);
      *prop = newprop;
      if( print_loaded )
      {
        std::cout << "   FACETOPOLOGY" << " : " << propsize << " : " << propertyname << std::endl;
        for( size_t i = 0; i < newprop->size(); ++i )
        {
          for( size_t j = 0; j < (*newprop)[i].size(); ++j ) std::cout << "     " << (*newprop)[i][j].idx() << " ";
          std::cout << std::endl;
        }        
      }
      break;
    }
    default:
    {
      std::cerr << "\033[31;1mERROR IN TOPOLOGICALOBJECTSERIALIZER:\033[m Attempt to load invalid property type. Exiting." << std::endl;
      exit(1);
      break;
    }
  }    
}
  

}


