/**
 * \file SerializationUtils.hh
 *
 * \author smith@cs.columbia.edu
 * \date 08/24/2010
 */

// TODO: This entire collection of functions is kind of hacky, could make it a lot cleaner
// and maintainable using some basic templates. But it works for now.

// TODO: Could add extra error checks to input/output.

#ifndef SERIALIZATIONUTILS_HH
#define SERIALIZATIONUTILS_HH

#ifdef WETA
#include "../Core/EigenIncludes.hh"
#include "../Core/Property.hh"
#include "../Core/Option.hh"
#include "../Physics/ElasticRods/ElasticRod.hh"
#include "../Physics/ElasticRods/RodStretchingForce.hh"
#include "../Physics/ElasticRods/RodTwistingForceSym.hh"
#include "../Physics/ElasticRods/RodBendingForceSym.hh"
#else
#include "BASim/src/Core/Property.hh"
#include "BASim/src/Core/Option.hh"
//#include "BASim/src/Physics/ElasticRods/ElasticRod.hh"
//#include "BASim/src/Physics/ElasticRods/RodStretchingForce.hh"
//#include "BASim/src/Physics/ElasticRods/RodTwistingForceSym.hh"
//#include "BASim/src/Physics/ElasticRods/RodBendingForceSym.hh"
#endif

#include <queue>
#include <fstream>

namespace BASim 
{

enum ObjectID { UNRECOGNIZEDOBJECTID=-1, RODID=0, TRIMESHID=1 };

enum PropertyTypeID { OBJECTPROPERTYID=1, VERTEXPROPERTYID=2, EDGEPROPERTYID=3, FACEPROPERTYID=4 };

enum PropertyID 
{ 
  INT=5, BOOL=6, SCALAR=7, ELASTICRODREFFRAMETYPE=8, 
  RODBOUNDARYCONDITIONBCLIST=9, ELASTICRODRODFORCES=10, 
  DOFMAP=11, VEC2D=12, VEC3D=13, VECXD=14, MAT2D=15, MATXD=16, 
  VECTORPAIRSCALARSCALAR=17, VECTORVEC3D=18, PAIRMATXDMATXD=19, 
  VERTEXTOPOLOGY=20, PAIRSCALARSCALAR=21, PAIRVEC3DVEC3D=22, 
  EDGETOPOLOGY=23, FACETOPOLOGY=24
};

enum ForceID
{
  RODSTRETCHINGFORCE = 25, RODTWISTINGFORCESYM = 26, RODBENDINGFORCESYM = 27,
  RODSTRETCHINGFORCEVISCOUS = 28, RODBENDINGFORCESYMVISCOUS = 29, RODTWISTINGFORCESYMVISCOUS = 30
};


///////////////////////
// Integer Functions
void serializeVal( std::ofstream& of, const int& val );
void loadVal( std::ifstream& ifs, int& val );

void serializeProperty( std::ofstream& of, Property<int>* prop );
void loadProperty( std::ifstream& ifs, Property<int>* prop );  

///////////////////////
// Vec3d Functions
void printProperty( Property<Vec3d>* prp );
bool propertiesEqual( Property<Vec3d>* prpA, Property<Vec3d>* prpB, Scalar eps = 1.0e-6 );

void serializeVec3d( std::ofstream& of, const Vec3d& val );
void loadVec3d( std::ifstream& ifs, Vec3d& val );

void serializeProperty( std::ofstream& of, Property<Vec3d>* prop );
void loadProperty( std::ifstream& ifs, Property<Vec3d>* prop );

//////////////////////////
// Scalar Functions
void printProperty( Property<Scalar>* prp );
bool propertiesEqual( Property<Scalar>* prpA, Property<Scalar>* prpB, Scalar eps = 1.0e-6 );

void serializeScalar( std::ofstream& of, const Scalar& val );
void loadScalar( std::ifstream& ifs, Scalar& val );

void serializeProperty( std::ofstream& of, Property<Scalar>* prop );
void loadProperty( std::ifstream& ifs, Property<Scalar>* prop );  

//////////////////////////
// String Functions
void serializeString( std::ofstream& of, const std::string& name );
void loadString( std::ifstream& ifs, std::string& name );

//////////////////////////
// PropertyID Functions
void serializePropertyID( std::ofstream& of, const PropertyID& id );

//////////////////////////
// Force Functions
//ForceID getForceID( RodForce* rforce );
//RodForce* getRodForce( const ForceID& fid, ElasticRod* erod );
//
//void serializeElasticRodRodForces( std::ofstream& of, const ElasticRod::RodForces& val );
//void loadElasticRodRodForce( std::ifstream& ifs,  ElasticRod::RodForces& val, ElasticRod* erod );
//
//void serializeProperty( std::ofstream& of, Property<ElasticRod::RodForces>* prop );
//void loadProperty( std::ifstream& ifs, Property<ElasticRod::RodForces>* prop, ElasticRod* erod );
  
//////////////////////////
// Bool Functions
void serializeBool( std::ofstream& of, const bool& val );  
void loadBool( std::ifstream& ifs, bool& val );

void serializeProperty( std::ofstream& of, Property<bool>* prop );
void loadProperty( std::ifstream& ifs, Property<bool>* prop );
  
//////////////////////////
// ElasticRodRefFrame Functions
//void serializeElasticRodRefFrameType( std::ofstream& of, const ElasticRod::RefFrameType& val );  
//void loadElasticRodRefFrameType( std::ifstream& ifs, ElasticRod::RefFrameType& val );
//  
//void serializeProperty( std::ofstream& of, Property<ElasticRod::RefFrameType>* prop );
//void loadProperty( std::ifstream& ifs, Property<ElasticRod::RefFrameType>* prop );
//  
////////////////////////////
//// RodBoundaryConditionBCList Functions
//void serializeRodBoundaryConditionBCList( std::ofstream& of, const RodBoundaryCondition::BCList& val );
//void loadRodBoundaryConditionBCList( std::ifstream& ifs, RodBoundaryCondition::BCList& val );
//
//void serializeProperty( std::ofstream& of, Property<RodBoundaryCondition::BCList>* prop );
//void loadProperty( std::ifstream& ifs, Property<RodBoundaryCondition::BCList>* prop );

//////////////////////////
// DofHandleType Functions
void serializeDofHandleType( std::ofstream& of, const DofHandle::Type& val );
void loadDofHandleType( std::ifstream& ifs, DofHandle::Type& val );
  
//////////////////////////
// DofHandle Functions
void serializeDofHandle( std::ofstream& of, const DofHandle& val );  
void loadDofHandle( std::ifstream& ifs, DofHandle& val );

//////////////////////////
// DOFMap Functions
void serializeDOFMap( std::ofstream& of, DOFMap& val );  
void loadDOFMap( std::ifstream& ifs, DOFMap& val );

void serializeProperty( std::ofstream& of, Property<DOFMap>* prop );  
void loadProperty( std::ifstream& ifs, Property<DOFMap>* prop );

//////////////////////////
// Vec2d Functions
void serializeVec2d( std::ofstream& of, const Vec2d& val );
void loadVec2d( std::ifstream& ifs, Vec2d& val );

void serializeProperty( std::ofstream& of, Property<Vec2d>* prop );  
void loadProperty( std::ifstream& ifs, Property<Vec2d>* prop );
  
//////////////////////////
// VecXd Functions
void serializeVecXd( std::ofstream& of, const VecXd& val );
void loadVecXd( std::ifstream& ifs, VecXd& val );
  
void serializeProperty( std::ofstream& of, Property<VecXd>* prop );
void loadProperty( std::ifstream& ifs, Property<VecXd>* prop );

//////////////////////////
// Mat2d Functions
void serializeMat2d( std::ofstream& of, const Mat2d& val );  
void loadMat2d( std::ifstream& ifs, Mat2d& val );

void serializeProperty( std::ofstream& of, Property<Mat2d>* prop );
void loadProperty( std::ifstream& ifs, Property<Mat2d>* prop );

//////////////////////////
// MatXd Functions
void serializeMatXd( std::ofstream& of, const MatXd& val );  
void loadMatXd( std::ifstream& ifs, MatXd& val );

void serializeProperty( std::ofstream& of, Property<MatXd>* prop );  
void loadProperty( std::ifstream& ifs, Property<MatXd>* prop );

//////////////////////////
// std::pair<Scalar,Scalar> Functions
void serializePairScalarScalar( std::ofstream& of, const std::pair<Scalar,Scalar>& val );  
void loadPairScalarScalar( std::ifstream& ifs, Util::pair<Scalar,Scalar>& val );
  
//////////////////////////
// std::vector< std::pair<Scalar,Scalar> > Functions
void serializeVectorPairScalar( std::ofstream& of, const std::vector< std::pair<Scalar,Scalar> >& val );
void loadVectorPairScalar( std::ifstream& ifs, std::vector< std::pair<Scalar,Scalar> >& val );
  
void serializeProperty( std::ofstream& of, Property<std::vector< std::pair<Scalar,Scalar> > >* prop );  
void loadProperty( std::ifstream& ifs, Property<std::vector< std::pair<Scalar,Scalar> > >* prop );

//////////////////////////
// std::vector<Vec3d> Functions
void serializeVectorVec3d( std::ofstream& of, const std::vector<Vec3d>& val );
void loadVectorVec3d( std::ifstream& ifs, std::vector<Vec3d>& val );  
  
void serializeProperty( std::ofstream& of, Property<std::vector<Vec3d> >* prop );  
void loadProperty( std::ifstream& ifs, Property<std::vector<Vec3d> >* prop );

//////////////////////////
// std::vector<bool> Functions
void serializeVectorBool( std::ofstream& of, const std::vector<bool>& val );
void loadVectorBool( std::ifstream& ifs, std::vector<bool>& val );  

//////////////////////////
// std::vector<unsigned int> Functions
void serializeVectorUint( std::ofstream& of, const std::vector<unsigned int>& val );
void loadVectorUint( std::ifstream& ifs, std::vector<unsigned int>& val );  

//////////////////////////
// std::pair<MatXd,MatXd> Functions
void serializePairMatXdMatXd( std::ofstream& of, const std::pair<MatXd,MatXd>& val );  
void loadPairMatXdMatXd( std::ifstream& ifs, std::pair<MatXd,MatXd>& val );

void serializeProperty( std::ofstream& of, Property<std::pair<MatXd,MatXd> >* prop );  
void loadProperty( std::ifstream& ifs, Property<std::pair<MatXd,MatXd> >* prop );
  
//////////////////////////
// Util::pair<Scalar,Scalar> Functions
void serializePairScalarScalar( std::ofstream& of, const Util::pair<Scalar,Scalar>& val );  
void loadPairScalarScalar( std::ifstream& ifs, Util::pair<Scalar,Scalar>& val );

void serializeProperty( std::ofstream& of, Property<Util::pair<Scalar,Scalar> >* prop );
void loadProperty( std::ifstream& ifs, Property<Util::pair<Scalar,Scalar> >* prop );

//////////////////////////
// Util::pair<Vec3d,Vec3d> Functions
void serializePairVec3dVec3d( std::ofstream& of, const Util::pair<Vec3d,Vec3d>& val );
void loadPairVec3dVec3d( std::ifstream& ifs, Util::pair<Vec3d,Vec3d>& val );

void serializeProperty( std::ofstream& of, Property<Util::pair<Vec3d,Vec3d> >* prop );
void loadProperty( std::ifstream& ifs, Property<Util::pair<Vec3d,Vec3d> >* prop );

//////////////////////////
// VertexTopology<TopologicalObject> Functions
void serializeVeretexTopology( std::ofstream& of, const VertexTopology<TopologicalObject>& val );  
void loadVeretexTopology( std::ifstream& ifs, VertexTopology<TopologicalObject>& val );
  
void serializeProperty( std::ofstream& of, Property<VertexTopology<TopologicalObject> >* prop );  
void loadProperty( std::ifstream& ifs, Property<VertexTopology<TopologicalObject> >* prop );

//////////////////////////
// EdgeTopology<TopologicalObject> Functions
void serializeEdgeTopology( std::ofstream& of, const EdgeTopology<TopologicalObject>& val );
void loadEdgeTopology( std::ifstream& ifs, EdgeTopology<TopologicalObject>& val );

void serializeProperty( std::ofstream& of, Property<EdgeTopology<TopologicalObject> >* prop );
void loadProperty( std::ifstream& ifs, Property<EdgeTopology<TopologicalObject> >* prop );

//////////////////////////
// FaceTopology<TopologicalObject> Functions
void serializeFaceTopology( std::ofstream& of, const FaceTopology<TopologicalObject>& val );  
void loadFaceTopology( std::ifstream& ifs, FaceTopology<TopologicalObject>& val );

void serializeProperty( std::ofstream& of, Property<FaceTopology<TopologicalObject> >* prop );  
void loadProperty( std::ifstream& ifs, Property<FaceTopology<TopologicalObject> >* prop );

//////////////////////////
// std::queue<double> Functions
void serializeDoubleQueue( std::ofstream& of, const std::queue<double>& val );  
void loadDoubleQueue( std::ifstream& ifs, std::queue<double>& val );

//////////////////////////
// std::vector< std::vector<int> > Functions
void serializeVectorVectorInt( std::ofstream& of, const std::vector< std::vector<int> >& val );
void loadVectorVectorInt( std::ifstream& ifs, std::vector< std::vector<int> > & val );


//////////////////////////
// Option Functions
void serializeOptionType( std::ofstream& of, const Option::Type& val );
void loadOptionType( std::ifstream& ifs, Option::Type& val );

void serializeOption( std::ofstream& of, const Option& val );
void loadOption( std::ifstream& ifs, Option& val );

void serializeMapStringOption( std::ofstream& of, const std::map<std::string,Option>& val );
void loadMapStringOption( std::ifstream& ifs, std::map<std::string,Option>& val );


//////////////////////////
// General Functions
void serializeProperty( std::ofstream& of, PropertyBase* prop );
void loadProperty( std::ifstream& ifs, PropertyBase** prop, ObjectBase* topobj, int propsize );

}

#endif
