// RodRodPenaltyForce.hh
//

#ifndef RODRodPenaltyForce_HH
#define RODRodPenaltyForce_HH

#include "ElasticRod.hh"
#include "RodExternalForce.hh"
#include "../../Collisions/CollisionMeshData.hh"

//#include "RodRodExternalForce.hh"

//#include <ext/hash_map>
#include <map>

namespace BASim {

struct ObjectCollisionInfo
{
    CollisionMeshData* cmData;
    int triangle;
    Vec3d normal;
    double distance;
};

//typedef __gnu_cxx::hash_map<int, ObjectCollisionInfo> VertexObjectMap;
typedef std::map<int, ObjectCollisionInfo> VertexObjectMap;
typedef VertexObjectMap::iterator VertexObjectMapIterator;

//typedef __gnu_cxx::hash_map<int, std::pair<ElasticRod *, int> > EdgeRodMap;
typedef std::map<int, std::pair<ElasticRod *, int> > EdgeRodMap;
typedef EdgeRodMap::iterator EdgeRodMapIterator;

typedef std::map<int, std::pair<ElasticRod *, int> > VertexRodMap;
typedef EdgeRodMap::iterator VertexRodMapIterator;

typedef std::vector<Scalar> RealArray;

typedef std::multimap<int, std::pair<Vec3d, double>> VertexPositionMap;
typedef VertexPositionMap::iterator VertexPositionMapIterator;

// class defined by XX.
class RodVertexConstraint
{
public:
    RodVertexConstraint(  Vec3d target, double stiff, double dis, short type = 0 )
            : m_target( target ),
              m_stiff( stiff ),
              m_restDistance( dis ),
              m_type( type )
            {}

    ~RodVertexConstraint(){}

    enum RodVertexConstraintType {kFix, kRest, kDistance};

    Vec3d  m_target;
    double m_restDistance;
    double m_stiff;
    short  m_type;
};

typedef std::multimap<int, RodVertexConstraint> VertexConstraintMap;
typedef VertexConstraintMap::iterator  VertexConstraintMapIter;


class RodPenaltyForce : public RodExternalForce
{
public:
  RodPenaltyForce();
  ~RodPenaltyForce();

	void clearCollisions();
	
  //void computeEnergy(Real& e)

  virtual void computeForce(const ElasticRod& rod, VecXd& F);

  virtual void computeForceDX(int baseindex, const ElasticRod& rod, Scalar scale, MatrixBase& J);
  virtual void computeForceDV(int baseindex, const ElasticRod& rod, Scalar scale, MatrixBase& J);
  
  void addRodPenaltyForce(int vertex, CollisionMeshData *cmData, int triangle, Collision& collision);
  void addRodPenaltyForce(int edge, ElasticRod *rod, int otherEdge);

  void clearPenaltyForces();

  void computeForceDX(const ElasticRod& rod, Scalar scale, MatrixBase& J);
  void computeForceDV(const ElasticRod& rod, Scalar scale, MatrixBase& J) {}

  VertexConstraintMapIter setVertexPositionPenalty(const ElasticRod* rod,
    int vertex_id, Vec3d& target_position, double stiffness, short type = RodVertexConstraint::kFix);

  // id vertex_id = -1, delete all
  void clearVertexPositionPenalty(int vertex_id = -1);
  
  void addRodClumpingForce(int vertex, ElasticRod *rod, int otherVertex);

  void setClumping(bool flag, Real coeff = 0.0) 
  {
    clumping_enbld = flag;
    clumping_coeff = coeff;
  }
	
  //void clearBulkSprings() {
//    m_bulk_springs.clear();
//  }
  
//  void activateBulkSpring( RodGroupManager::RodRodSpring spr ) {
//    m_bulk_springs.push_back(spr);
//  }

	
protected:
  VertexObjectMap _vertexObjects;
  EdgeRodMap _edgeRods;
  
  VertexRodMap _clumpingVerts;

  Real getClosestPointsVertexTriangle(const Vec3d& v0, const Vec3d& v1,
                                      const Vec3d& v2, const Vec3d& v3,
                                      Real &t1, Real &t2, Real &t3) const;

  Real getClosestPointsEdgeEdge(const Vec3d& e11, const Vec3d& e12,
                                const Vec3d& e21, const Vec3d& e22,
                                Real &s, Real &t) const;
                                
  void localJacobian(MatXd& J, const Scalar stiffness, const Vec3d& normal);

  bool clumping_enbld;
  Real clumping_coeff;
  
  //std::vector <RodGroupManager::RodRodSpring> m_bulk_springs;
                                
  // position based spring penalty force
  //VertexPositionMap m_vertex_position_penalties;
  VertexConstraintMap m_vertex_position_penalties;
};


 
/*struct _RodRodBulkSpringPenalty {
  uint start_v;
  uint nv;
  double rest_length;
  Vec3d opposite_end;
};

typedef struct _RodRodBulkSpringPenalty RodRodBulkSpringPenalty;

  
class RodGroupManager
{
  struct _RodRodSpring {
    uint g;
    uint id1;
    uint id2;
    ElasticRod* rod1;
    ElasticRod* rod2;
    uint start_v;
    uint nv;
    double rest_length;
    double stiff;
    Vec3d opposite_end;
};
  
  typedef struct _RodRodSpring RodRodSpring;
  
  struct RodRodSpringCompare {
    bool operator() (RodRodSpring i,RodRodSpring j) { return (i.rest_length < j.rest_length);}
};
  
  public:
    
    RodGroupManager( std::vector<ElasticRod*> & rods, int ngs, double gap, double stf );
    
    void Setup();
    int ClusterRods( int K );
    void SetupSpring();
    void PrecomputeSpring();

    void ActivateSpring(double dt);

    
  
  
  public:  
    int nr;
    int ns;
    int ng;
    
    double layer_gap;
    double stiffness;
    
    std::vector<ElasticRod*> m_rods;
    std::vector<RodRodSpring> m_springs;
    std::vector<Vec3d> m_rod_locs;
    std::vector<int> m_rod_groups;
    std::vector<int> m_group_nrs;
    
};


//typedef std::vector < std::pair<ElasticRod*, int> > RodVertexVec;
typedef std::vector < std::pair<int, int> > RodVertexVec;
typedef RodVertexVec::iterator RodVertexVecIterator;

class RodGroupSpringForce : public RodRodExternalForce
{
  public:
    RodGroupSpringForce(std::vector<ElasticRod*>& rods);
    ~RodGroupSpringForce();
    
    virtual void computeForce( VecXd& force );
    virtual void computeForceDX( Scalar scale, MatrixBase& J );
    virtual void computeForceDV( Scalar scale, MatrixBase& J ) {}
        
    void checkActivatingCondition();
    int binaryClustering(std::vector<int>& idx, int gid);
    
    void setStiffness(double s) { stiffness = s; }
    double getStiffness() { return stiffness; }
    
    int getNumberOfSprings() { return ns; }
    
    Vec3d& computeCenterPosition( RodVertexVec& rv );
    
  public:
    // spring direction : groupB ==> groupA
    
    int ns; // total number of springs
    int nr;
    int ng;
    std::vector<RodVertexVec> groupA;
    std::vector<RodVertexVec> groupB;
    std::vector<double> restL;
    std::vector<bool> is_activated;
    
    double stiffness;
    
    // each rod info
    std::vector<ElasticRod*> m_rods;
    std::vector<Vec3d> m_rod_locs;     // representative position of each rod for clustering, e.g. each root
    std::vector<int> m_base_id;    // as in global force vector
  
    std::vector<int> m_group_id;
    
    std::vector < std::vector<int> > m_group_rods;  // { {1,3,5}, {2,4}, {7,8,9}, ... }
    
  
};

*/


}

#endif
