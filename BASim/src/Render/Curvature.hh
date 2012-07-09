#ifndef CURVATURE_HH
#define CURVATURE_HH


//Adapted from Etienne/Yotam's code for computing mesh curvature data

#include <vector>
#include "BASim/src/Core/TopologicalObject/TopologicalObject.hh"
#include "BASim/src/Core/EigenIncludes.hh"
#include "BASim/src/Core/TopologicalObject/TopObjProperty.hh"

namespace BASim {

struct Frame
{
    Frame(const Eigen::Vector3d &normal);
    Eigen::Vector3d normal;
    Eigen::Vector3d u;
    Eigen::Vector3d v;
};

struct ShapeOperator
{
    ShapeOperator() : frame(Vec3d(1,0,0)) {S.setZero();}
    ShapeOperator(Frame frame, const Eigen::Matrix2d &S) : frame(frame), S(S) {}
    ShapeOperator(Frame frame) : frame(frame) {S.setZero();}
    Frame frame;
    Eigen::Matrix2d S;
    bool notConverged;
};

struct CurvatureInfo
{
    CurvatureInfo();
    double principalCurvature[2];
    Eigen::Vector3d curvatureDir[2];
    bool notConverged;
};

class MeshCurvature
{
public:
    MeshCurvature(TopologicalObject& obj, const VertexProperty<Vec3d>& vertices);
    
    void renderCurvatureDirs();
    
    double gaussianCurvature(const VertexHandle& vh);
    double meanCurvature(const VertexHandle& vh);
    double curvatureSpread(const VertexHandle& vh);
    double vertexArea(const VertexHandle& vh);

    void totalSquaredGaussianCurvature(std::pair<double, double> &values, double cutoff);

private:
    VertexProperty<CurvatureInfo> curvature_;

    double curCutoff_;
    double aboveSqGaussCurvature_;
    double belowSqGaussCurvature_;

    
    void computeCurvatures();
    void computeCurvature(const ShapeOperator &shapeOperator, CurvatureInfo &curvature);
    void computeShapeOperators(VertexProperty<ShapeOperator> &operators);
    

    void recomputeSquaredGaussianCurvature(double cutoff);
    ShapeOperator computeShapeOperator(const VertexHandle& vh);
    
    ShapeOperator estimateShapeOperator(const VertexHandle& vh, Eigen::Vector3d &normal);
    ShapeOperator estimateShapeOperatorNew(const VertexHandle& vh, Eigen::Vector3d &normal);

    Eigen::Matrix2d transportShapeOperator(const ShapeOperator &source, const ShapeOperator &dest);
    void smoothShapeOperators(const VertexProperty<ShapeOperator> &oldOperators, VertexProperty<ShapeOperator> &newOperators);

    TopologicalObject& topology_;
    const VertexProperty<Vec3d>& vertices_;
};

}

#endif // CURVATURE_HH
