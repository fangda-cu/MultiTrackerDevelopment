#include "BASim/src/Render/Curvature.hh"

#include "BASim/src/Render/OpenGLHeaders.hh"

#include <set>
#include <Eigen/Cholesky>
#include <Eigen/Dense>
#include <fstream>

using namespace std;
using namespace Eigen;

namespace BASim {

CurvatureInfo::CurvatureInfo()
{
    principalCurvature[0] = 0;
    principalCurvature[1] = 0;
    curvatureDir[0].setZero();
    curvatureDir[1].setZero();
}

MeshCurvature::MeshCurvature(TopologicalObject& topo, const VertexProperty<Vec3d> &v): topology_(topo), vertices_(v), curvature_(&topo)
{
    computeCurvatures();
}

void MeshCurvature::computeCurvatures()
{

    VertexProperty<ShapeOperator> shapeOperators(&topology_);
    computeShapeOperators(shapeOperators);
    VertexProperty<ShapeOperator> smoothedOperators(&topology_);
    
    //Smooth twice?
    smoothShapeOperators(shapeOperators, smoothedOperators);
    shapeOperators = smoothedOperators;
    smoothShapeOperators(shapeOperators, smoothedOperators);
    
    for(VertexIterator vit = topology_.vertices_begin(); vit != topology_.vertices_end(); ++vit)
    {
        VertexHandle vh = *vit;
        CurvatureInfo c;
        computeCurvature(smoothedOperators[vh], c);
        curvature_[vh] = c;
        curvature_[vh].notConverged = smoothedOperators[vh].notConverged;
    }
    
}

Frame::Frame(const Eigen::Vector3d &normal) : normal(normal)
{
    Vector3d seed(0,0,0);
    if(fabs(normal[0]) < 0.5)
        seed[0] = 1.0;
    else
        seed[1] = 1.0;
    u = normal.cross(seed);
    u.normalize();
    v = normal.cross(u);
}

void MeshCurvature::computeCurvature(const ShapeOperator &shapeOperator, CurvatureInfo &curvature)
{
    SelfAdjointEigenSolver<Matrix2d> eigensolver(shapeOperator.S);
    for(int i=0; i<2; i++)
    {
        curvature.principalCurvature[i] = eigensolver.eigenvalues()[i];
        curvature.curvatureDir[i] = eigensolver.eigenvectors().coeffRef(0,i)*shapeOperator.frame.u + eigensolver.eigenvectors().coeffRef(1,i)*shapeOperator.frame.v;
        curvature.curvatureDir[i].normalize();
    }
    if(fabs(curvature.principalCurvature[1]) < fabs(curvature.principalCurvature[0]))
    {
        swap(curvature.principalCurvature[0],curvature.principalCurvature[1]);
        swap(curvature.curvatureDir[0],curvature.curvatureDir[1]);
    }
}

void MeshCurvature::renderCurvatureDirs()
{


    glPolygonOffset(1.0, 1.0);

    glLineWidth(4.0);
    glBegin(GL_LINES);
    for(VertexIterator vit = topology_.vertices_begin(); vit != topology_.vertices_end(); ++vit)
    {
        {
            VertexHandle vh = *vit;
            Vec3d pt = vertices_[vh];
            /*if(curvature_[vh].notConverged) 
                continue;*/

            Vector3d ept(pt[0],pt[1],pt[2]);
            
            VertexEdgeIterator veit = topology_.ve_iter(vh);
            EdgeHandle eh = *veit;
            double radius = (vertices_[topology_.fromVertex(eh)] - vertices_[topology_.toVertex(eh)]).norm();

            double strength = max(fabs(curvature_[vh].principalCurvature[0]), fabs(curvature_[vh].principalCurvature[1])) / 0.2;
            //double strength = 0.5*(curvature_[vh].principalCurvature[0]+ curvature_[vh].principalCurvature[1]) / 0.06;

            strength = max(strength, 0.0);
            strength = min(strength, 1.0);
            glColor3f(0.0, (float)strength, 0.0);
            
            //std::cout << curvature_[vh].principalCurvature[0] << " ";

            Vector3d pt1 = ept - 0.5*radius*curvature_[vh].curvatureDir[0];
            Vector3d pt2 = ept + 0.5*radius*curvature_[vh].curvatureDir[0];

            glVertex3d(pt1[0], pt1[1], pt1[2]);
            glVertex3d(pt2[0], pt2[1], pt2[2]);
        }
    }
    glEnd();
    
    glPointSize(10);
    glBegin(GL_POINTS);
    glColor3f(1.0f,0.0f,1.0f);
    for(VertexIterator vit = topology_.vertices_begin(); vit != topology_.vertices_end(); ++vit)
    {
        VertexHandle vh = *vit;
        Vec3d pt = vertices_[vh];
        if(curvature_[vh].notConverged)  {
            glVertex3d(pt[0], pt[1], pt[2]);
        }
    }
    glEnd();
    glPointSize(1);
    glPolygonOffset(0.0, 0.0);

}


double MeshCurvature::gaussianCurvature(const VertexHandle& vh)
{
    assert(vh.isValid());
    return curvature_[vh].principalCurvature[0]*curvature_[vh].principalCurvature[1];
}


double MeshCurvature::meanCurvature(const VertexHandle& vh)
{
    assert(vh.isValid());
    return 0.5*(curvature_[vh].principalCurvature[0]+curvature_[vh].principalCurvature[1]);
}

double MeshCurvature::vertexArea(const VertexHandle& vh) {
    assert(vh.isValid());
    double totalArea = 0;
    for(VertexFaceIterator vfit = topology_.vf_iter(vh); vfit; ++vfit) {
        FaceHandle fh = *vfit;
        FaceVertexIterator fvit = topology_.fv_iter(fh);
        VertexHandle v0 = *fvit; ++fvit;
        VertexHandle v1 = *fvit; ++fvit;
        VertexHandle v2 = *fvit;
        Vec3d pos0 = vertices_[v0];
        Vec3d pos1 = vertices_[v1];
        Vec3d pos2 = vertices_[v2];
        Vec3d areaVector = 0.5*(pos1 - pos0).cross(pos2-pos0);
        totalArea += areaVector.norm() / 3.0;
    }
    return totalArea;
}


double MeshCurvature::curvatureSpread(const VertexHandle& vh)
{
    assert(vh.isValid());
    return fabs(fabs(curvature_[vh].principalCurvature[0])-fabs(curvature_[vh].principalCurvature[1]));
}

void MeshCurvature::computeShapeOperators(VertexProperty<ShapeOperator> &operators)
{

    for(VertexIterator vit = topology_.vertices_begin(); vit != topology_.vertices_end(); ++vit)
    {
        VertexHandle vh = *vit;
        operators[vh] = computeShapeOperator(vh);
    }
}

ShapeOperator MeshCurvature::computeShapeOperator(const VertexHandle& vh)
{
    //Estimate initial vertex normal for vh, as plain old average of incident face normals
    Vector3d oldnormal(0,0,0); 
    int count = 0;
    for(VertexFaceIterator vfit = topology_.vf_iter(vh); vfit; ++vfit) {
        FaceHandle fh = *vfit;
        FaceVertexIterator fvit = topology_.fv_iter(fh);
        VertexHandle v0 = *fvit; ++fvit;
        VertexHandle v1 = *fvit; ++fvit;
        VertexHandle v2 = *fvit;
        Vec3d pos0 = vertices_[v0];
        Vec3d pos1 = vertices_[v1];
        Vec3d pos2 = vertices_[v2];
        Vec3d faceNormal = (pos1 - pos0).cross(pos2-pos0);
        faceNormal.normalize();
        oldnormal += faceNormal;
        ++count;
    }
    oldnormal/=count;
    oldnormal.normalize();

    Vector3d normal = oldnormal;
    ShapeOperator result = estimateShapeOperator(vh, normal);
    int iters = 1;
    //while( (oldnormal-normal).norm() > 1e-4  && iters < 40)
    while( acos(oldnormal.dot(normal)) > 0.1*3.141592653589793/180 && iters < 10) //while it differs by more than 1/10 of a degree
    {
        oldnormal = normal;
        result = estimateShapeOperator(vh, normal);
        iters++;
    }
    result.notConverged = (iters == 40);
    
    return result;
}

ShapeOperator MeshCurvature::estimateShapeOperator(const VertexHandle& vh, Eigen::Vector3d &normal)
{
    Frame frame(normal);

    set<int> twoneighbors;
    twoneighbors.insert(vh.idx());
    for(int i=0; i<2; i++)
    {
        set<int> toadd;
        for(set<int>::iterator it = twoneighbors.begin(); it != twoneighbors.end(); ++it)
        {
            VertexHandle curVert(*it);
            for(VertexEdgeIterator veit = topology_.ve_iter(curVert); veit; ++veit) 
            {
                EdgeHandle eh = *veit;
                VertexHandle fv = topology_.fromVertex(eh);
                if(fv == curVert)
                    toadd.insert(topology_.toVertex(eh).idx());
                else
                    toadd.insert(fv.idx());
            }
        }
        for(set<int>::iterator it = toadd.begin(); it != toadd.end(); ++it)
            twoneighbors.insert(*it);
    }

    int numneighbors = twoneighbors.size();
    if(numneighbors < 6)
        return ShapeOperator(frame);

    //a u^2 + b v^2 + c uv + d u + e v + f

    MatrixXd Mat(numneighbors, 6);
    VectorXd rhs(numneighbors);
    Vec3d centpt = vertices_[vh];
    int row=0;
    for(set<int>::iterator it = twoneighbors.begin(); it != twoneighbors.end(); ++it)
    {
        VertexHandle curVert(*it);
        Vec3d adjpt = vertices_[curVert];
        Vector3d diff;
        for(int i=0; i<3; i++)
            diff[i] = adjpt[i]-centpt[i];
        double uval = diff.dot(frame.u);
        double vval = diff.dot(frame.v);
        rhs[row] = diff.dot(normal);
        Mat.coeffRef(row,0) = uval*uval;
        Mat.coeffRef(row,1) = vval*vval;
        Mat.coeffRef(row,2) = uval*vval;
        Mat.coeffRef(row,3) = uval;
        Mat.coeffRef(row,4) = vval;
        Mat.coeffRef(row,5) = 1;
        row++;
    }
    assert(row == numneighbors);
    /*VectorXd nrhs = Mat.transpose()*rhs;
    MatrixXd MTM = Mat.transpose()*Mat;
    VectorXd coeffs = MTM.ldlt().solve(nrhs);*/

    JacobiSVD<MatrixXd> svd(Mat, ComputeThinU | ComputeThinV);
    VectorXd coeffs = svd.solve(rhs);

    Vector3d newnormal = (normal - coeffs[4]*frame.v - coeffs[3]*frame.u)/sqrt(1+coeffs[3]*coeffs[3]+coeffs[4]*coeffs[4]);

    double E = 1 + coeffs[3]*coeffs[3];
    double F = coeffs[3]*coeffs[4];
    double G = 1 + coeffs[4]*coeffs[4];

    double fac = normal.dot(newnormal);
    double L = 2*coeffs[0]*fac;
    double M = coeffs[2]*fac;
    double N = 2*coeffs[1]*fac;

    double det = E*G-F*F;
    if(fabs(det) < 1e-6)
        return ShapeOperator(frame);

    Matrix2d shapeOperator;
    shapeOperator << F*M-G*L, F*L-E*M, F*N-G*M, F*M-E*N;
    shapeOperator /= det;

    Matrix2d g;
    g << 1 + coeffs[4]*coeffs[4], -coeffs[3]*coeffs[4], -coeffs[3]*coeffs[4], 1 + coeffs[3]*coeffs[3];
    g /= det;

    shapeOperator = g*shapeOperator;

    assert(fabs(shapeOperator.coeff(0,1)-shapeOperator.coeff(1,0)) < 1e-6);

    normal = newnormal;

    return ShapeOperator(frame, shapeOperator);
}

//This version follows more closely the local height function formulas in Jiao & Zha 2008 / Wang et al. 2009
//Note Wang et al. incorrectly state l = (1 +fu^2+fv^2) instead of l = sqrt(1 +fu^2+fv^2)
ShapeOperator MeshCurvature::estimateShapeOperatorNew(const VertexHandle& vh, Eigen::Vector3d &normal)
{
    Frame frame(normal);

    set<int> twoneighbors;
    twoneighbors.insert(vh.idx());
    for(int i=0; i<2; i++)
    {
        set<int> toadd;
        for(set<int>::iterator it = twoneighbors.begin(); it != twoneighbors.end(); ++it)
        {
            VertexHandle curVert(*it);
            for(VertexEdgeIterator veit = topology_.ve_iter(curVert); veit; ++veit) 
            {
                EdgeHandle eh = *veit;
                VertexHandle fv = topology_.fromVertex(eh);
                if(fv == curVert)
                    toadd.insert(topology_.toVertex(eh).idx());
                else
                    toadd.insert(fv.idx());
            }
        }
        for(set<int>::iterator it = toadd.begin(); it != toadd.end(); ++it)
            twoneighbors.insert(*it);
    }

    int numneighbors = twoneighbors.size();
    if(numneighbors < 6)
        return ShapeOperator(frame);

    //a u^2 + b v^2 + c uv + d u + e v + f

    MatrixXd Mat(numneighbors, 6);
    VectorXd rhs(numneighbors);
    Vec3d centpt = vertices_[vh];
    int row=0;
    for(set<int>::iterator it = twoneighbors.begin(); it != twoneighbors.end(); ++it)
    {
        VertexHandle curVert(*it);
        Vec3d adjpt = vertices_[curVert];
        Vector3d diff;
        for(int i=0; i<3; i++)
            diff[i] = adjpt[i]-centpt[i];
        double uval = diff.dot(frame.u);
        double vval = diff.dot(frame.v);
        rhs[row] = diff.dot(normal);
        Mat.coeffRef(row,0) = uval*uval;
        Mat.coeffRef(row,1) = vval*vval;
        Mat.coeffRef(row,2) = uval*vval;
        Mat.coeffRef(row,3) = uval;
        Mat.coeffRef(row,4) = vval;
        Mat.coeffRef(row,5) = 1;
        row++;
    }
    assert(row == numneighbors);
    VectorXd nrhs = Mat.transpose()*rhs;
    MatrixXd MTM = Mat.transpose()*Mat;
    VectorXd coeffs = MTM.ldlt().solve(nrhs);

    double areaElt = sqrt(1+coeffs[3]*coeffs[3]+coeffs[4]*coeffs[4]);
    Vector3d newnormal = (normal - coeffs[3]*frame.u - coeffs[4]*frame.v) / areaElt;

    double E = 1 + coeffs[3]*coeffs[3];
    double F = coeffs[3]*coeffs[4];
    double G = 1 + coeffs[4]*coeffs[4];

    //double fac = normal.dot(newnormal);
    double L = 2*coeffs[0]/areaElt;
    double M = coeffs[2]/areaElt;
    double N = 2*coeffs[1]/areaElt;

    double det = E*G-F*F;
    if(fabs(det) < 1e-6)
        return ShapeOperator(frame);

    Matrix2d H;
    H << L , M, M, N;
    
    Matrix3d U;
    U << frame.u[0], frame.v[0], frame.normal[0], 
        frame.u[1], frame.v[1], frame.normal[1], 
        frame.u[2], frame.v[2], frame.normal[2]; 

    Matrix<Scalar, 3, 2> J;
    J << 1, 0,
        0, 1, 
        coeffs[3], coeffs[4];
    J = U*J;
    

    Matrix2d shapeOperator;
    
    //shapeOperator << F*M-G*L, F*L-E*M, F*N-G*M, F*M-E*N;
    //shapeOperator /= det;

    /*Matrix2d g;
    g << 1 + coeffs[4]*coeffs[4], -coeffs[3]*coeffs[4], -coeffs[3]*coeffs[4], 1 + coeffs[3]*coeffs[3];
    g /= det;

    shapeOperator = g*shapeOperator;
*/
    shapeOperator = 1/areaElt * (J.transpose()*J).inverse() * H;

    assert(fabs(shapeOperator.coeff(0,1)-shapeOperator.coeff(1,0)) < 1e-6);

    normal = newnormal;

    return ShapeOperator(frame, shapeOperator);
}


Matrix2d MeshCurvature::transportShapeOperator(const ShapeOperator &source, const ShapeOperator &dest)
{
    MatrixXd sourceM(3,2);
    sourceM << source.frame.u, source.frame.v;

    MatrixXd destM(3,2);
    destM << dest.frame.u, dest.frame.v;

    Matrix2d M = sourceM.transpose()*sourceM;
    Matrix2d T = M.inverse() * sourceM.transpose()*destM;
    return T.transpose()*source.S*T;
}

void MeshCurvature::smoothShapeOperators(const VertexProperty<ShapeOperator> &oldOperators, VertexProperty<ShapeOperator> &newOperators)
{

    for(VertexIterator vit = topology_.vertices_begin(); vit != topology_.vertices_end(); ++vit)
    {
        VertexHandle curVert = *vit;
        int denom = 1;
        Matrix2d S = oldOperators[curVert].S;
        for(VertexEdgeIterator veit = topology_.ve_iter(curVert); veit; ++veit) 
        {
            EdgeHandle eh = *veit;
            VertexHandle nbrVert = topology_.fromVertex(eh);
            if(nbrVert == curVert)
                nbrVert = topology_.toVertex(eh);
            
            S += transportShapeOperator(oldOperators[nbrVert], oldOperators[nbrVert]);
            denom++;
        }
        S /= denom;
        newOperators[curVert] = ShapeOperator(oldOperators[curVert].frame, S);
    }
}

void MeshCurvature::recomputeSquaredGaussianCurvature(double cutoff)
{
    double totabove = 0;
    double totbelow = 0;
    for(VertexIterator vit = topology_.vertices_begin(); vit != topology_.vertices_end(); ++vit)
    {
        VertexHandle vh = *vit;
        double cur = gaussianCurvature(vh);
        double area = vertexArea(vh);
        if(fabs(cur) > cutoff)
            totabove += area*cur*cur;
        else
            totbelow += area*cur*cur;
    }
    curCutoff_ = cutoff;
    belowSqGaussCurvature_ = totbelow;
    aboveSqGaussCurvature_ = totabove;
}

void MeshCurvature::totalSquaredGaussianCurvature(pair<double, double> &values, double cutoff)
{
    if(cutoff != curCutoff_)
    {
        recomputeSquaredGaussianCurvature(cutoff);
    }
    values.first = belowSqGaussCurvature_;
    values.second = aboveSqGaussCurvature_;
}

}