// LevelSet.cc
//

#include "LevelSet.hh"
#include "../Util/TextLog.hh"
//#include <omp.h>
//#include <GL/gl.h>
//#include <OpenGL/OPENGL_COLOR_RAMP.h>
//#include <OpenGL/OPENGL_COLOR.h>

using namespace std;

namespace BASim {

LevelSet::LevelSet()
    : _initialized(false)
{
    m_realRequestPositions.clear();
    m_transformedRequestPositions.clear();
    m_grad.clear();
    m_gradPosition.clear();
}

LevelSet::~LevelSet()
{
}
/*
void LevelSet::buildLevelSet(const std::vector<bridson::Vec3ui> &triangles,
                             const std::vector<unsigned int> &triIndices,
                             const std::vector<bridson::Vec3f>  &x,
                             const std::vector<bridson::Vec3f>  &v,
                             const bridson::Vec3f &origin, Real length[3],
                             Real dx, int nx, int ny, int nz, int nbrTriangles)
*/
void LevelSet::buildLevelSet(const Vec3Indices &triangles,
                             const Indices &triIndices,
                             const std::vector<bridson::Vec3f>  &x,
                             const std::vector<bridson::Vec3f>  &v,
                             const bridson::Vec3f &origin, Real length[3],
                             Real dx, int nx, int ny, int nz, int nbrTriangles,
                             Eigen::Matrix4f& transformMatrix)
{
    for (unsigned int i=0; i<3; ++i)
        _origin[i] = origin[i];

    _dx = dx;

    //TraceStream(g_log, "LevelSet::buildLevelSet") << "Level set dimensions: (" << nx << "," << ny << "," << nz << ") " << _dx << "\n";
    buildLevelSet(triangles, triIndices, x, v, _origin, _dx, nx, ny, nz, _phi, _phiVel);
    _initialized = true;
    
    // Maya matrices are flipped from Eigen matrices
    // TODO: Transpose these back in Maya land
    m_transformMatrixAtCreation = transformMatrix.transpose();
}

void LevelSet::setTransformationMatrix( Eigen::Matrix4f& i_matrix )
{
    // Maya matrices are flipped from Eigen matrices
    // TODO: Transpose these back in Maya land
    m_currentTransformMatrix = i_matrix.transpose();
}

Real LevelSet::getLevelSetValue(Vec3<Real> x)
{    
    Eigen::Vector4f samplePoint( x[ 0 ], x[ 1 ], x[ 2 ], 1 );        
    
    //m_realRequestPositions.push_back( samplePoint );
    
    Eigen::Matrix4f concatMatrix = ( m_transformMatrixAtCreation * m_currentTransformMatrix.inverse() );
    samplePoint = concatMatrix * samplePoint;
    
    //m_transformedRequestPositions.push_back( samplePoint );;
    
    Real fi = ((double)samplePoint[0] - _origin[0]) / _dx;
    Real fj = ((double)samplePoint[1] - _origin[1]) / _dx;
    Real fk = ((double)samplePoint[2] - _origin[2]) / _dx;

    //TraceStream(g_log, "LevelSet::getLevelSetValue") << "this = " << this << " Array " << &_phi << " dim " << _phi.ni << ", " << _phi.nj << ", " << _phi.nk << "\n";    
   
    int i, j, k;
    bridson::get_barycentric(fi, i, fi, 0, _phi.ni);
    bridson::get_barycentric(fj, j, fj, 0, _phi.nj);
    bridson::get_barycentric(fk, k, fk, 0, _phi.nk);

    //TraceStream(g_log, "LevelSet::getLevelSetValue") << "Array dim " << _phi.ni << ", " << _phi.nj << ", " << _phi.nk  << "  ijk = " << i << ", " << j << ", " << k << "  fijk = " << fi << ", " << fj << ", " << fk << "\n";    

    assert(_initialized);

    Real dist = bridson::trilerp(_phi(i    , j    , k    ),
                                 _phi(i + 1, j    , k    ),
                                 _phi(i    , j + 1, k    ),
                                 _phi(i + 1, j + 1, k    ),
                                 _phi(i    , j    , k + 1),
                                 _phi(i + 1, j    , k + 1),
                                 _phi(i    , j + 1, k + 1),
                                 _phi(i + 1, j + 1, k + 1),
                                 fi, fj, fk);

    return dist;
}

Real LevelSet::getLevelSetValueVelocity(Vec3<Real> &x, Vec3<Real> &v)
{
    Eigen::Vector4f samplePoint( x[ 0 ], x[ 1 ], x[ 2 ], 1 );        
    
    //m_realRequestPositions.push_back( samplePoint );
    
    Eigen::Matrix4f concatMatrix = ( m_transformMatrixAtCreation * m_currentTransformMatrix.inverse() );
    samplePoint = concatMatrix * samplePoint;
    
    Real fi = ((double)samplePoint[0] - _origin[0]) / _dx;
    Real fj = ((double)samplePoint[1] - _origin[1]) / _dx;
    Real fk = ((double)samplePoint[2] - _origin[2]) / _dx;
   
    int i, j, k;
    bridson::get_barycentric(fi, i, fi, 0, _phi.ni);
    bridson::get_barycentric(fj, j, fj, 0, _phi.nj);
    bridson::get_barycentric(fk, k, fk, 0, _phi.nk);
   
    Real dist = bridson::trilerp(_phi(i    , j    , k    ),
                                 _phi(i + 1, j    , k    ),
                                 _phi(i    , j + 1, k    ),
                                 _phi(i + 1, j + 1, k    ),
                                 _phi(i    , j    , k + 1),
                                 _phi(i + 1, j    , k + 1),
                                 _phi(i    , j + 1, k + 1),
                                 _phi(i + 1, j + 1, k + 1),
                                 fi, fj, fk);

    bridson::Vec3f vel = bridson::trilerp(_phiVel(i    , j    , k    ),
                                          _phiVel(i + 1, j    , k    ),
                                          _phiVel(i    , j + 1, k    ),
                                          _phiVel(i + 1, j + 1, k    ),
                                          _phiVel(i    , j    , k + 1),
                                          _phiVel(i + 1, j    , k + 1),
                                          _phiVel(i    , j + 1, k + 1),
                                          _phiVel(i + 1, j + 1, k + 1),
                                          fi, fj, fk);

    for (int i=0; i<3; ++i)
        v[i] = vel[i];
    
    return dist;
}

void LevelSet::getGradient(Vec3<Real> &x, Vec3<Real> &grad)
{
    Eigen::Vector4f samplePoint( x[ 0 ], x[ 1 ], x[ 2 ], 1 );        

    //m_gradPosition.push_back( samplePoint );

    Eigen::Matrix4f concatMatrix = ( m_transformMatrixAtCreation * m_currentTransformMatrix.inverse() );
    samplePoint = concatMatrix * samplePoint;
    
    
    float i = ((double)samplePoint[0] - _origin[0]) / _dx;
    float j = ((double)samplePoint[1] - _origin[1]) / _dx;
    float k = ((double)samplePoint[2] - _origin[2]) / _dx;

    int p, q, r;
    bridson::get_barycentric(i, p, i, 1, _phi.ni);
    bridson::get_barycentric(j, q, j, 1, _phi.nj);
    bridson::get_barycentric(k, r, k, 1, _phi.nk);

    Real phiP1 = bridson::trilerp(_phi(p    , q    , r    ),
                                  _phi(p + 1, q    , r    ),
                                  _phi(p    , q + 1, r    ),
                                  _phi(p + 1, q + 1, r    ),
                                  _phi(p    , q    , r + 1),
                                  _phi(p + 1, q    , r + 1),
                                  _phi(p    , q + 1, r + 1),
                                  _phi(p + 1, q + 1, r + 1),
                                  1.0f, j, k);

    Real phiP0 = bridson::trilerp(_phi(p    , q    , r    ),
                                  _phi(p + 1, q    , r    ),
                                  _phi(p    , q + 1, r    ),
                                  _phi(p + 1, q + 1, r    ),
                                  _phi(p    , q    , r + 1),
                                  _phi(p + 1, q    , r + 1),
                                  _phi(p    , q + 1, r + 1),
                                  _phi(p + 1, q + 1, r + 1),
                                  0.0f, j, k);

    Real phiQ1 = bridson::trilerp(_phi(p    , q    , r    ),
                                  _phi(p + 1, q    , r    ),
                                  _phi(p    , q + 1, r    ),
                                  _phi(p + 1, q + 1, r    ),
                                  _phi(p    , q    , r + 1),
                                  _phi(p + 1, q    , r + 1),
                                  _phi(p    , q + 1, r + 1),
                                  _phi(p + 1, q + 1, r + 1),
                                  i, 1.0f, k);

    Real phiQ0 = bridson::trilerp(_phi(p    , q    , r    ),
                                  _phi(p + 1, q    , r    ),
                                  _phi(p    , q + 1, r    ),
                                  _phi(p + 1, q + 1, r    ),
                                  _phi(p    , q    , r + 1),
                                  _phi(p + 1, q    , r + 1),
                                  _phi(p    , q + 1, r + 1),
                                  _phi(p + 1, q + 1, r + 1),
                                  i, 0.0f, k);

    Real phiR1 = bridson::trilerp(_phi(p    , q    , r    ),
                                  _phi(p + 1, q    , r    ),
                                  _phi(p    , q + 1, r    ),
                                  _phi(p + 1, q + 1, r    ),
                                  _phi(p    , q    , r + 1),
                                  _phi(p + 1, q    , r + 1),
                                  _phi(p    , q + 1, r + 1),
                                  _phi(p + 1, q + 1, r + 1),
                                  i, j, 1.0f);

    Real phiR0 = bridson::trilerp(_phi(p    , q    , r    ),
                                  _phi(p + 1, q    , r    ),
                                  _phi(p    , q + 1, r    ),
                                  _phi(p + 1, q + 1, r    ),
                                  _phi(p    , q    , r + 1),
                                  _phi(p + 1, q    , r + 1),
                                  _phi(p    , q + 1, r + 1),
                                  _phi(p + 1, q + 1, r + 1),
                                  i, j, 0.0f);

    Eigen::Vector4f eigenGrad;
    eigenGrad[0] = (float)(phiP1 - phiP0);
    eigenGrad[1] = (float)(phiQ1 - phiQ0);
    eigenGrad[2] = (float)(phiR1 - phiR0);
    eigenGrad[3] = 1.0f;
    
    // Remove the transform at level set creation time and transform the gradient vector by 
    // the current transformation matrix
    Eigen::Matrix4f M1 = m_currentTransformMatrix * m_transformMatrixAtCreation.inverse();
    Eigen::Matrix4f M2 = M1.inverse();
    Eigen::Matrix4f M3 = M1.transpose();
    eigenGrad = M3 * eigenGrad;
    
    grad[0] = eigenGrad[0];
    grad[1] = eigenGrad[1];
    grad[2] = eigenGrad[2];
    //assert(approxEq(grad.norm(), 1.0, 1e-8));
    //m_grad.push_back( eigenGrad );
}

/*
void LevelSet::draw()
{
//    PhysBAM::OPENGL_COLOR_RAMP<float> colorramp  = *(PhysBAM::OPENGL_COLOR_RAMP<float>::Levelset_Color_Linear_Ramp(PhysBAM::OPENGL_COLOR::Red(), PhysBAM::OPENGL_COLOR::Blue(),(float)(_phi.ni+_phi.nj+_phi.nk)));

    glPushAttrib(GL_ENABLE_BIT | GL_POINT_BIT);

    glDisable(GL_LIGHTING);

    //std::cerr<< "level set draw..." << std::endl;

    glPolygonMode(GL_FRONT_AND_BACK, GL_LINE);

    glBegin(GL_QUADS);
    for (unsigned int i=0; i<(_phi.ni-1); ++i)
        for (unsigned int j=0; j<(_phi.nj-1); ++j)
            for (unsigned int k=0; k<(_phi.nk-1); ++k)
            {
                
      //          std::cerr << "_phi( " << i << ", " << j << ", " << k << ") = " << _phi(i,j,k) << std::endl;
                
                if ( _phi( i, j, k ) < 0 )
                {
                    glColor3d( 1.0, 0.0, 0.0 );
                }
                else
                {
                    glColor3d( 1.0, _phi(i,j,k) / 10.0, 1.0 );
                }
//		colorramp.Lookup(_phi(i,j,k)).Send_To_GL_Pipeline();

                float lbb[3] = { _origin[0] + (i    ) * _dx,
                                 _origin[1] + (j    ) * _dx,
                                 _origin[2] + (k    ) * _dx };
                float rbb[3] = { _origin[0] + (i + 1) * _dx,
                                 _origin[1] + (j    ) * _dx,
                                 _origin[2] + (k    ) * _dx };
                float rtb[3] = { _origin[0] + (i + 1) * _dx,
                                 _origin[1] + (j + 1) * _dx,
                                 _origin[2] + (k    ) * _dx };
                float ltb[3] = { _origin[0] + (i    ) * _dx,
                                 _origin[1] + (j + 1) * _dx,
                                 _origin[2] + (k    ) * _dx };
                float lbf[3] = { _origin[0] + (i    ) * _dx,
                                 _origin[1] + (j    ) * _dx,
                                 _origin[2] + (k + 1) * _dx };
                float rbf[3] = { _origin[0] + (i + 1) * _dx,
                                 _origin[1] + (j    ) * _dx,
                                 _origin[2] + (k + 1) * _dx };
                float rtf[3] = { _origin[0] + (i + 1) * _dx,
                                 _origin[1] + (j + 1) * _dx,
                                 _origin[2] + (k + 1) * _dx };
                float ltf[3] = { _origin[0] + (i    ) * _dx,
                                 _origin[1] + (j + 1) * _dx,
                                 _origin[2] + (k + 1) * _dx };

                glVertex3fv(lbb);
                glVertex3fv(rbb);
                glVertex3fv(rtb);
                glVertex3fv(ltb);
                
                glVertex3fv(lbb);
                glVertex3fv(ltb);
                glVertex3fv(ltf);
                glVertex3fv(lbf);
                
                glVertex3fv(lbf);
                glVertex3fv(ltf);
                glVertex3fv(rtf);
                glVertex3fv(rbf);

                glVertex3fv(rbf);
                glVertex3fv(rtf);
                glVertex3fv(rtb);
                glVertex3fv(rbb);

                glVertex3fv(lbb);
                glVertex3fv(lbf);
                glVertex3fv(rbf);
                glVertex3fv(rbb);
                
                glVertex3fv(rtb);
                glVertex3fv(rtf);
                glVertex3fv(ltf);
                glVertex3fv(ltb);
		
            }
    glEnd();

    glPolygonMode(GL_FRONT_AND_BACK, GL_FILL);
    
//    glPointSize( 2.0 );
//    glEnable( GL_POINT_SMOOTH );
//    glBegin( GL_POINTS );
//    glColor3d( 1.0, 0.0, 0.0 );
//    for (size_t p = 0; p < m_transformedRequestPositions.size(); ++p )
//    {
//        glVertex3d( m_transformedRequestPositions[ p ][ 0 ],
//                    m_transformedRequestPositions[ p ][ 1 ],
//                    m_transformedRequestPositions[ p ][ 2 ] );
//    }
//
//    glColor3d( 0.0, 1.0, 0.0 );
//    for (size_t p = 0; p < m_realRequestPositions.size(); ++p )
//    {
//        glVertex3d( m_realRequestPositions[ p ][ 0 ],
//                    m_realRequestPositions[ p ][ 1 ],
//                    m_realRequestPositions[ p ][ 2 ] );
//    }
//    glEnd();
//    glDisable( GL_POINT_SMOOTH );
//
//    glBegin( GL_LINES );
//    glColor3d( 0.0, 0.0, 1.0 );
//    for (size_t p = 0; p < m_grad.size(); ++p )
//    {
//        glVertex3d( m_gradPosition[ p ][ 0 ],
//                    m_gradPosition[ p ][ 1 ],
//                    m_gradPosition[ p ][ 2 ] );
//
//        Eigen::Vector4f end = m_gradPosition[ p ] + m_grad[ p ];
//        glVertex3d( end[ 0 ],
//                    end[ 1 ],
//                    end[ 2 ] );
//    }
//    glEnd();
    
    glPopAttrib();
}

*/

void LevelSet::writeFile(std::fstream &levelSetFile)
{
    levelSetFile.write(reinterpret_cast<char *>(&_origin[0]), sizeof(float) * 3);
    levelSetFile.write(reinterpret_cast<char *>(&_dx),        sizeof(Real));
    levelSetFile.write(reinterpret_cast<char *>(&_phi.ni),    sizeof(int));
    levelSetFile.write(reinterpret_cast<char *>(&_phi.nj),    sizeof(int));
    levelSetFile.write(reinterpret_cast<char *>(&_phi.nk),    sizeof(int));
    levelSetFile.write(reinterpret_cast<char *>(&_phi(0,0,0)),sizeof(float) * _phi.ni * _phi.nj * _phi.nk);

    for (int i=0; i<_phi.ni; ++i)
        for (int j=0; j<_phi.nj; ++j)
            for (int k=0; k<_phi.nk; ++k)
                levelSetFile.write(reinterpret_cast<char *>(&_phiVel(i,j,k)[0]), sizeof(float) * 3);
}

void LevelSet::loadFile(std::fstream &levelSetFile)
{
    levelSetFile.read(reinterpret_cast<char *>(&_origin[0]), sizeof(float) * 3);
    levelSetFile.read(reinterpret_cast<char *>(&_dx),        sizeof(Real));

    int nx, ny, nz;
    levelSetFile.read(reinterpret_cast<char *>(&nx), sizeof(int));
    levelSetFile.read(reinterpret_cast<char *>(&ny), sizeof(int));
    levelSetFile.read(reinterpret_cast<char *>(&nz), sizeof(int));

    _phi.resize(nx, ny, nz);
    levelSetFile.read(reinterpret_cast<char *>(&_phi(0,0,0)), sizeof(float) * _phi.ni * _phi.nj * _phi.nk);

    _phiVel.resize(nx, ny, nz);
    for (int i=0; i<_phi.ni; ++i)
        for (int j=0; j<_phi.nj; ++j)
            for (int k=0; k<_phi.nk; ++k)
                levelSetFile.read(reinterpret_cast<char *>(&_phiVel(i,j,k)[0]), sizeof(float) * 3);
}

void LevelSet::buildLevelSet(const Vec3Indices &triangles,
                             const std::vector<unsigned int> &triIndices,
                                     const std::vector<bridson::Vec3f>  &x,
                                     const std::vector<bridson::Vec3f>  &v,
                                     const bridson::Vec3f &origin,
                                     float dx,
                                     int ni, int nj, int nk,
                                     bridson::Array3f &phi,
			     bridson::Array3<bridson::Vec3f, bridson::Array1<bridson::Vec3f> > &phiVel)
{    
    int exact_band = 1;
    phi.resize(ni, nj, nk);
    phiVel.resize(ni, nj, nk);
    phi.assign((ni+nj+nk)*dx); // upper bound on distance

    //TraceStream(g_log, "LevelSet::buildLevelSet") << "this = " << this << " nijk = " << ni << ", " << nj << ", " << nk << ", Array " << &_phi << " dim " << _phi.ni << ", " << _phi.nj << ", " << _phi.nk << "\n";    

    bridson::Array3i closest_tri(ni, nj, nk, -1);
    bridson::Array3i intersection_count(ni, nj, nk, 0); // intersection_count(i,j,k) is # of tri intersections in (i-1,i]x{j}x{k}
    // we begin by initializing distances near the mesh, and figuring out intersection counts
    //
    bridson::Vec3f ijkmin, ijkmax;
//    for(unsigned int t=0; t<triangles.size(); ++t){
    for (std::vector<unsigned int>::const_iterator tItr=triIndices.begin(); tItr!=triIndices.end(); ++tItr){
        unsigned int t = *tItr;
//        unsigned int p, q, r; assign(triangles[t], p, q, r);
	unsigned int p = triangles[t][0], q= triangles[t][1], r = triangles[t][2];

        // coordinates in grid to high precision
        double fip=((double)x[p][0]-origin[0])/dx, fjp=((double)x[p][1]-origin[1])/dx, fkp=((double)x[p][2]-origin[2])/dx;
        double fiq=((double)x[q][0]-origin[0])/dx, fjq=((double)x[q][1]-origin[1])/dx, fkq=((double)x[q][2]-origin[2])/dx;
        double fir=((double)x[r][0]-origin[0])/dx, fjr=((double)x[r][1]-origin[1])/dx, fkr=((double)x[r][2]-origin[2])/dx;
//        // do distances nearby
        int i0=bridson::clamp(int(bridson::min3(fip,fiq,fir))-exact_band, 0, ni-1), i1=bridson::clamp(int(bridson::max3(fip,fiq,fir))+exact_band+1, 0, ni-1);
        int j0=bridson::clamp(int(bridson::min3(fjp,fjq,fjr))-exact_band, 0, nj-1), j1=bridson::clamp(int(bridson::max3(fjp,fjq,fjr))+exact_band+1, 0, nj-1);
        int k0=bridson::clamp(int(bridson::min3(fkp,fkq,fkr))-exact_band, 0, nk-1), k1=bridson::clamp(int(bridson::max3(fkp,fkq,fkr))+exact_band+1, 0, nk-1);
//        for (int k=0; k<=(nk-1); ++k) for (int j=0; j<=(nj-1); ++j) for (int i=0; i<=(ni-1); ++i){
        for(int k=k0; k<=k1; ++k) for(int j=j0; j<=j1; ++j) for(int i=i0; i<=i1; ++i){
            bridson::Vec3f gx(i*dx+origin[0], j*dx+origin[1], k*dx+origin[2]);
            float t1, t2, t3;
            float d=point_triangle_distance(gx, x[p], x[q], x[r], t1, t2, t3);
            if(d<phi(i,j,k)){
                phi(i,j,k)=d;
                phiVel(i,j,k) = (v[p] * t1 + v[q] * t2 + v[r] * t3);
                closest_tri(i,j,k)=t;
            }
        }
        // and do intersection counts
        j0=bridson::clamp((int)std::ceil(bridson::min3(fjp,fjq,fjr)), 0, nj-1);
        j1=bridson::clamp((int)std::floor(bridson::max3(fjp,fjq,fjr)), 0, nj-1);
        k0=bridson::clamp((int)std::ceil(bridson::min3(fkp,fkq,fkr)), 0, nk-1);
        k1=bridson::clamp((int)std::floor(bridson::max3(fkp,fkq,fkr)), 0, nk-1);

        for(int k=k0; k<=k1; ++k) for(int j=j0; j<=j1; ++j){
            double a, b, c;
            if(point_in_triangle_2d(j, k, fjp, fkp, fjq, fkq, fjr, fkr, a, b, c)){
                double fi=a*fip+b*fiq+c*fir; // intersection i coordinate
                int i_interval=int(std::ceil(fi)); // intersection is in (i_interval-1,i_interval]
                if(i_interval<0) ++intersection_count(0, j, k); // we enlarge the first interval to include everything to the -x direction
                else if(i_interval<ni) ++intersection_count(i_interval,j,k);
                // we ignore intersections that are beyond the +x side of the grid
            }
        }
    }

    // and now we fill in the rest of the distances with fast sweeping
   // omp_set_num_threads(4);
    for(unsigned int pass=0; pass<2; ++pass)
    {
#pragma omp parallel sections
        {
#pragma omp section
            {
                sweep(triangles, x, v, phi, phiVel, closest_tri, origin, dx, +1, +1, +1);
            }
#pragma omp section
            {
                sweep(triangles, x, v, phi, phiVel, closest_tri, origin, dx, -1, -1, -1);
            }
#pragma omp section
            {
                sweep(triangles, x, v, phi, phiVel, closest_tri, origin, dx, +1, +1, -1);
            }
#pragma omp section
            {
                sweep(triangles, x, v, phi, phiVel, closest_tri, origin, dx, -1, -1, +1);
            }
#pragma omp section
            {
                sweep(triangles, x, v, phi, phiVel, closest_tri, origin, dx, +1, -1, +1);
            }
#pragma omp section
            {
                sweep(triangles, x, v, phi, phiVel, closest_tri, origin, dx, -1, +1, -1);
            }
#pragma omp section
            {
                sweep(triangles, x, v, phi, phiVel, closest_tri, origin, dx, +1, -1, -1);
            }
#pragma omp section
            {
                sweep(triangles, x, v, phi, phiVel, closest_tri, origin, dx, -1, +1, +1);
            }
        }
    }

    // then figure out signs (inside/outside) from intersection counts
    for(int k=0; k<nk; ++k) for(int j=0; j<nj; ++j){
        int total_count=0;
        for(int i=0; i<ni; ++i){
            total_count+=intersection_count(i,j,k);
            if(total_count%2==1){ // if parity of intersections so far is odd,
                phi(i,j,k)=-phi(i,j,k); // we are inside the mesh
            }
        }
    }

    //TraceStream(g_log, "LevelSet::buildLevelSet") << "Completed with array " << &_phi << " dim " << _phi.ni << ", " << _phi.nj << ", " << _phi.nk << "\n";    
}

float LevelSet::point_triangle_distance(const bridson::Vec3f &p,
                                                const bridson::Vec3f &a,
                                                const bridson::Vec3f &b,
                                                const bridson::Vec3f &c,
                                                float &t1, float &t2, float &t3)
{
    float ab[3], ac[3], ap[3], bp[3];

    ab[0] = b[0] - a[0];
    ab[1] = b[1] - a[1];
    ab[2] = b[2] - a[2];

    ac[0] = c[0] - a[0];
    ac[1] = c[1] - a[1];
    ac[2] = c[2] - a[2];

    ap[0] = p[0] - a[0];
    ap[1] = p[1] - a[1];
    ap[2] = p[2] - a[2];

    float d1 = ab[0]*ap[0] + ab[1]*ap[1] + ab[2]*ap[2];
    float d2 = ac[0]*ap[0] + ac[1]*ap[1] + ac[2]*ap[2];

    if ((d1 <= 0.0f) && (d2 <= 0.0f))
    {
        t1 = 1.0f;
        t2 = 0.0f;
        t3 = 0.0f;

        return std::sqrt((p[0]-a[0])*(p[0]-a[0]) + (p[1]-a[1])*(p[1]-a[1]) + (p[2]-a[2])*(p[2]-a[2]));
    }

    bp[0] = p[0] - b[0];
    bp[1] = p[1] - b[1];
    bp[2] = p[2] - b[2];

    float d3 = ab[0]*bp[0] + ab[1]*bp[1] + ab[2]*bp[2];
    float d4 = ac[0]*bp[0] + ac[1]*bp[1] + ac[2]*bp[2];

    if ((d3 >= 0.0f) && (d4 <= d3))
    {
        t1 = 0.0f;
        t2 = 1.0f;
        t3 = 0.0f;

        return std::sqrt((p[0]-b[0])*(p[0]-b[0]) + (p[1]-b[1])*(p[1]-b[1]) + (p[2]-b[2])*(p[2]-b[2]));
    }

    float vc = d1*d4 - d3*d2;

    if ((vc <= 0.0f) && (d1 >= 0.0f) && (d3 <= 0.0f))
    {
        float v = d1 / (d1 - d3);

        t1 = 1-v;
        t2 = v;
        t3 = 0;

        float vec[3];
        vec[0] = p[0] - (a[0]+v*ab[0]);
        vec[1] = p[1] - (a[1]+v*ab[1]);
        vec[2] = p[2] - (a[2]+v*ab[2]);

        return std::sqrt(vec[0]*vec[0] + vec[1]*vec[1] + vec[2]*vec[2]);
    }

    float cp[3];
    cp[0] = p[0] - c[0];
    cp[1] = p[1] - c[1];
    cp[2] = p[2] - c[2];

    float d5 = ab[0]*cp[0] + ab[1]*cp[1] + ab[2]*cp[2];
    float d6 = ac[0]*cp[0] + ac[1]*cp[1] + ac[2]*cp[2];

    if ((d6 >= 0.0f) && (d5 <= d6))
    {
        t1 = 0;
        t2 = 0;
        t3 = 1;

        return std::sqrt((p[0]-c[0])*(p[0]-c[0]) + (p[1]-c[1])*(p[1]-c[1]) + (p[2]-c[2])*(p[2]-c[2]));
    }

    float vb = d5*d2 - d1*d6;

    if ((vb <= 0.0f) && (d2 >= 0.0f) && (d6 <= 0.0f))
    {
        float w = d2 / (d2 - d6);

        t1 = 1-w;
        t2 = 0;
        t3 = w;

        float vec[3];
        vec[0] = p[0] - (a[0]+w*ac[0]);
        vec[1] = p[1] - (a[1]+w*ac[1]);
        vec[2] = p[2] - (a[2]+w*ac[2]);

        return std::sqrt(vec[0]*vec[0] + vec[1]*vec[1] + vec[2]*vec[2]);
    }

    float va = d3*d6 - d5*d4;

    if ((va <= 0.0f) && ((d4-d3) >= 0.0f) && ((d5-d6) >= 0.0f))
    {
        float w = (d4 - d3) / ((d4 - d3) + (d5 - d6));

        t1 = 0;
        t2 = 1-w;
        t3 = w;

        float vec[3];
        vec[0] = p[0] - (b[0]+w*(c[0]-b[0]));
        vec[1] = p[1] - (b[1]+w*(c[1]-b[1]));
        vec[2] = p[2] - (b[2]+w*(c[2]-b[2]));

        return std::sqrt(vec[0]*vec[0] + vec[1]*vec[1] + vec[2]*vec[2]);
    }

    float denom = 1.0f / (va + vb + vc);
    float v = vb * denom;
    float w = vc * denom;
    float u = 1.0f - v - w;

    t1 = u;
    t2 = v;
    t3 = w;

    float vec[3];
    vec[0] = p[0] - (u*a[0] + v*b[0] + w*c[0]);
    vec[1] = p[1] - (u*a[1] + v*b[1] + w*c[1]);
    vec[2] = p[2] - (u*a[2] + v*b[2] + w*c[2]);

    return std::sqrt(vec[0]*vec[0] + vec[1]*vec[1] + vec[2]*vec[2]);
}

void LevelSet::check_neighbour(const Vec3Indices &tri,
                                       const std::vector<bridson::Vec3f> &x,
                                       const std::vector<bridson::Vec3f> &v,
                                       bridson::Array3f &phi,
                                     bridson::Array3<bridson::Vec3f,bridson::Array1<bridson::Vec3f> > &phi_vel,
                                       bridson::Array3i &closest_tri,
                                       const bridson::Vec3f &gx,
                                       int i0, int j0, int k0, int i1, int j1, int k1)
{
    if(closest_tri(i1,j1,k1)>=0){
//        unsigned int p, q, r; assign(tri[closest_tri(i1,j1,k1)], p, q, r);
	unsigned int p = tri[closest_tri(i1,j1,k1)][0];
	unsigned int q = tri[closest_tri(i1,j1,k1)][1];
	unsigned int r = tri[closest_tri(i1,j1,k1)][2];

        float t1, t2, t3;
        float d=point_triangle_distance(gx, x[p], x[q], x[r], t1, t2, t3);
        if(d<phi(i0,j0,k0)){
            phi(i0,j0,k0)=d;
            phi_vel(i0,j0,k0) = (v[p] * t1 + v[q] * t2 + v[r] * t3);
            closest_tri(i0,j0,k0)=closest_tri(i1,j1,k1);
        }
    }
}

void LevelSet::sweep(const Vec3Indices &tri,
                             const std::vector<bridson::Vec3f> &x,
                             const std::vector<bridson::Vec3f> &v,
                             bridson::Array3f &phi,
                             bridson::Array3<bridson::Vec3f, bridson::Array1<bridson::Vec3f> > &phi_vel,
                             bridson::Array3i &closest_tri,
                             const bridson::Vec3f &origin,
                             float dx,
                             int di, int dj, int dk)
{
    int i0, i1;
    if(di>0){ i0=1; i1=phi.ni; }
    else{ i0=phi.ni-2; i1=-1; }
    int j0, j1;
    if(dj>0){ j0=1; j1=phi.nj; }
    else{ j0=phi.nj-2; j1=-1; }
    int k0, k1;
    if(dk>0){ k0=1; k1=phi.nk; }
    else{ k0=phi.nk-2; k1=-1; }
    for(int k=k0; k!=k1; k+=dk) for(int j=j0; j!=j1; j+=dj) for(int i=i0; i!=i1; i+=di){
        bridson::Vec3f gx(i*dx+origin[0], j*dx+origin[1], k*dx+origin[2]);
        check_neighbour(tri, x, v, phi, phi_vel, closest_tri, gx, i, j, k, i-di, j,    k);
        check_neighbour(tri, x, v, phi, phi_vel, closest_tri, gx, i, j, k, i-di, j-dj, k);
        check_neighbour(tri, x, v, phi, phi_vel, closest_tri, gx, i, j, k, i,    j,    k-dk);
        check_neighbour(tri, x, v, phi, phi_vel, closest_tri, gx, i, j, k, i-di, j,    k-dk);
        check_neighbour(tri, x, v, phi, phi_vel, closest_tri, gx, i, j, k, i,    j-dj, k-dk);
        check_neighbour(tri, x, v, phi, phi_vel, closest_tri, gx, i, j, k, i-di, j-dj, k-dk);
    }
}

int LevelSet::orientation(double x1, double y1, double x2, double y2, double &twice_signed_area)
{
    twice_signed_area=y1*x2-x1*y2;
    if(twice_signed_area>0) return 1;
    else if(twice_signed_area<0) return -1;
    else if(y2>y1) return 1;
    else if(y2<y1) return -1;
    else if(x1>x2) return 1;
    else if(x1<x2) return -1;
    else return 0; // only true when x1==x2 and y1==y2
}

bool LevelSet::point_in_triangle_2d(double x0, double y0, 
                                            double x1, double y1, double x2, double y2, double x3, double y3,
                                            double& a, double& b, double& c)
{
    x1-=x0; x2-=x0; x3-=x0;
    y1-=y0; y2-=y0; y3-=y0;
    int signa=orientation(x2, y2, x3, y3, a);
    if(signa==0) return false;
    int signb=orientation(x3, y3, x1, y1, b);
    if(signb!=signa) return false;
    int signc=orientation(x1, y1, x2, y2, c);
    if(signc!=signa) return false;
    double sum=a+b+c;
    assert(sum!=0); // if the SOS signs match and are nonzero, there's no way all of a, b, and c are zero.
    a/=sum;
    b/=sum;
    c/=sum;
    return true;
}

}
