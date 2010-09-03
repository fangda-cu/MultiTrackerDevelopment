// CollisionMeshData.cc
//

#include "CollisionMeshData.hh"

#include <sstream>
#include <map>
#include <GL/gl.h>

namespace BASim {

CollisionMeshData::CollisionMeshData()
    :/*_levelsetDx(0.0),*/  _thickness(1.0), _friction(0.0), _fullCollisions(false), _initialized(false),_fps(24)
{
    //_bvTree(TRIANGLE_TREE),
    //std::cout<<"Constructing... ";
 
    //_phiPrevious = new LevelSet;
    //_phiCurrent = new LevelSet;
    //_phiPrevious = new AdaptiveLevelSet;
    //_phiCurrent = new AdaptiveLevelSet;
    //std::cout<<"complete"<<std::endl;
  
  
  // To export mesh data into files so as to run in BASimulator
  recordToFile = false;
  recordFrames = 0;
  
  //recordFilename = "/local1/jjoo/scenes/sac2/sac_mesh.";
  recordFilename = "/local1/jjoo/scenes/sak_full/mesh.";

}

CollisionMeshData::~CollisionMeshData()
{
    //delete _phiPrevious;
    //delete _phiCurrent;
}

void CollisionMeshData::initialize()
{
  
  std::cout << "CollisionMeshData::initialize()...\n";
   // if (!_initialized)
    {
        _nbrTriangles = triangleIndices.size() / 3;
 //       _x.resize(_nbrTriangles);
 //       _v.resize(_nbrTriangles);

        // Centroids of the triangles are used to build a bounding
        // volume hierarchy
        //
        Positions centroids(_nbrTriangles);

        // BVH code needs a list of all triangles contained in
        // the BVH
        //
        _triIndices.resize(_nbrTriangles);

        // A copy of triangleIndices, but isntead of a flat array
        // it is stored as a vector of Vec3ui, as needed by the level set code
        //
        _tri.resize(_nbrTriangles);

        // Stores the indices of the edges that make up each triangle
        //
        _triangleEdgeIndices.resize(_nbrTriangles);

        std::map<std::pair<uint, uint>, uint> edgeMap;
        for (size_t i=0; i<_nbrTriangles; ++i)
        {
            Vec3d v1 = triangleVertex(i, 0);
            Vec3d v2 = triangleVertex(i, 1);
            Vec3d v3 = triangleVertex(i, 2);

            centroids[i] = (v1 + v2 + v3) / 3.0;

            _triIndices[i] = i;

            _tri[i][0] = triangleIndices[3 * i];
            _tri[i][1] = triangleIndices[3 * i + 1];
            _tri[i][2] = triangleIndices[3 * i + 2];

            // For each edge in this triangle, see if its already been referenced
            // otherwise increment the index counter
            //
            for (uint j=0; j<3; ++j)
            {
                // The two vertices that make up this edge, have they been
                // seen before?
                //
                uint idx1 = _tri[i][j];
                uint idx2 = _tri[i][(j+1)%3];

                if (idx2 < idx1)
                    std::swap(idx1, idx2);

                std::map<std::pair<uint,uint>, uint>::iterator itr = edgeMap.find(std::make_pair(idx1,idx2));
                if (itr == edgeMap.end())
                {
                    // This edge is new, aincrement the counter and flag it
                    //
                    _edgeIndices.push_back(idx1);
                    _edgeIndices.push_back(idx2);
                    uint idx = (_edgeIndices.size() / 2) - 1;
                    edgeMap.insert(std::make_pair(std::make_pair(idx1, idx2), idx));

                    _triangleEdgeIndices[i][j] = idx;
                }
                else
                    _triangleEdgeIndices[i][j] = itr->second;
            }
        }
	
	
//	buildLevelSet();
//      _bvTree.buildTree(centroids, _triIndices);

        _initialized = true;
        
        recordFrames = 0;
        writeMeshesToFile();
        
    }
   // std::cout<<"CollisionMeshData::initialize()...complete"<<std::endl;
}

void CollisionMeshData::writeMeshesToFile() {
  if (recordToFile) {
    std::cout << "recording meshes\n";
    
    std::stringstream ss;
    std::string filename = recordFilename;
    
    ss << recordFilename;
    ss << recordFrames;
    ss >> filename;
    ss >> ".mesh";
    
    std::cout << filename << "\n";
    
    FILE *fp;
    fp = fopen( filename.c_str(), "w" );
    if ( fp == NULL )
    {
      std::cout << "can't open file\n";
      return;
    }

    size_t numberOfVerts = currPositions.size();
    size_t numberOfFaces = _nbrTriangles;
    fwrite( &numberOfVerts, sizeof( size_t ), 1, fp );
    fwrite( &numberOfFaces, sizeof( size_t ), 1, fp );

    for ( size_t r=0; r<numberOfVerts; r++ )
    {
      double pos[3];
      
      Vec3d vertex = oldPositions[r];

      pos[ 0 ] = vertex.x();
      pos[ 1 ] = vertex.y();
      pos[ 2 ] = vertex.z();

      fwrite( &pos[0], sizeof( double ), 3, fp );
    }
    
    for ( size_t r=0; r<numberOfFaces; r++ )
    {
      uint idx[3];

      idx[ 0 ] = _tri[r][0];
      idx[ 1 ] = _tri[r][1];
      idx[ 2 ] = _tri[r][2];

      fwrite( &idx[0], sizeof( uint ), 3, fp );
    }
    
    fclose ( fp );
    
    recordFrames++;
  }
}

void CollisionMeshData::clearAll()
{
    triangleIndices.clear();
    currPositions.clear();
    velocities.clear();
}

void CollisionMeshData::reset(vector<Vec3d>& points)
{

   // if (!_initialized)
        initialize();

    for (size_t currVertex=0; currVertex<points.size(); ++currVertex)
    {
        for (size_t i=0; i<3; ++i)
        {
            oldPositions[currVertex][i]  = newPositions[currVertex][i] =
            currPositions[currVertex][i] = points[currVertex][i];
            velocities[currVertex][i]    = 0.0;
        }
    } 

    updateGrid(points);     
//    _bvTree.fitTree(oldPositions, triangleIndices);
    //buildLevelSet();

}

void CollisionMeshData::update(vector<Vec3d>& points, std::string filename, int currFrame)
{
    // This should probably be an assertion
    //
    if (!_initialized)
    {
        cerr << "Cannot update collision data before it has been initialised, please rewind simulation to start time\n"; 
        //initialize();
        return;
    }
        

    for (size_t currVertex=0; currVertex<points.size(); ++currVertex)
    {
        for (size_t i=0; i<3; ++i)
        {
            oldPositions[currVertex][i] = currPositions[currVertex][i] = newPositions[currVertex][i];
            newPositions[currVertex][i] = points[currVertex][i];

            // TODO: Replace hard-coded 24 with fps as given by user
            //
            velocities[currVertex][i] = (newPositions[currVertex][i] - oldPositions[currVertex][i]) * _fps;
        }
    }    
    
    updateGrid(points, filename);
    
    writeMeshesToFile();
    
 //   buildLevelSet();
}

void CollisionMeshData::updateGrid(vector<Vec3d>& points, std::string filename)
{
    // If the file doesn't exist, it either means we're not caching or this cache
    // file doesn't exist, either way we need to compute the collision data
    //
    if (!_grid.readFile(filename))
    {
        //_bvTree.fitTree(oldPositions, newPositions, triangleIndices, _thickness);

        // All this does it basically find the bounds for the grid and the "optimal"
        // cell size
        //
        std::vector<Vec3d> xmins(_nbrTriangles);
        std::vector<Vec3d> xmaxs(_nbrTriangles);
        for (uint i=0; i<_nbrTriangles; ++i)
        {
            Vec3i tri = _tri[i];

            minmax(oldPositions[tri[2]],
                   newPositions[tri[2]],
                   oldPositions[tri[1]],
                   newPositions[tri[1]],
                   oldPositions[tri[0]],
                   newPositions[tri[0]],
                   xmins[i], xmaxs[i]);
        }

        Vec3d xmax = xmaxs[0];
        Vec3d xmin = xmins[0];
        Real maxdistance = 0.0;

        unsigned int n = xmins.size();
        for (unsigned int i=0; i<n; i++)
        {
            update_minmax(xmins[i], xmin, xmax);
            update_minmax(xmaxs[i], xmin, xmax);
            maxdistance = std::max(maxdistance, (xmaxs[i] - xmins[i]).norm());
        }

        for (unsigned int i=0; i<3; i++)
        {
            xmin[i] -= 2 * maxdistance + 1e-6;
            xmax[i] += 2 * maxdistance + 1e-6;
        }

        Vec3i dims(1,1,1);

        unsigned int maxdim = (unsigned int)std::sqrt(n);

        if ((xmax-xmin).norm() > 1e-6)
        {
            double elements_per_cell = 0.1;
            double volume = (xmax[0] - xmin[0]) * (xmax[1] - xmin[1]) * (xmax[2] - xmin[2]);
            double volume_per_cell = volume * elements_per_cell / n;
            double ideal_cell_size = std::pow(volume_per_cell, 1.0/3.0);

            for (unsigned int i=0; i<3; i++)
            {
                unsigned int d = (unsigned int)std::ceil((xmax[i] - xmin[i]) / ideal_cell_size);
                if (d < 1)
                    d = 1;
                if (d > maxdim)
                    d = maxdim;

                dims[i] = d;
            }
        }

        _grid.set(dims, xmin, xmax);

        // Insert each triangle into the grid
        //
        for (unsigned int i=n; i>0; i--)
            _grid.addElement(i-1, xmins[i-1], xmaxs[i-1]);

        if (!filename.empty())
        {
            std::cout << "Writing " << filename << "...";
            std::cout.flush();
            _grid.writeFile(filename);
            std::cout << "Done" << std::endl;
        }
    }
}


void CollisionMeshData::interpolate(Real percentage)
{    
  //  allPositions.resize(allPositions.size()+1);
   // cerr << "interpolating and allPositions.size() = " << allPositions.size() << endl;
    for (size_t i=0; i<currPositions.size(); ++i)
    {
        prevPositions[i] = currPositions[i];
        currPositions[i] = (1.0 - percentage) * oldPositions[i] + percentage * newPositions[i];

     //   allPositions[allPositions.size()-1][i] = currPositions[i];
    }
//    _bvTree.fitTree(prevPositions, currPositions, triangleIndices);

    _percent = percentage;
}

void CollisionMeshData::draw()
{
    // Useful for debugging
    //
     //if(_phiCurrent->isInitialized())
	//_phiCurrent->draw();
         
    //for ( size_t p=0; p<currPositions.size(); p++ )     
    {
      glColor3f(0, 1, 0);
      glBegin(GL_TRIANGLES);
      for (size_t i=0; i<_nbrTriangles; ++i)
      {    
          float c = (float)(i)/(float)(_nbrTriangles);          
          
          glVertex3dv(currPositions[triangleIndices[(3 * i)    ]].data());
          glVertex3dv(currPositions[triangleIndices[(3 * i) + 1]].data());
          glVertex3dv(currPositions[triangleIndices[(3 * i) + 2]].data());
      }
      glEnd();
      
      glLineWidth(5.0);
      glColor3f(0, 0, 0);          
      glBegin(GL_LINES);
      for (size_t i=0; i<_nbrTriangles; ++i)
      {    
          float c = (float)(i)/(float)(_nbrTriangles);
          
          glVertex3dv(currPositions[triangleIndices[(3 * i)    ]].data());
          glVertex3dv(currPositions[triangleIndices[(3 * i) + 1]].data());
          
          glVertex3dv(currPositions[triangleIndices[(3 * i) + 1]].data());
          glVertex3dv(currPositions[triangleIndices[(3 * i) + 2]].data());
          
          glVertex3dv(currPositions[triangleIndices[(3 * i) + 2]].data());
          glVertex3dv(currPositions[triangleIndices[(3 * i)    ]].data());
      }
      glEnd();
      glLineWidth(1.0);
    }
}

void CollisionMeshData::sizeLevelSet(Vec3d &origin,Vec3i &dims, Real &dx, Real length[3])
{
   /* Vec3d xmin,xmax,dX;
    _grid.getGridDims(xmin,xmax,dX);
    Real mindx = dX[0];
    for(uint i=0; i<3; i++)
    {
	origin[i]=xmin[i];
	mindx =std:: min(mindx, dX[i]);
	length[i]=xmax[i]-xmin[i];
    }
    if(_levelsetDx)
	dx = _levelsetDx;
    else
	dx = mindx; 

    for(uint i=0; i<3; i++)
	dims[i] = (int)ceil(length[i] / dx);*/
}

void CollisionMeshData::buildLevelSet()
{   
    /*std::cout<<"building level set ...";

    for (int currVertex=0; currVertex<currPositions.size(); ++currVertex)
    {
        for (int i=0; i<3; ++i)
        {
            _x[currVertex][i] = currPositions[currVertex][i];
            _v[currVertex][i] = velocities[currVertex][i];
        }
    }

    bridson::Vec3f origin;
    
    Real dx;
    Vec3ui dims;
    Real length[3];
    sizeLevelSet(origin,dims,dx,length);
    

    if(_phiCurrent->isInitialized())
    {
	delete _phiPrevious;
	_phiPrevious = _phiCurrent;
	_phiCurrent = new LevelSet;
    }


    _phiCurrent->buildLevelSet(_tri, _triIndices, _x, _v, origin, length, dx, dims[0], dims[1], dims[2], _nbrTriangles);
    std::cout<<"Complete!"<<std::endl;
*/
}

Real CollisionMeshData::getLevelSetValue(Vec3d& x, Vec3d& v)
{
    /*Real distPrev;
    Real distCurr;
    distCurr = _phiCurrent->getLevelSetValueVelocity(x, v);
    if (_phiPrevious->isInitialized())
        distPrev = _phiPrevious->getLevelSetValue(x);
    else
        distPrev = distCurr;

    return (bridson::lerp(distPrev, distCurr, _percent));*/
    return 0;
}


Vec3d& CollisionMeshData::vertex(uint i)
{
    return (currPositions[i]);
}

const Vec3d& CollisionMeshData::vertex(uint i) const
{
    return (currPositions[i]);
}

Vec3d& CollisionMeshData::velocity(uint i)
{
    return (velocities[i]);
}

const Vec3d& CollisionMeshData::velocity(uint i) const
{
    return (velocities[i]);
}

Vec3d& CollisionMeshData::triangleVertex(uint faceID, uint i)
{
    return vertex(triangleIndices[3 * faceID + i]);
}

const Vec3d& CollisionMeshData::triangleVertex(uint faceID, uint i) const
{
    return vertex(triangleIndices[3 * faceID + i]);
}

Vec3d& CollisionMeshData::triangleVelocity(uint faceID, uint i)
{
    return velocity(triangleIndices[3 * faceID + i]);
}

const Vec3d& CollisionMeshData::triangleVelocity(uint faceID, uint i) const
{
    return velocity(triangleIndices[3 * faceID + i]);
}

}

