#include "VolumetricCollisionsGPU.hh"
#include "VolumetricCollisionsGPU.cuh"

#include <sys/time.h>

#include <ocuequation/sol_projectvardiv3d.h>

VolumetricCollisionsGPU::VolumetricCollisionsGPU()
    :_initialised(false),_gridSortBits(18),printdbg(true),_radius(1.0)
{}

VolumetricCollisionsGPU::VolumetricCollisionsGPU(Rods &rods)
    :_initialised(false),_gridSortBits(18),printdbg(true),_radius(1.0)
{
    reinitialiseRods(rods);
}

VolumetricCollisionsGPU::~VolumetricCollisionsGPU()
{

    deleteRods();
}

void VolumetricCollisionsGPU::reinitialiseRods(Rods &rods)
{
    deleteRods();
    uint rodsn = rods.size();

    int addedVertices = 0;
    // compute number of edges and vertices
    //
    for(RodsIterator rItr = rods.begin(); rItr!=rods.end(); ++rItr)
    {
        Rod *rod = *rItr;
        _numEdges += rod->ne();
        _numVertices += rod->nv();
        
        for(uint i=0; i< rod->ne(); i++)
        {
            if(rod->edgeLen(i) > _radius)
                addedVertices += ceil(rod->edgeLen(i) / _radius) -1;
        }
    }
    _numOriginalVertices = _numVertices;
    _numVertices += addedVertices;
    
    // allocate all edge-based arrays
    //
    _hPos = new float[_numVertices*4];
    _hVel = new float[_numVertices*4];
    
    allocateArray((void**)&_dPos, _numVertices*4*sizeof(float));
    allocateArray((void**)&_dVel, _numVertices*4*sizeof(float));
    allocateArray((void**)&_dSortedPos, _numVertices*4*sizeof(float));
    allocateArray((void**)&_dSortedVel, _numVertices*4*sizeof(float));
    
    allocateArray((void**)&_dGridHash, _numVertices*sizeof(uint));
    allocateArray((void**)&_dGridIndex, _numVertices*sizeof(uint));
    
    if(printdbg)
        std::cout<<"\nFound "<<_numVertices<<" vertices and including "<<addedVertices<<" added vertices.\n"<<std::endl;
    _initialised = true;
}



void VolumetricCollisionsGPU::deleteRods()
{
    if(_initialised)
    {
        delete [] _hPos;
        delete [] _hVel;

        freeArray(_dPos);
        freeArray(_dVel);
        freeArray(_dSortedPos);
        freeArray(_dSortedVel);
        freeArray(_dGridHash);
        freeArray(_dGridIndex);
    }
    _numEdges = 0;
    _numVertices = 0;
    _numOriginalVertices = 0;
    _initialised = false;
}

void VolumetricCollisionsGPU::respondVolumetricCollisions(Rods &rods, Real targetEdgeDensity, 
                                                       Real volumetricRadius, Real gridDx,
                                                       Vec3<Real> separationCondition,
                                                       CollisionMeshDataHashMap &collisionMeshes)
{
    if(rods.size()==0)
        return;

    bool printdbg = true;
    double rsq = volumetricRadius*volumetricRadius;
    double oneOverR=1.0/volumetricRadius;
    double oneOverR4 = 1.0/(rsq*rsq);
    
    timeval t1,t2;
    double elapsedTime;

    //volume parameters
    float minCorner[3];
    uint gridSize[3];    
    
    gettimeofday(&t1,NULL);
    if(_radius != volumetricRadius || !_initialised)
        reinitialiseRods(rods);
   
    _radius = volumetricRadius;   
 
    // update positions/velocities in curves (use saved positions)
    // find bounding box for hairs
    //
    int vid = 0;
    
    RodsIterator rItr=rods.begin();
    Rod *rod = *rItr;
    Vec3r minc= rod->getStartPositions()[0];
    Vec3r maxCorner = minc;

    float flip = rod->getFlip(); //assume constant over all rods
    float slip = rod->getSlip(); //assume constant over all rods

    if(slip == 1)
    { 
        std::cout<<"slip = 1. exiting volume\n";
        return;
    }

    for(uint i=0;i<3; i++)
        minCorner[i] = minc[i];
    
    for (; rItr!=rods.end(); ++rItr)
    {
        rod = *rItr;
        for (uint i=0; i<rod->nv(); ++i)
        {
            Vec3d x = rod->getStartPositions()[i];
            Vec3d v = rod->getVelocities()[i];
            for(uint k=0; k<3; ++k)
            {
                _hPos[vid*4 + k] = x[k];
                _hVel[vid*4 + k] = v[k];
                minCorner[k] = min(minCorner[k], (float)x[k]);
                maxCorner[k] = max(maxCorner[k], x[k]);
            }
            vid++;
        }
        
        for(uint i=0; i<rod->ne(); ++i)
        {
            if(rod->edgeLen(i) > _radius)
            {
                int addedVertices = ceil(rod->edgeLen(i) / _radius);
                float  b = 1.0 / (float)addedVertices;
                addedVertices--;
                Vec3r x0 = rod->getStartPositions()[rod->edgeIndex(i,0)];
                Vec3r x1mx0 = rod->getStartPositions()[rod->edgeIndex(i,1)] - x0;
                Vec3r v0 = rod->getVelocities()[rod->edgeIndex(i,0)];
                Vec3r v1mv0 = rod->getVelocities()[rod->edgeIndex(i,1)] - v0; 
                for( uint j = 1; j<= addedVertices; j++)
                {
                    for(uint k=0; k<3; ++k)
                    {
                        _hPos[vid*4 + k] = x0[k] + (float)j * b * x1mx0[k];
                        _hVel[vid*4 + k] = v0[k] + (float)j * b * v1mv0[k];
                    }
                    vid++;
                }
            }
        }       
    }
    

    if(printdbg)
        std::cout<<"\nread "<<vid<<" vertices.\n minCorner = ["<<minCorner[0]<<", "<<minCorner[1]<<", "<<minCorner[2]<<"] maxCorner = "<<maxCorner<<std::endl;

    
    // copy _hPos and _hVel to device and bind to appropriate textures
    //
    copyArrayToDevice(_dPos,_hPos,0,_numVertices*4*sizeof(float));
    copyArrayToDevice(_dVel,_hVel,0,_numVertices*4*sizeof(float));
    setPosTexFromDevice(_dPos,_numVertices);
    setVelTexFromDevice(_dVel,_numVertices);


    // compute grid dims
    //
    double oneOverDx = 1.0/gridDx;
    for(int d=0; d<3; d++)
    {
        minCorner[d] = std::floor((minCorner[d]-volumetricRadius)*oneOverDx - 2.0);
        maxCorner[d] = std::floor((maxCorner[d]+volumetricRadius)*oneOverDx + 2.0);
        gridSize[d] = (uint)(maxCorner[d] - minCorner[d]);
        minCorner[d] *= gridDx;
        maxCorner[d] *= gridDx;
    }
    
    // allocate grid based arrays
    //

    std::vector<ocu::Grid3DDimension> faceDimensions(3);
    faceDimensions[0].init(gridSize[0]+1,gridSize[1], gridSize[2],1,1,1);
    faceDimensions[1].init(gridSize[0],gridSize[1]+1, gridSize[2],1,1,1);
    faceDimensions[2].init(gridSize[0],gridSize[1], gridSize[2]+1,1,1,1);
    ocu::Grid3DDimension::pad_for_congruence(faceDimensions);
 
    int numCells = gridSize[0]*gridSize[1]*gridSize[2];
    if(printdbg)
        std::cout<<"There are "<<numCells<<" cells in grid"<<std::endl;

    if(numCells > std::pow(2.0,(int)_gridSortBits)){
        std::cout<<"Grid has too many cells for sort.  Increase gridDx or increase _gridSortBits"<<std::endl;
        return;
    }


    std::vector<Grid3DDeviceF*> dfaceWeights;
    std::vector<Grid3DDeviceF*> dinitialVelocities;
    std::vector<Grid3DDeviceF*> dfinalVelocities;

    // use these to read in collision object data
    //
    std::vector<ocu::Grid3DHostF*> hinitialVelocities;
    std::vector<ocu::Grid3DHostB*> hpsiN;

    for(uint i=0; i<3; i++){
        Grid3DDeviceF* fw = new Grid3DDeviceF;
        Grid3DDeviceF* ivel = new Grid3DDeviceF;
        Grid3DDeviceF* fvel = new Grid3DDeviceF;

        ocu::Grid3DHostF* hvel = new ocu::Grid3DHostF;
        ocu::Grid3DHostB* hp = new ocu::Grid3DHostB;

        fw->init_congruent(faceDimensions[i]);
        ivel->init_congruent(faceDimensions[i]);
        fvel->init_congruent(faceDimensions[i]);

        hvel->init_congruent(faceDimensions[i]);
        hp->init_congruent(faceDimensions[i]);
        
        fw->clear_zero();
        ivel->clear_zero();
        fvel->clear_zero();

        hvel->clear_zero();
        hp->clear_zero();
        
        dfaceWeights.push_back(fw);
        dinitialVelocities.push_back(ivel);
        dfinalVelocities.push_back(fvel);
        
        hinitialVelocities.push_back(hvel);
        hpsiN.push_back(hp);
    }

    ocu::Sol_ProjectVariableDivergence3DDeviceF projection;
    projection.initialize_storage(gridSize[0], gridSize[1], gridSize[2], 
                                  gridDx, gridDx, gridDx, 
                                  dfinalVelocities[0],dfinalVelocities[1],dfinalVelocities[2],false);

    Grid3DDeviceF* dCellWeights = new Grid3DDeviceF;
    dCellWeights->init_congruent(projection.divergence);
    dCellWeights->clear_zero();

    std::vector<Grid3DDeviceF*> dCellVelocities;
    for(uint i=0;i<3;i++){
        Grid3DDeviceF* cv = new Grid3DDeviceF;
        cv->init_congruent(projection.divergence);
        cv->clear_zero();
        dCellVelocities.push_back(cv);
    }

    uint *dCellStart, *dCellEnd;
    allocateArray((void**)&dCellStart,numCells*sizeof(uint));
    allocateArray((void**)&dCellEnd,numCells*sizeof(uint));

    if(printdbg)
        std::cout<<"Initialized "<<gridSize[0]<<" by "<<gridSize[1]<<" by "<<gridSize[2]<<" grid"<<std::endl;

     // store parameters to constant memory
    //
    setVolumeParams(minCorner, gridSize, gridDx, volumetricRadius, targetEdgeDensity);

    gettimeofday(&t2, NULL);
    elapsedTime = ( ( t2.tv_sec * 1000000 + t2.tv_usec ) - ( t1.tv_sec * 1000000 + t1.tv_usec) ) / 1000.0;
    if(printdbg)
        std::cout<<"Total initialization time: "<<elapsedTime<<std::endl;

    gettimeofday(&t1,NULL);

    // read in collision object data
    //
    for(int axis = 0; axis < 3; axis++)
        for(int i=0; i<faceDimensions[axis].nx(); i++)
            for(int j=0; j<faceDimensions[axis].ny(); j++)
                for(int k=0; k<faceDimensions[axis].nz(); k++)
                {
                    Vec3r face,v;
                    face[0]=(i+.5)*gridDx + minCorner[0];
                    face[1]=(j+.5)*gridDx + minCorner[1];
                    face[2]=(k+.5)*gridDx + minCorner[2];
                    face[axis]-=.5*gridDx;
                    
                    bool insideCollisionObject=false;
                    CollisionMeshData *cmData;
                    for(CollisionMeshDataHashMapIterator cmItr = collisionMeshes.begin();
                        cmItr != collisionMeshes.end(); cmItr++)
                    {
                        cmData = cmItr->second;
                        if(cmData->getLevelSetValue(face,v) <= 0.0)
                        {
                            insideCollisionObject = true;
                            break;
                        }
                    }
                    if(insideCollisionObject)
                    {
                        hpsiN[axis]->at(i,j,k)=true;
                        hinitialVelocities[axis]->at(i,j,k) = v[axis];
                    }
                }


    // copy collision data to device
    //
    for(int i=0; i<3; i++)
    {
        dinitialVelocities[i]->copy_all_data(*hinitialVelocities[i]);
        projection.psiN[i]-> copy_all_data(*hpsiN[i]);
    }
    
      
    // calculate grid hash
    //
    calcVertexHash(_dGridHash, _dGridIndex,_numVertices);
    
    
    // use radix sort to sort grid and edge indices
    //
    radixSort(_dGridHash, _dGridIndex, _numVertices, _gridSortBits);
    
    // reorder particle arrays into sorted order and find start and end of each cell
    // bind pos and vel tex to sorted arrays
    //
    reorderVertexDataAndFindCellStart(dCellStart, dCellEnd, 
                                      _dSortedPos, _dSortedVel, 
                                      _dGridHash, _dGridIndex,
                                      _numVertices, gridSize,numCells);
    
    // these have been copied and bound to textures, so can free them
    //
    freeArray(dCellStart);
    freeArray(dCellEnd);
    
    // rasterize to grid
    //
    rasterizeVerticesToCells(dCellWeights, dCellVelocities, &projection.divergence, &projection.divergenceMultiplier,projection.psiD,gridSize, numCells);
    
    rasterizeVerticesToFaces(dfaceWeights, dinitialVelocities, projection.psiN, gridSize, numCells);

    // apply separation condition
    //
    float separationF[3];
    for(int i=0;i<3;i++)
        separationF[i] = separationCondition[i];
    applySeparationCondition(&projection.divergence, &projection.divergenceMultiplier, projection.psiD, projection.psiN, dCellWeights, dCellVelocities, dinitialVelocities, separationF,gridSize, numCells);

    if(printdbg)
    {

        float maxWeight;
        dCellWeights->reduce_max(maxWeight);
        std::cout<<"Target cell weight: " <<targetEdgeDensity<<std::endl;
        std::cout<<"Maximum cell weight: "<<maxWeight<<std::endl;

    }


    for(uint i=0; i<3; i++)
        dfinalVelocities[i]->copy_all_data(*dinitialVelocities[i]);

    gettimeofday(&t2, NULL);
    elapsedTime = ( ( t2.tv_sec * 1000000 + t2.tv_usec ) - ( t1.tv_sec * 1000000 + t1.tv_usec) ) / 1000.0;
    if(printdbg)
        std::cout<<"Total rasterization time: "<<elapsedTime<<std::endl;
            
    // project
    //
    projection.solve((float)1e-5,(float)1e-7, 1000);


    gettimeofday(&t1,NULL);

    // find change in velocities (write over dinitialVelocities)    
    for(uint i=0; i<3; i++)
        dinitialVelocities[i]->linear_combination(-1.0f, *dinitialVelocities[i], 1.0f, *dfinalVelocities[i]);


    // interpolate velocities back
    //
    interpolateToVertices(_dVel, _dPos, slip, flip, _numOriginalVertices, dfinalVelocities, dinitialVelocities,projection.psiD);

    copyArrayFromDevice(_hVel,_dVel,4*_numOriginalVertices*sizeof(float));

    Vec3d maxV(0,0,0);
    Vec3d maxdV(0,0,0);
    vid=0;
    
    for(RodsIterator rItr=rods.begin(); rItr!=rods.end(); ++rItr)
    {
        Rod *rod = *rItr;
        for(uint j=0; j<rod->nv(); ++j)
        {
            if(!rod->vertFixed(j))
            {
                Vec3d& v = rod->getVelocities()[j];
                for(uint i=0; i<3; ++i)
                {
                    double oldv = v[i];
                    v[i] = _hVel[vid*4+i];
                    maxV[i] = std::max(maxV[i], abs(v[i]));
                    maxdV[i] = std::max(maxdV[i],abs(v[i]-oldv));
                }
            }
            vid++;
        }
    }

    if(printdbg){
        std::cout<<"MaxV = "<<maxV<<std::endl;
        std::cout<<"MaxdV = "<<maxdV<<std::endl;
    }
                
                
    gettimeofday(&t2,NULL);
    elapsedTime = ( ( t2.tv_sec * 1000000 + t2.tv_usec ) - ( t1.tv_sec * 1000000 + t1.tv_usec) ) / 1000.0;
    if(printdbg)
        std::cout<<"Total interpolation time: "<<elapsedTime<<std::endl;
    
    // free dynamic arrays
    //
    if(printdbg)
        std::cout<<"freeing dynamic arrays"<<std::endl;
    unbindAndFreePosTex();
    unbindAndFreeVelTex();
    unbindAndFreeCellStartTex();
    unbindAndFreeCellEndTex();

    delete dCellWeights;
    for(uint i=0; i<3;i++)
    {
        delete dfaceWeights[i];
        delete dinitialVelocities[i];
        delete dfinalVelocities[i];
        delete dCellVelocities[i];

        delete hinitialVelocities[i];
        delete hpsiN[i];
    }
    
    if(printdbg)
        checkDeviceMemory();


}


void VolumetricCollisionsGPU::draw(bool displayGrid, Real displayGridVelocitiesMultiplier, 
                                Real maxDisplayDensity, bool displayPsiN, bool displayPsiD)
{
}
