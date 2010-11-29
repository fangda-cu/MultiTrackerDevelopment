#include "VolumetricCollisionsCPU.hh"

#include <sys/time.h>

#include <Arrays/ARRAYS_3D.h>
#include <Arrays/LIST_ARRAY.h>
#include <Collisions_And_Interactions/SEGMENT_HIERARCHY.h>
#include <Geometry/BOX.h>
#include <Geometry/SEGMENTED_CURVE.h>
#include <Geometry/SEGMENT_3D.h>
#include <Grids/GRID_3D.h>
#include <Grids/SEGMENT_MESH.h>
#include <Grid_Based_Fields/FACE_ARRAYS_3D.h>
#include <Incompressible_Flows/PROJECTION_UNIFORM.h>
#include <Interpolation/LINEAR_INTERPOLATION_UNIFORM.h>
#include <OpenGL/OPENGL_CONSTANT_COLOR_MAP.h>
#include <OpenGL/OPENGL_FACE_SCALAR_FIELD_3D.h>
#include <OpenGL/OPENGL_MAC_VELOCITY_FIELD_3D.h>
#include <OpenGL/OPENGL_GRID_3D.h>
#include <OpenGL/OPENGL_SCALAR_FIELD_3D.h>
#include <OpenGL/OPENGL_COLOR_RAMP.h>
#include <Particles/SOLIDS_PARTICLE.h>
#include <Matrices_And_Vectors/VECTOR.h>

VolumetricCollisionsCPU::VolumetricCollisionsCPU()
    :_initialised(false), _hairCurve(NULL),_projection(NULL), m_flip( 0.95 ), m_slip( 0.1 )
{
    m_stepNumber = 0; 
}

VolumetricCollisionsCPU::VolumetricCollisionsCPU( RodDataMap& i_rodDataMap )
    :_initialised(false), _hairCurve(NULL),_projection(NULL), m_flip( 0.95 ), m_slip( 0.1 )
{
    reinitialiseRods( i_rodDataMap );
    cerr << "VolumetricCollisionsCPU reinitialiseRods. \n" << endl;
}

VolumetricCollisionsCPU::~VolumetricCollisionsCPU()
{
    deleteRods();
}

ElasticRod* VolumetricCollisionsCPU::initialiseRodMap( RodDataMap& i_rodDataMap )
{
    m_rodDataMap = &i_rodDataMap;
    m_rodDataMapIterator = m_rodDataMap->begin();

    WmFigRodGroup* pRodGroup = m_rodDataMapIterator->second;
    m_rodIndex = 0;

    // This assumes every group has at least one rod...
    return pRodGroup->elasticRod( m_rodIndex );
}

ElasticRod* VolumetricCollisionsCPU::nextRod()
{
    // This will add in placeholder rods, need to remove those...

    if ( m_rodDataMapIterator == m_rodDataMap->end() )
    {
        return NULL;
    }

    WmFigRodGroup* pRodGroup = m_rodDataMapIterator->second;

    if ( m_rodIndex <  (pRodGroup->numberOfRods() - 1 ) )
    {
        ++m_rodIndex;        
    }
    else // End of this rod group, move onto the next
    {
        ++m_rodDataMapIterator;
        if ( m_rodDataMapIterator == m_rodDataMap->end() )
        {
            return NULL;
        }
        // otherwise we still have another group to go so reset the index
        m_rodIndex = 0;
    }

    // This assumes every group has at least one rod...
    return pRodGroup->elasticRod( m_rodIndex );
}

void VolumetricCollisionsCPU::reinitialiseRods( RodDataMap& i_rodDataMap )
{
    deleteRods();
    _hairCurve = PhysBAM::SEGMENTED_CURVE<TV>::Create();
    _hairCurve->particles.Store_Velocity();

    int numVertices = 0;
    ElasticRod* rod = initialiseRodMap( i_rodDataMap );    
    while ( rod != NULL )
    {
        numVertices += rod->nv();
        rod = nextRod();
    }

    _hairCurve->particles.Add_Particles(numVertices);

    int id=1;
    rod = initialiseRodMap( i_rodDataMap );
    while ( rod != NULL )
    {
        int firstcv = id;
        
        for (uint j=0; j<rod->nv(); ++j)
        {

            BASim::Vec3d x = rod->getStartPositions()[j];
            BASim::Vec3d v = rod->getVelocities()[j];
            for (uint i=0; i<3; ++i)
            {
                _hairCurve->particles.X(id)[i+1]=x[i];
                _hairCurve->particles.V(id)[i+1]=v[i];
            }

            if(j<rod->ne())
            {
                _hairCurve->mesh.elements.Append(PhysBAM::VECTOR<int,2>(id,id+1));
            }
            id++;
        }
        
        rod = nextRod();
    }
    _projection = new PhysBAM::PROJECTION_UNIFORM<GRID_3D>(_grid);
    _initialised = true;
}

void VolumetricCollisionsCPU::deleteRods()
{
    if(_initialised)
    {
        delete _hairCurve;
        delete _projection;
    }
    _initialised = false;
}

void VolumetricCollisionsCPU::respondVolumetricCollisions( RodDataMap& i_rodDataMap, Real targetEdgeDensity, 
                                                       Real volumetricRadius, Real gridDx,
                                                       BASim::Vec3d separationCondition, double flip,
                                                       double slip,
                                                       CollisionMeshDataHashMap &collisionMeshes)
{
    m_flip = flip;
    m_slip = slip;

    bool printdbg = true;
    double rsq = volumetricRadius*volumetricRadius;
    double oneOverR=1.0/volumetricRadius;
    double oneOverR4 = 1.0/(rsq*rsq);
    timeval t1, t2;
    double elapsedTime;

    gettimeofday(&t1,NULL);
    
    // update positions/velocities in curves (use saved positions)
    //
    int id=0;
    ElasticRod* rod = initialiseRodMap( i_rodDataMap );    
    while ( rod != NULL )
    {
        for (uint j=0; j<rod->nv(); ++j)
        {
            id++;
            BASim::Vec3d x = rod->getStartPositions()[j];
            BASim::Vec3d v = rod->getVelocities()[j];
            for (uint i=0; i<3; ++i)
            {
                _hairCurve->particles.X(id)[i+1]=x[i];
                _hairCurve->particles.V(id)[i+1]=v[i];
            }
        }
        rod = nextRod();
    }

    if(id==0)
    {
        std::cout<<"\nNo rods found. exiting volume"<<std::endl;
        return;
    }

    double oneOverDx = 1.0/gridDx;
    
    PhysBAM::BOX<TV> boundingBox = PhysBAM::BOX<TV>::Bounding_Box(_hairCurve->particles.X);
    TV rodMin = boundingBox.Minimum_Corner();
    TV rodMax = boundingBox.Maximum_Corner();

    // Create grid
    //
    PhysBAM::VECTOR<int,3> dims;
    for(int d=1; d<=3; d++)
    {
        rodMin(d) = std::floor((rodMin(d) - volumetricRadius)* oneOverDx - 1.0);
        rodMax(d) = std::floor((rodMax(d) + volumetricRadius) * oneOverDx + 1.0);
        dims(d) = (int)(rodMax(d) - rodMin(d));
        rodMin(d) *= gridDx;
        rodMax(d) *= gridDx;
    }

    boundingBox=PhysBAM::BOX<TV>(rodMin,rodMax);
    _grid.Initialize(dims,boundingBox,true);

    if(printdbg)
        std::cout<<"\nInitialized "<<dims(1)<<" by "<<dims(2)<<" by "<<dims(3)<<" grid"<<std::endl;

    if(dims.Max()>100)
        return;


    _projection->elliptic_solver->pcg.Show_Results(printdbg);
    _projection->elliptic_solver->pcg.cg_restart_iterations = 0;
    _projection->elliptic_solver->Set_Relative_Tolerance((double)1e-5);
    _projection->elliptic_solver->Set_Absolute_Tolerance((double)1e-7);
    _projection->elliptic_solver->pcg.Set_Maximum_Iterations(500);

    // initialize everything on the grid
    //
    _projection->Initialize_Grid(_grid);
    _faceWeights.Resize(_grid,1,false,false);_faceWeights.Fill(0);
    _cellWeights.Resize(_grid,1,false,false);_cellWeights.Fill(0);
    _cellVelocities.Resize(_grid,1,false,false);_cellVelocities.Fill(TV());
    _separatingCells.Resize(_grid,1,false,false);_separatingCells.Fill(false);
    _initialVelocities.Resize(_grid,1,false,false);_initialVelocities.Fill(0);
    _projection->face_velocities = _initialVelocities;

    // set up projection for divergence control
    //
    _projection->Use_Non_Zero_Divergence(true);
    _projection->Use_Divergence_Multiplier(true);
    _projection->divergence_multiplier.Fill(0);
    _projection->divergence.Fill(0);
    _projection->elliptic_solver->psi_N.Fill(false);
    _projection->elliptic_solver->psi_D.Fill(false);

    gettimeofday(&t2, NULL);
    elapsedTime = ( ( t2.tv_sec * 1000000 + t2.tv_usec ) - ( t1.tv_sec * 1000000 + t1.tv_usec) ) / 1000.0;
    if(printdbg)
        std::cout<<"Total initialization time: "<<elapsedTime<<" ms"<<std::endl;


    gettimeofday(&t1, NULL);

    if(printdbg)
        std::cout << "initializing segment hierarchy... \n";
    _hairCurve->Initialize_Hierarchy();
    
    PhysBAM::SEGMENT_MESH& guideMesh = _hairCurve->Get_Segment_Mesh();
    PhysBAM::SEGMENT_HIERARCHY<TV>& guideHierarchy = *(_hairCurve->hierarchy);

    // Rasterize face weights and velocities
    //
    if(printdbg)
        std::cout << "rasterizing face data... \n";

    PhysBAM::LIST_ARRAY<int> intersectionList;
    for(pbFaceIterator faceItr(_grid,1); faceItr.Valid(); faceItr.Next())
    {
        PhysBAM::VECTOR<int,3> face = faceItr.Face_Index();
        int axis = faceItr.Axis();
        TV faceCenter = _grid.Face(axis,face);
        
        BASim::Vec3d x,v;
        for(int d=0; d<3; d++)
            x[d] = faceCenter[d+1];

        bool insideCollisionObject = false;
        CollisionMeshData *cmData;
        for (CollisionMeshDataHashMapIterator cmItr = collisionMeshes.begin();
             cmItr != collisionMeshes.end(); cmItr++)
        {
            cmData = cmItr->second;
            if(cmData->getLevelSetValue(x,v) <= 0.0)
            {
                insideCollisionObject = true;
                break;
            }
        }
        if(insideCollisionObject)
        {
            _projection->elliptic_solver->psi_N(axis,face)=true;
            _projection->face_velocities(axis,face) = v[axis-1];
            _initialVelocities(axis,face)=_projection->face_velocities(axis,face);
        }else
        {
        
            intersectionList.Remove_All();
            guideHierarchy.Intersection_List(faceCenter,intersectionList,volumetricRadius);
            
            for(int k=1; k<=intersectionList.m; k++)
            {
                int id1 = guideMesh.elements(intersectionList(k))[1];
                int id2 = guideMesh.elements(intersectionList(k))[2];
                PhysBAM::SEGMENT_3D<double> segment(_hairCurve->particles.X(id1), _hairCurve->particles.X(id2));
                PhysBAM::VECTOR<double,2> bweights = segment.Clamped_Barycentric_Coordinates(faceCenter);
                TV closestPoint = segment.x1 + bweights(2)*(segment.x2-segment.x1);
                double dsq = (closestPoint-faceCenter).Magnitude_Squared();
                
                if(dsq < rsq)
                {
                    double dist = sqrt(dsq);
//                  double weight = 1-dist*oneOverR;
                    double weight = oneOverR4*(volumetricRadius-dist);
                    _faceWeights(axis,face) += weight;
                    _initialVelocities(axis,face) += weight *
                        (bweights(1)*_hairCurve->particles.V(id1)[axis] + 
                         bweights(2)*_hairCurve->particles.V(id2)[axis]);
                }
            }
            if(_faceWeights(axis,face))
                _initialVelocities(axis,face) /= _faceWeights(axis,face);
            _projection->face_velocities(axis,face) = _initialVelocities(axis,face);
        }
    }

    // Rasterize cell weights and 
    // set up projection density control
    //
    double targetWeightPerCell = targetEdgeDensity;
    double maxCellWeight = 0.0;
    std::cout << "TargetWeightPerCell = " << targetWeightPerCell << std::endl;

    if(printdbg)
        std::cout << "rasterizing cell data... ";

    gettimeofday(&t1, NULL);

    for(pbCellIterator cellItr(_grid,1); cellItr.Valid(); cellItr.Next())
    {
        PhysBAM::VECTOR<int,3> cell = cellItr.Cell_Index();
        TV cellCenter = _grid.Center(cell);

        BASim::Vec3d x,v;
        for(int d=0; d<3; d++)
            x[d] = cellCenter[d+1];

        bool insideCollisionObject = false;
        CollisionMeshData *cmData;
        for (CollisionMeshDataHashMapIterator cmItr = collisionMeshes.begin();
             cmItr != collisionMeshes.end(); cmItr++)
        {
            cmData = cmItr->second;
            if(cmData->getLevelSetValue(x,v) <= 0.0)
            {
                insideCollisionObject = true;
                break;
            }
        }

        if(insideCollisionObject)
        {
            _cellWeights(cell)=targetWeightPerCell;
            if (_projection->divergence.Valid_Index(cell))
            {
                _projection->divergence_multiplier(cell)=1.0;
                _projection->divergence(cell)=0.0;
            }
            for(int d=0; d<3; d++)
                _cellVelocities(cell)[d+1]=v[d];
        }
        else
        {
            
            intersectionList.Remove_All();
            guideHierarchy.Intersection_List(cellCenter,intersectionList,volumetricRadius);
            
            for(int k=1; k<=intersectionList.m; k++)
            {
                int id1 = guideMesh.elements(intersectionList(k))[1];
                int id2 = guideMesh.elements(intersectionList(k))[2];
                PhysBAM::SEGMENT_3D<double> segment(_hairCurve->particles.X(id1), _hairCurve->particles.X(id2));
                PhysBAM::VECTOR<double,2> bweights = segment.Clamped_Barycentric_Coordinates(cellCenter);
                TV closestPoint = segment.x1 + bweights(2)*(segment.x2-segment.x1);
                double dsq = (closestPoint-cellCenter).Magnitude_Squared();
                
                if(dsq < rsq)
                {
                    double dist = sqrt(dsq);
//                  double weight = 1-dist*oneOverR;
                    double weight = oneOverR4*(volumetricRadius-dist);
                    _cellWeights(cell) += weight;
                    _cellVelocities(cell) += weight*
                        (bweights(1)*_hairCurve->particles.V(id1) + 
                         bweights(2)*_hairCurve->particles.V(id2));
                }
            }
            maxCellWeight = max(maxCellWeight, _cellWeights(cell));
            
            if(_projection->divergence.Valid_Index(cell))
                if(_cellWeights(cell))
                {
                    _cellVelocities(cell) /= _cellWeights(cell);
                    double normalizedWeight = _cellWeights(cell) / targetWeightPerCell;
                    _projection->divergence_multiplier(cell) = min(1.0, normalizedWeight);
                    _projection->divergence(cell) = max(0.0, normalizedWeight-1.0);
                } else {
                    _projection->elliptic_solver->psi_D(cell) = true;
                    _projection->elliptic_solver->u(cell) = 0;
                    _projection->divergence_multiplier(cell) = 1.0;
                    _projection->divergence(cell) = 0.0;
                }
        }
    }

    
    std::cout<<"Max cell weight = "<<maxCellWeight<<std::endl;

    gettimeofday(&t2, NULL);
    elapsedTime = ( ( t2.tv_sec * 1000000 + t2.tv_usec ) - ( t1.tv_sec * 1000000 + t1.tv_usec) ) / 1000.0;
    if(printdbg)
        std::cout<<"Total rasterization time: "<<elapsedTime<<" ms"<<std::endl;


    // Apply separation condition
    //
    vector<int> separationAxis; //vector of axes with a separation condition
    for(int a=0; a<3; a++)
        if(separationCondition[a]>=0.0)
            separationAxis.push_back(a+1);
    int separationsTriggered = 0;
    if(printdbg)
        std::cout<<"SeparationCondition: "<<separationCondition<<std::endl;

    if(separationAxis.size())
    {
        if(printdbg)
            std::cout<<"applying separation condition...\n";

        for(pbCellIterator cellItr(_grid,1); cellItr.Valid(); cellItr.Next())
        {
            PhysBAM::VECTOR<int,3> cell = cellItr.Cell_Index();
            for(int a=0; a<separationAxis.size(); a++) //compare with right, above, front
            {
                int axis = separationAxis[a];
                PhysBAM::VECTOR<int,3> face = cellItr.Second_Face_Index(axis);
                PhysBAM::VECTOR<int,3> neighbor = cellItr.Cell_Neighbor(2*axis);

                // if separating face is inside collision object, 
                // only apply separation condition to cell which is separating 
                // from collision object
                //
                if(_projection->elliptic_solver->psi_N(axis,face))
                {
                    if (_projection->divergence.Valid_Index(neighbor) && 
                        !_projection->elliptic_solver->psi_D(neighbor) && 
                        _cellVelocities(neighbor)[axis] - _initialVelocities(axis,face) >separationCondition[axis-1] )
                    {
                        _projection->elliptic_solver->psi_D(neighbor)=true;
                        _projection->elliptic_solver->u(neighbor)=0;
                        _projection->divergence_multiplier(neighbor) = 1.0;
                        _projection->divergence(neighbor)=0.0;
                        _separatingCells(neighbor)=true;
                        separationsTriggered++;
                    }
                    if (!_projection->elliptic_solver->psi_D(cell) 
                        && _initialVelocities(axis,face) - _cellVelocities(cell)[axis] > separationCondition[axis -1])
                    {
                        _projection->elliptic_solver->psi_D(cell)=true;
                        _projection->elliptic_solver->u(cell)=0;
                        _projection->divergence_multiplier(cell) = 1.0;
                        _projection->divergence(cell)=0.0;
                        _separatingCells(cell)=true;
                        separationsTriggered++;
                    }
                    
                } 
                
                // if cell velocities are separating
                //
                else if(_projection->divergence.Valid_Index(neighbor) && 
                        _cellWeights(cell) > 2.0 && _cellWeights(neighbor) > 2.0 &&
                   _cellVelocities(neighbor)[axis] - _cellVelocities(cell)[axis] > separationCondition[axis-1])
                {
                    _projection->elliptic_solver->psi_D(cell)=true;
                    _projection->elliptic_solver->u(cell)=0;
                    _projection->divergence_multiplier(cell) = 1.0;
                    _projection->divergence(cell)=0.0;
                    separationsTriggered++;
                
                    _projection->elliptic_solver->psi_D(neighbor)=true;
                    _projection->elliptic_solver->u(neighbor)=0;
                    _projection->divergence_multiplier(neighbor) = 1.0;
                    _projection->divergence(neighbor)=0.0;
                    separationsTriggered++;
                    
                    _separatingCells(cell)=true;
                    _separatingCells(neighbor)=true;
                }
            }
        }
    }
    if(printdbg)
        std::cout<<"Separations Triggered: "<<separationsTriggered<<std::endl;


    // Solve for divergence free grid velocity
    //
    if(printdbg)
        std::cout<<"making divergence free... \n";

    gettimeofday(&t1,NULL);
    _projection->Make_Divergence_Free(1.0,1.0);
    gettimeofday(&t2,NULL);
    elapsedTime = ( ( t2.tv_sec * 1000000 + t2.tv_usec ) - ( t1.tv_sec * 1000000 + t1.tv_usec) ) / 1000.0;
    if(printdbg)
        std::cout<<"Total projection time: "<<elapsedTime<<" ms"<<std::endl;

    // Compute change in velocity
    //
    if(printdbg)
        std::cout<<"computing deltas ...\n";

    gettimeofday(&t1,NULL);
    _deltas.Resize(_grid,1,false,false);
    for(pbFaceIterator faceItr(_grid,1); faceItr.Valid(); faceItr.Next())
    {
        PhysBAM::VECTOR<int,3> face = faceItr.Face_Index();
        int axis = faceItr.Axis();
        _deltas(axis,face) = _projection->face_velocities(axis,face)-_initialVelocities(axis,face);
    }

    // Interpolate back to particles
    //
    if(printdbg)
        std::cout<<"interpolating back to particles... \n";
    
    PhysBAM::LINEAR_INTERPOLATION_UNIFORM<GRID_3D,double> interpolation;
    id=0;
    BASim::Vec3d maxV(0,0,0);
    BASim::Vec3d startMaxV(0,0,0);
    BASim::Vec3d maxdV(0,0,0);
    double oldVi;

    rod = initialiseRodMap( i_rodDataMap );
    while ( rod != NULL )
    {
        Real slip = m_slip, omslip = 1.0-slip;
        Real flip = m_flip, omflip = 1.0-flip;
        if(slip == 1.0)
            id += rod->nv();
        else
        {
            for (uint j=0; j<rod->nv(); ++j)
            {
                id++;
                if(!rod -> vertFixed(j))
                {
                    TV cell_double = (_hairCurve->particles.X(id) - rodMin) * oneOverDx + 1;
                    PhysBAM::VECTOR<int,3> cell;
                    for(int d=1; d<=3; d++)
                        cell(d) = (int)cell_double(d);
                    
                    // Don't modify velocity if in nearly empty cell
                    // or in cell marked as separating
                    //
                    if(_cellWeights(cell) > 2.0 && !_separatingCells(cell)){

                        BASim::Vec3d& v = rod->getVelocities()[j];
                     //   cerr << "veclocity before (" << j << ") = " << v << endl;
                      
                        TV divFreeV = interpolation.Clamped_To_Array_Face(_grid,
                                                                          _projection->face_velocities,
                                                                          _hairCurve->particles.X(id));
                        TV deltaV = interpolation.Clamped_To_Array_Face(_grid,
                                                                        _deltas,
                                                                        _hairCurve->particles.X(id));
                        
                        for (uint i=0; i<3; ++i)
                        {
                            oldVi=v[i];
                            startMaxV[i] = std::max(startMaxV[i],abs(v[i]));
                            v[i] = omslip * (flip * (deltaV[i+1] + v[i]) + 
                                             omflip * divFreeV[i+1]) + slip * v[i];
                            maxV[i] = std::max(maxV[i],abs(v[i]));
                            maxdV[i]=std::max(maxdV[i],abs(v[i]-oldVi));
                        }
                    //    cerr << "veclocity after (" << j << ") = " << v << endl;
                    }
                }
            }
        }
        rod = nextRod();
    }
    if(printdbg){
        std::cout<<"startMaxV = "<<startMaxV<<std::endl;
        std::cout<<"MaxV = "<<maxV<<std::endl;
        std::cout<<"MaxdV = "<<maxdV<<std::endl;
    }
    gettimeofday(&t2,NULL);
    elapsedTime = ( ( t2.tv_sec * 1000000 + t2.tv_usec ) - ( t1.tv_sec * 1000000 + t1.tv_usec) ) / 1000.0;
    if(printdbg)
        std::cout<<"Total interpolation time: "<<elapsedTime<<" ms"<<std::endl;

    m_stepNumber++;
}


void VolumetricCollisionsCPU::draw(bool displayGrid, Real displayGridVelocitiesMultiplier, 
                                Real maxDisplayDensity, bool displayPsiN, bool displayPsiD)
{
    if(!_initialised)
        return;

    if(displayGrid)
        PhysBAM::OPENGL_GRID_3D<double>(_grid).Display();
    if(displayGridVelocitiesMultiplier && _grid.Is_MAC_Grid())
    {
        PhysBAM::FACE_ARRAYS_3D<double> velocities = _projection->face_velocities;
        velocities*=displayGridVelocitiesMultiplier;

        PhysBAM::OPENGL_MAC_VELOCITY_FIELD_3D<double>
            (_grid, velocities.u, velocities.v, velocities.w).Display();
    }
    
    if(maxDisplayDensity)
        PhysBAM::OPENGL_SCALAR_FIELD_3D<double,double>(_grid, _cellWeights,
                                                       PhysBAM::OPENGL_COLOR_RAMP<double>::
                                                       Matlab_Jet(0, maxDisplayDensity)).Display();

    if(displayPsiN ){
        PhysBAM::FACE_ARRAYS_3D<bool> &psiN = _projection->elliptic_solver->psi_N;
        PhysBAM::OPENGL_CONSTANT_COLOR_MAP<bool> psiNcm(PhysBAM::OPENGL_COLOR::Cyan());
        PhysBAM::OPENGL_FACE_SCALAR_FIELD_3D<double,bool> psiNgl(_grid, psiN.u, psiN.v, psiN.w, &psiNcm);
        psiNgl.Update();
        psiNgl.Display();
    }
    
    if(displayPsiD){
//      PhysBAM::ARRAYS_3D<bool> &psiD = _projection->elliptic_solver->psi_D;
        PhysBAM::ARRAYS_3D<bool> &psiD = _separatingCells;
        PhysBAM::OPENGL_CONSTANT_COLOR_MAP<bool> psiDcm(PhysBAM::OPENGL_COLOR::Magenta());
        PhysBAM::OPENGL_SCALAR_FIELD_3D<double,bool> psiDgl(_grid, psiD, &psiDcm,PhysBAM::OPENGL_SCALAR_FIELD_3D<double,bool>::DRAW_POINTS);
        psiDgl.Update();
        psiDgl.Display();
    }
}
