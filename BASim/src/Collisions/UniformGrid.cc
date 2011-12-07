// UniformGrid.cc
//

#include "UniformGrid.hh"
//#include "BridsonUtil.hh"
//#include <fstream>
//
//namespace BASim {
//
//UniformGrid::UniformGrid()
//{
//    Vec3i dims(1,1,1);
//    Vec3d xmin(0,0,0), xmax(1,1,1);
//    set(dims, xmin, xmax);
//}
//
//UniformGrid::~UniformGrid()
//{
//    clear();
//}
//    
//void UniformGrid::getProximities(ElasticRods &rods, Collisions &collisions)
//{
//    Vec3d rodMin, rodMax;
//    for (uint i=0; i<3; ++i)
//    {
//        rodMin[i] =  std::numeric_limits<Real>::max();
//        rodMax[i] = -std::numeric_limits<Real>::max();
//    }
//
//    for (ElasticRodsIterator rItr=rods.begin(); rItr!=rods.end(); ++rItr)
//    {
//        ElasticRod *rod = *rItr;
//        for (int j=0; j<rod->nv(); ++j)
//        {
//            Vec3d x = rod->getVertex(j);
//            for (uint i=0; i<3; ++i)
//            {
//                // This is assuming the rod has a circular cross section
//                rodMin[i] = std::min(rodMin[i], x[i] - rod->radius());
//                rodMax[i] = std::max(rodMax[i], x[i] + rod->radius());
//            }
//        }
//    }
//
//    Real dx = std::max(rods.front()->getEdgeLength(0), 0.5);
//
//    Vec3i dims(1,1,1);
//    if ((rodMax-rodMin).norm() > 1e-6)
//    {
//        for (uint i=0; i<3; ++i)
//        {
//            uint d = (unsigned int)(std::ceil(rodMax[i] - rodMin[i]) / dx);
//            if (d < 1)
//                d = 1;
//
//            dims[i] = d;
//        }
//    }
//
//    set(dims, rodMin, rodMax);
//
//    EdgeIndexMap edgeIndexMap;
//    int currTotal = 0;
//    for (ElasticRodsIterator rItr=rods.begin(); rItr!=rods.end(); ++rItr)
//    {
//        ElasticRod *rod = *rItr;
//
//        for (int i=0; i<rod->ne(); ++i)
//            edgeIndexMap.insert(std::make_pair(currTotal + i, std::make_pair(currTotal, rod)));
//
//        currTotal += rod->ne();
//    }
//
//    currTotal = 0;
//    std::vector<uint> cands;
//    for (ElasticRodsIterator rItr=rods.begin(); rItr!=rods.end(); ++rItr)
//    {
//        ElasticRod *rod = *rItr;
//
//        std::vector<Vec3d> edgeMins(rod->ne()), edgeMaxs(rod->ne());
//        for (int i=0; i<rod->ne(); ++i)
//        {
//            if (!rod->vertFixed(i) && !rod->vertFixed((i+1)%rod->nv()))
//            {
//                minmax(*(Vec3d *)&rod->getVertex(i),
//                       *(Vec3d *)&rod->getVertex(((i+1)%rod->nv())), edgeMins[i], edgeMaxs[i]);
//                addScalertoVec3d(edgeMins[i], -rod->radius());
//                addScalertoVec3d(edgeMaxs[i], rod->radius());
//
//                cands.clear();
//                findOverlappingElements(edgeMins[i], edgeMaxs[i], cands);
//                for (std::vector<uint>::iterator uiItr=cands.begin(); uiItr!=cands.end(); ++uiItr)
//                {
//                    EdgeIndexMapIterator eimItr = edgeIndexMap.find(*uiItr);
//                    assert(eimItr != edgeIndexMap.end());
//
//                    ElasticRod *rod2 = eimItr->second.second;
//                    int idx2 = *uiItr - eimItr->second.first;
//                    CandidateCollision cand(rod, i, rod2, idx2, EDGE_EDGE);
//                    Collision collision;
//                    if (cand.getProximity(collision))
//                        collisions.push_back(collision);
//                }
//            }
//        }
//
//        // Go backwards so the grid isn't constantly resized
//        //
//        for (uint i=rod->ne(); i>0; --i)
//            addElement(i-1 + currTotal, edgeMins[i-1], edgeMaxs[i-1]);
//
//        currTotal += rod->ne();
//    }
//}
//
//void UniformGrid::getContinuousTimeCollisions(ElasticRods &rods, Real dt, Collisions &collisions)
//{
//    Vec3d rodMin, rodMax;
//    for (uint i=0; i<3; ++i)
//    {
//        rodMin[i] =  std::numeric_limits<Real>::max();
//        rodMax[i] = -std::numeric_limits<Real>::max();
//    }
//
//    for (ElasticRodsIterator rItr=rods.begin(); rItr!=rods.end(); ++rItr)
//    {
//        ElasticRod *rod = *rItr;
//        for (int j=0; j<rod->nv(); ++j)
//        {
//            Vec3d x0 = rod->getStartPositions()[j];
//            Vec3d x1 = rod->getEndPositions()[j];
//            for (uint i=0; i<3; ++i)
//            {
//                rodMin[i] = std::min(rodMin[i], std::min(x0[i], x1[i]));
//                rodMax[i] = std::max(rodMax[i], std::max(x0[i], x1[i]));
//            }
//        }
//    }
//
//    Real dx = std::max(rods.front()->getEdgeLength(0), 0.5);
//
//    Vec3i dims(1,1,1);
//    if ((rodMax-rodMin).norm() > 1e-6)
//    {
//        for (int i=0; i<3; ++i)
//        {
//            int d = (unsigned int)(std::ceil(rodMax[i] - rodMin[i]) / dx);
//            if (d < 1)
//                d = 1;
//
//            dims[i] = d;
//        }
//    }
//
//    UniformGrid grid;
//    grid.set(dims, rodMin, rodMax);
//
//    // TODO: This is a big mess. The UniformGrid stores indices to primitives in
//    // the cells, but with multiple rods we have multiple edges with the same
//    // index, so we need some way to differentiate. So we append the index counts,
//    // so the first edge of the 2nd rod has index (# edges in rod1) + 1. We still
//    // need a way to retrieve the rod that the edge belongs to, so we also store
//    // that in this structure
//    //
//    // A better solution would be to store <Rod *, uint> pairs in the grid cells
//    // It would be cleaner, not sure of the impact on speed
//    //
//    EdgeIndexMap edgeIndexMap;
//    int currTotal = 0;
//    for (ElasticRodsIterator rItr=rods.begin(); rItr!=rods.end(); ++rItr)
//    {
//        ElasticRod *rod = *rItr;
//
//        for (int i=0; i<rod->ne(); ++i)
//            edgeIndexMap.insert(std::make_pair(currTotal + i, std::make_pair(currTotal, rod)));
//
//        currTotal += rod->ne();
//    }
//
//    currTotal = 0;
//    std::vector<uint> cands;
//    for (ElasticRodsIterator rItr=rods.begin(); rItr!=rods.end(); ++rItr)
//    {
//        ElasticRod *rod = *rItr;
//
//        std::vector<Vec3d> edgeMins(rod->ne()), edgeMaxs(rod->ne());
//        for (int i=0; i<rod->ne(); ++i)
//        {
//            if (!rod->vertFixed(i) && !rod->vertFixed((i+1)%rod->nv()))
//            {
//                minmax(rod->getStartPositions()[i],
//                       rod->getEndPositions()[i],
//                       rod->getStartPositions()[(i+1)%rod->nv()],
//                       rod->getEndPositions()[(i+1)%rod->nv()],
//                       edgeMins[i], edgeMaxs[i]);
//                addScalertoVec3d(edgeMins[i], -1e-6);
//                addScalertoVec3d(edgeMaxs[i], 1e-6);
//
//                cands.clear();
//                grid.findOverlappingElements(edgeMins[i], edgeMaxs[i], cands);
//                for (std::vector<uint>::iterator uiItr=cands.begin(); uiItr!=cands.end(); ++uiItr)
//                {
//                    EdgeIndexMapIterator eimItr = edgeIndexMap.find(*uiItr);
//                    assert(eimItr != edgeIndexMap.end());
//
//                    ElasticRod *rod2 = eimItr->second.second;
//                    CandidateCollision candColl(rod, i, rod2, *uiItr - eimItr->second.first, EDGE_EDGE);
//                    candColl.getContinuousTime(dt, collisions);
//                }
//            }
//        }
//
//        // Go backwards so the grid isn't constantly resized
//        //
//        for (uint i=rod->ne(); i>0; --i)
//            grid.addElement(i-1 + currTotal, edgeMins[i-1], edgeMaxs[i-1]);
//
//        currTotal += rod->ne();
//    }
//}
//
//void UniformGrid::boundsToIndices(Vec3d& xmin, Vec3d& xmax, Vec3i& xmini, Vec3i& xmaxi)
//{
//    xmini = floor( divideVec3dByVec3d( (xmin - gridxmin), cellsize) );
//    xmaxi = floor( divideVec3dByVec3d( (xmax - gridxmin), cellsize) );
//
//    if (xmini[0] < 0) xmini[0] = 0;
//    if (xmini[1] < 0) xmini[1] = 0;
//    if (xmini[2] < 0) xmini[2] = 0;
//
//    if (xmaxi[0] >= cells.nx) xmaxi[0] = cells.nx-1;
//    if (xmaxi[1] >= cells.ny) xmaxi[1] = cells.ny-1;
//    if (xmaxi[2] >= cells.nz) xmaxi[2] = cells.nz-1;
//}
//
//void UniformGrid::addElement(unsigned int idx, Vec3d& xmin, Vec3d& xmax)
//{
//    if (elementquery.size() <= idx)
//    {
////        elementidxs.resize(idx+1);
//        elementxmins.resize(idx+1);
//        elementxmaxs.resize(idx+1);
//        elementquery.resize(idx+1);
//    }
//
//    elementxmins[idx] = xmin;
//    elementxmaxs[idx] = xmax;
//    elementquery[idx] = 0;
//
//    Vec3i xmini, xmaxi;
//    boundsToIndices(xmin, xmax, xmini, xmaxi);
//
//    for (int i = xmini[0]; i <= xmaxi[0]; i++)
//        for (int j = xmini[1]; j <= xmaxi[1]; j++)
//            for (int k = xmini[2]; k <= xmaxi[2]; k++)
//            {
//                std::vector<unsigned int>*& cell = cells(i, j, k);
//                if(!cell)
//                    cell = new std::vector<unsigned int>();
//
//                cell->push_back(idx);
////                elementidxs[idx].push_back(Vec3ui(i, j, k));
//            }
//}
//
///*
//void UniformGrid::removeElement(unsigned int idx)
//{
//    for(unsigned int c = 0; c < elementidxs[idx].size(); c++)
//    {
//        Vec3ui cellcoords = elementidxs[idx][c];
//        std::vector<unsigned int>* cell = cells(cellcoords[0], cellcoords[1], cellcoords[2]);
//
//        std::vector<unsigned int>::iterator it = cell->begin();
//        while(*it != idx)
//            it++;
//        cell->erase(it);
//    }
//
//    elementidxs[idx].clear();
//}
//*/
//
///*
//void UniformGrid::updateElement(unsigned int idx, Vec3d& xmin, Vec3d& xmax)
//{
//    removeElement(idx);
//    addElement(idx, xmin, xmax);
//}
//*/
//
//void UniformGrid::set(Vec3i& dims, Vec3d& xmin, Vec3d& xmax)
//{
//    clear();
//
//    gridxmin = xmin;
//    gridxmax = xmax;
//    for (unsigned int i = 0; i < 3; i++)
//        cellsize[i] = (gridxmax[i] - gridxmin[i]) / dims[i];
//
//    cells.resize(dims[0], dims[1], dims[2]);    
//    for (unsigned int i = 0; i < cells.a.size(); i++)
//        cells[i] = 0;
//}
//
//void UniformGrid::clear()
//{
//    for (unsigned int i = 0; i < cells.a.size(); i++)
//    {
//        std::vector<unsigned int>*& cell = cells[i];
//
//        if (cell)
//        {
//            delete cell;
//            cell = 0;
//        }
//    }
//
////    elementidxs.clear();
//    elementxmins.clear();
//    elementxmaxs.clear();
//    elementquery.clear();
//    lastquery = 0;
//}
//
//void UniformGrid::findOverlappingElements(Vec3d& xmin, Vec3d& xmax, std::vector<uint> &idxs)
//{
//    if (lastquery == std::numeric_limits<unsigned int>::max())
//    {
//        for (unsigned int i = 0; i < elementquery.size(); i++)
//            elementquery[i] = 0;
//
//        lastquery = 0;
//    }
//    lastquery++;
//
//    Vec3i xmini, xmaxi;
//    boundsToIndices(xmin, xmax, xmini, xmaxi);
//
//    for (int i = xmini[0]; i <= xmaxi[0]; i++)
//        for (int j = xmini[1]; j <= xmaxi[1]; j++)
//            for (int k = xmini[2]; k <= xmaxi[2]; k++)
//            {
//                std::vector<unsigned int>* cell = cells(i, j, k);
//
//                if (cell)
//                {
//                    for(unsigned int c = 0; c < cell->size(); c++)
//                    {
//                        unsigned int oidx = cell->at(c);
//                        if(elementquery[oidx] < lastquery)
//                        {
//                            elementquery[oidx] = lastquery;
//
//                            Vec3d& oxmin = elementxmins[oidx];
//                            Vec3d& oxmax = elementxmaxs[oidx];
//                            if( isLessThanOrEqual(xmin, oxmax) && isGreaterThanOrEqual(xmax, oxmin))
//                                idxs.push_back(oidx);
//                        }
//                    }
//                }
//            }
//}
//
//void UniformGrid::writeFile(std::string filename)
//{
//    std::fstream gridFile(filename.c_str(), std::ios::out | std::ios::binary);
//
//    gridFile.write(reinterpret_cast<char *>(gridxmin.data()), sizeof(double) * 3);
//    gridFile.write(reinterpret_cast<char *>(gridxmax.data()), sizeof(double) * 3);
//    gridFile.write(reinterpret_cast<char *>(cellsize.data()), sizeof(double) * 3);
//
//    size_t nbrElements = elementquery.size();
//    gridFile.write(reinterpret_cast<char *>(&nbrElements), sizeof(size_t));
//
//    if (nbrElements != 0)
//    {
//        gridFile.write(reinterpret_cast<char *>(elementxmins[0].data()), sizeof(double) * 3 * nbrElements);
//        gridFile.write(reinterpret_cast<char *>(elementxmaxs[0].data()), sizeof(double) * 3 * nbrElements);
//    }
//
////    for (size_t i=0; i<nbrElements; ++i)
////    {
////        size_t nbrElementCells = elementidxs[i].size();
////        gridFile.write(reinterpret_cast<char *>(&nbrElementCells), sizeof(size_t));
////        if (nbrElementCells != 0)
////            gridFile.write(reinterpret_cast<char *>(elementidxs[i][0].data()), sizeof(uint) * 3 * nbrElementCells);
////    }
//
//    gridFile.write(reinterpret_cast<char *>(&cells.nx), sizeof(int));
//    gridFile.write(reinterpret_cast<char *>(&cells.ny), sizeof(int));
//    gridFile.write(reinterpret_cast<char *>(&cells.nz), sizeof(int));
//
//    for (unsigned int i = 0; i < cells.a.size(); i++)
//    {
//        std::vector<unsigned int>*& cell = cells[i];
//        if (cell)
//        {
//            size_t nbrElementsInCell = cell->size();
//            gridFile.write(reinterpret_cast<char *>(&nbrElementsInCell), sizeof(size_t));
//            if (nbrElementsInCell != 0)
//                gridFile.write(reinterpret_cast<char *>(&(*cell)[0]), sizeof(uint) * nbrElementsInCell);
//        }
//        else
//        {
//            size_t nbrElementsInCell = 0;
//            gridFile.write(reinterpret_cast<char *>(&nbrElementsInCell), sizeof(size_t));
//        }
//    }
//}
//
//bool UniformGrid::readFile(std::string filename)
//{
//    std::fstream gridFile(filename.c_str(), std::ios::in | std::ios::binary);
//
//    if (gridFile.fail())
//    {
//        gridFile.clear(std::ios::failbit);
//        return false;
//    }
//
//    gridFile.read(reinterpret_cast<char *>(gridxmin.data()), sizeof(double) * 3);
//    gridFile.read(reinterpret_cast<char *>(gridxmax.data()), sizeof(double) * 3);
//    gridFile.read(reinterpret_cast<char *>(cellsize.data()), sizeof(double) * 3);
//
//    size_t nbrElements;
//    gridFile.read(reinterpret_cast<char *>(&nbrElements), sizeof(size_t));
//
////    elementidxs.resize(nbrElements);
//    elementxmins.resize(nbrElements);
//    elementxmaxs.resize(nbrElements);
//    elementquery.resize(nbrElements);
//   
//    if (nbrElements != 0)
//    {
//        gridFile.read(reinterpret_cast<char *>(elementxmins[0].data()), sizeof(double) * 3 * nbrElements);
//        gridFile.read(reinterpret_cast<char *>(elementxmaxs[0].data()), sizeof(double) * 3 * nbrElements);
//    }
//
//    std::fill(elementquery.begin(), elementquery.end(), 0);
//
////    for (size_t i=0; i<nbrElements; ++i)
////    {
////        size_t nbrElementCells;
////        gridFile.read(reinterpret_cast<char *>(&nbrElementCells), sizeof(size_t));
////        
////        elementidxs[i].resize(nbrElementCells);
////        if (nbrElementCells != 0)
////            gridFile.read(reinterpret_cast<char *>(elementidxs[i][0].data()), sizeof(uint) * 3 * nbrElementCells);
////
////        elementquery[i] = 0;
////    }
//   
//    int nx, ny, nz;
//    gridFile.read(reinterpret_cast<char *>(&nx), sizeof(int));
//    gridFile.read(reinterpret_cast<char *>(&ny), sizeof(int));
//    gridFile.read(reinterpret_cast<char *>(&nz), sizeof(int));
//   
//    cells.resize(nx, ny, nz);
//
//    for (unsigned int i = 0; i < cells.a.size(); i++)
//    {
//        std::vector<unsigned int>*& cell = cells[i];
//        size_t nbrElementsInCell;
//        gridFile.read(reinterpret_cast<char *>(&nbrElementsInCell), sizeof(size_t));
//        if (nbrElementsInCell != 0)
//        {
//            cell = new std::vector<unsigned int>();
//            cell->resize(nbrElementsInCell);
//            gridFile.read(reinterpret_cast<char *>(&(*cell)[0]), sizeof(uint) * nbrElementsInCell);
//        }
//        else
//            cells[i] = 0;
//    }
//
//    lastquery = 0;
//
//    return true;
//}
//
//}
