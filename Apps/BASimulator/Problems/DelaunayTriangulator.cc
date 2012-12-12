/**
 * \file DelaunayTriangulator.cc
 *
 * \author fang@cs.columbia.edu
 * \date Dec 10, 2012
 */

#include "DelaunayTriangulator.hh"
#include "BASim/src/Core/TopologicalObject/TopObjUtil.hh"

namespace DelaunayTriangulator
{

  DelaunayTriangulator::DelaunayTriangulator(Mesh * mesh, VertexProperty<Vec3d> & positions) :
    m_mesh(mesh),
    m_positions(mesh)
  {
    assert(checkDelaunay());
    
    m_positions = positions;
  }

  bool DelaunayTriangulator::insertVertex(Vec3d & v)
  {
    // algorithm in figure 9 of Ledoux paper
    VertexHandle p = m_mesh->addVertex();
    m_vertices.push_back(p);
    m_positions[p] = v;
    
    // find the tet containing v. could be optimized (using an algorithm like Walk in the paper)
    TetHandle th;
    for (TetIterator tit = m_mesh->tets_begin(); tit != m_mesh->tets_end(); ++tit)
    {
      th = *tit;
      TetVertexIterator tvit = m_mesh->tv_iter(th); assert(tvit);
      VertexHandle a = *tvit; ++tvit; assert(tvit);
      VertexHandle b = *tvit; ++tvit; assert(tvit);
      VertexHandle c = *tvit; ++tvit; assert(tvit);
      VertexHandle d = *tvit; ++tvit; assert(!tvit);
      
      if (predicateInTetrahedron(a, b, c, d, p))
        break;
    }
    
    // flip14
    FaceHandle oldf[4];
    TetFaceIterator tfit = m_mesh->tf_iter(th); assert(tfit);
    oldf[0] = *tfit; ++tfit; assert(tfit);
    oldf[1] = *tfit; ++tfit; assert(tfit);
    oldf[2] = *tfit; ++tfit; assert(tfit);
    oldf[3] = *tfit; ++tfit; assert(!tfit);
    
    m_mesh->deleteTet(th, false);
    
    std::set<EdgeHandle> oldeset;
    for (int i = 0; i < 4; i++)
      for (FaceEdgeIterator feit = m_mesh->fe_iter(oldf[i]); feit; ++feit)
        oldeset.insert(*feit);
    
    std::vector<EdgeHandle> olde;
    olde.assign(oldeset.begin(), oldeset.end());
    assert(olde.size() == 6);
    
    FaceHandle newf[6];
    for (int i = 0; i < 6; i++)
      newf[i] = m_mesh->addFace(m_mesh->fromVertex(olde[i]), m_mesh->toVertex(olde[i]), p);

    TetHandle newt[4];
    for (int i = 0; i < 4; i++)
    {
      std::vector<FaceHandle> newfj;
      for (FaceEdgeIterator feit = m_mesh->fe_iter(oldf[i]); feit; ++feit)
      {
        EdgeHandle oldej = *feit;
        int k;
        for (k = 0; k < 6; k++)
          if (olde[k] == oldej)
            break;
        
        assert(k < 6);
        newfj.push_back(newf[k]);
      }
      assert(newfj.size() == 3);
      
      newt[i] = m_mesh->addTet(oldf[i], newfj[0], newfj[1], newfj[2]);
    }
    
    // iterative flipping
    std::vector<TetHandle> stack;
    for (int i = 0; i < 4; i++)
      stack.push_back(newt[i]);
    
    while (stack.size() > 0)
    {
      TetHandle tau = stack.back();
      stack.pop_back();
      assert(tau.isValid());
      
      FaceHandle abc;
      for (TetFaceIterator tfit = m_mesh->tf_iter(tau); tfit; ++tfit)
      {
        abc = *tfit;
        VertexHandle vh;
        getTetFourthVertex(*m_mesh, tau, abc, vh);
        assert(vh.isValid());
        if (vh == p)
          break;
      }
      assert(abc.isValid());
      
      FaceVertexIterator fvit = m_mesh->fv_iter(abc); assert(fvit);
      VertexHandle a = *fvit; ++fvit; assert(fvit);
      VertexHandle b = *fvit; ++fvit; assert(fvit);
      VertexHandle c = *fvit; ++fvit; assert(!fvit);
      
      TetHandle taua = getFaceOtherTet(*m_mesh, abc, tau);
      if (!taua.isValid())
        continue;
      
      VertexHandle d;
      getTetFourthVertex(*m_mesh, taua, abc, d);
      assert(d.isValid());
      
      if (predicateInSphere(a, b, c, p, d) < 0)
      {
        // local non-delaunay, need flpping
        // algorithm in figure 12 of Ledoux paper        
        
        // determine the intersection of line pd with facd abc
        Vec3d nabc = (pos(b) - pos(a)).cross(pos(c) - pos(a));
        Vec3d cabc = (pos(a) + pos(b) + pos(c)) / 3;
        Scalar t = (cabc - pos(p)).dot(nabc) / (pos(d) - pos(p)).dot(nabc);
        assert(t == t);
        
//        bool inside = predicateInTetrahedron(Vec3d(0, 0, 0), pos(a), pos(b), pos(c), pos(p) + (pos(d) - pos(p)) * t);
        Vec3d o = cabc + nabc;
        Vec3d proj = pos(p) + (pos(d) - pos(p)) * t;
        Scalar oabc = predicateOriented(o, pos(a), pos(b), pos(c));
        Scalar opbc = predicateOriented(o, proj, pos(b), pos(c));
        Scalar oapc = predicateOriented(o, pos(a), proj, pos(c));
        Scalar oabp = predicateOriented(o, pos(a), pos(b), proj);
        
        bool inside = (oabc * opbc >= 0 && oabc * oapc >= 0 && oabc * oabp >= 0);
        
        // ignore the degenerate cases
        if (inside)
        {
          // flip23
          TetHandle pabd, pacd, pbcd;
          flip23(tau, taua, abc, pabd, pacd, pbcd);
          stack.push_back(pabd);
          stack.push_back(pacd);
          stack.push_back(pbcd);
        } else
        {
          TetHandle pabd;
          if (oabc * opbc < 0)
            pabd = findTet(*m_mesh, p, b, c, d);
          if (oabc * oapc < 0)
            pabd = findTet(*m_mesh, p, a, c, d);
          if (oabc * oabp < 0)
            pabd = findTet(*m_mesh, p, a, b, d);
          // ignore the cases where two are both negative

          if (pabd.isValid())
          {
            // flip32
            TetHandle pacd, pbcd;
            flip32(tau, taua, pabd, abc, pacd, pbcd);
            stack.push_back(pacd);
            stack.push_back(pbcd);
          } else
          {
            // nothing to do
          }
        }
      }
    }
    
    assert(checkDelaunay());
    
    return true;
  }
      
  void DelaunayTriangulator::flip23(TetHandle pabc, TetHandle abcd, FaceHandle abc, TetHandle & pabd, TetHandle & pacd, TetHandle & pbcd)
  {
    VertexHandle p;
    getTetFourthVertex(*m_mesh, pabc, abc, p);
    assert(p.isValid());
    
    VertexHandle d;
    getTetFourthVertex(*m_mesh, abcd, abc, d);
    assert(p.isValid());
    
    FaceVertexIterator fvit = m_mesh->fv_iter(abc); assert(fvit);
    VertexHandle a = *fvit; ++fvit; assert(fvit);
    VertexHandle b = *fvit; ++fvit; assert(fvit);
    VertexHandle c = *fvit; ++fvit; assert(!fvit);
    
    FaceHandle pab = getVertexOppositeFaceInTet(*m_mesh, pabc, c);
    FaceHandle pbc = getVertexOppositeFaceInTet(*m_mesh, pabc, a);
    FaceHandle pac = getVertexOppositeFaceInTet(*m_mesh, pabc, b);
    
    FaceHandle abd = getVertexOppositeFaceInTet(*m_mesh, abcd, c);
    FaceHandle bcd = getVertexOppositeFaceInTet(*m_mesh, abcd, a);
    FaceHandle acd = getVertexOppositeFaceInTet(*m_mesh, abcd, b);
    
    m_mesh->deleteTet(pabc, false);
    m_mesh->deleteTet(abcd, false);
    
    m_mesh->deleteFace(abc, false);
    
    FaceHandle pda = m_mesh->addFace(p, d, a);
    FaceHandle pdb = m_mesh->addFace(p, d, b);
    FaceHandle pdc = m_mesh->addFace(p, d, c);
    
    pabd = m_mesh->addTet(pda, pdb, pab, abd);
    pbcd = m_mesh->addTet(pdb, pdc, pbc, bcd);
    pacd = m_mesh->addTet(pda, pdc, pac, acd);
  }
  
  void DelaunayTriangulator::flip32(TetHandle pabc, TetHandle abcd, TetHandle pabd, FaceHandle abc, TetHandle & pacd, TetHandle & pbcd)
  {
    VertexHandle p;
    getTetFourthVertex(*m_mesh, pabc, abc, p);
    assert(p.isValid());
    
    VertexHandle d;
    getTetFourthVertex(*m_mesh, abcd, abc, d);
    assert(p.isValid());
    
    FaceVertexIterator fvit = m_mesh->fv_iter(abc); assert(fvit);
    VertexHandle a = *fvit; ++fvit; assert(fvit);
    VertexHandle b = *fvit; ++fvit; assert(fvit);
    VertexHandle c = *fvit; ++fvit; assert(!fvit);
    
    for (int i = 0; i < 3; i++)
    {
      if (tetContainsVertex(*m_mesh, pabd, a) && tetContainsVertex(*m_mesh, pabd, b))
        break;
      VertexHandle t = a;
      a = b;
      b = c;
      c = t;
    }
    assert(tetContainsVertex(*m_mesh, pabd, a) && tetContainsVertex(*m_mesh, pabd, b));
    
    FaceHandle pab = getVertexOppositeFaceInTet(*m_mesh, pabc, c);
    FaceHandle pbc = getVertexOppositeFaceInTet(*m_mesh, pabc, a);
    FaceHandle pac = getVertexOppositeFaceInTet(*m_mesh, pabc, b);
    
    FaceHandle abd = getVertexOppositeFaceInTet(*m_mesh, abcd, c);
    FaceHandle bcd = getVertexOppositeFaceInTet(*m_mesh, abcd, a);
    FaceHandle acd = getVertexOppositeFaceInTet(*m_mesh, abcd, b);
    
    FaceHandle pad = getVertexOppositeFaceInTet(*m_mesh, pabd, b);
    FaceHandle pbd = getVertexOppositeFaceInTet(*m_mesh, pabd, a);
    
    EdgeHandle ab = findEdge(*m_mesh, a, b);
    
    m_mesh->deleteTet(pabc, false);
    m_mesh->deleteTet(abcd, false);
    m_mesh->deleteTet(pabd, false);
    
    m_mesh->deleteFace(abc, false);
    m_mesh->deleteFace(pab, false);
    m_mesh->deleteFace(abd, false);
    
    m_mesh->deleteEdge(ab, false);
    
    FaceHandle pcd = m_mesh->addFace(p, c, d);
    
    pacd = m_mesh->addTet(pad, pcd, pac, acd);
    pbcd = m_mesh->addTet(pbd, pcd, pbc, bcd);
  }
  
  bool DelaunayTriangulator::checkDelaunay()
  {
    for (TetIterator tit = m_mesh->tets_begin(); tit != m_mesh->tets_end(); ++tit)
    {
      TetVertexIterator tvit = m_mesh->tv_iter(*tit); assert(tvit);
      VertexHandle a = *tvit; ++tvit; assert(tvit);
      VertexHandle b = *tvit; ++tvit; assert(tvit);
      VertexHandle c = *tvit; ++tvit; assert(tvit);
      VertexHandle d = *tvit; ++tvit; assert(!tvit);
      
      // could be optimized: no need to check all vertices
      for (VertexIterator vit = m_mesh->vertices_begin(); vit != m_mesh->vertices_end(); ++vit)
      {
        if (*vit == a || *vit == b || *vit == c || *vit == d)
          continue;

        if (predicateInSphere(a, b, c, d, *vit) < 0)
          return false;
      }
    }
  
    return true;
  }

  bool DelaunayTriangulator::extractVoronoiDiagram(Mesh * tomesh, VertexProperty<Vec3d> & vdpos, const Vec3d & bbmin, const Vec3d & bbmax)
  {
    /////////////////////////////////////////////////////////////////////////////////////////
    // Extract an unbounded Voronoi Diagram
    
    // tet in dt = vertex in vd
    // create a vd vertex for each dt tet, located at the circumcenter of the tet
    TetProperty<VertexHandle> dt_tet_2_vd_vert(m_mesh);
    VertexProperty<TetHandle> vd_vert_2_dt_tet(tomesh);    
    VertexProperty<std::vector<VertexHandle> > dt_vert_2_vd_verts(m_mesh);
    for (TetIterator tit = m_mesh->tets_begin(); tit != m_mesh->tets_end(); ++tit)
    {
      TetVertexIterator tvit = m_mesh->tv_iter(*tit); assert(tvit);
      VertexHandle a = *tvit; ++tvit; assert(tvit);
      VertexHandle b = *tvit; ++tvit; assert(tvit);
      VertexHandle c = *tvit; ++tvit; assert(tvit);
      VertexHandle d = *tvit; ++tvit; assert(!tvit);
      
      // find the circumcenter of tet abcd
      //  reference:
      //    http://people.sc.fsu.edu/~jburkardt/presentations/cg_lab_tetrahedrons.pdf
      Vec3d circumcenter;
      Eigen::Matrix<Scalar, 4, 5> mat;
      mat.setZero();
      mat.block<1, 3>(0, 1) = pos(a);
      mat.block<1, 3>(1, 1) = pos(b);
      mat.block<1, 3>(2, 1) = pos(c);
      mat.block<1, 3>(3, 1) = pos(d);
      mat(0, 0) = pos(a).squaredNorm();
      mat(1, 0) = pos(b).squaredNorm();
      mat(2, 0) = pos(c).squaredNorm();
      mat(3, 0) = pos(d).squaredNorm();
      mat.col(4).setOnes();
      
      Scalar alpha = mat.block<4, 4>(0, 1).determinant();
      Mat4d Dx, Dy, Dz;
      Dx.col(0) = mat.col(0);
      Dx.col(1) = mat.col(2);
      Dx.col(2) = mat.col(3);
      Dx.col(3) = mat.col(4);
      Dy.col(0) = mat.col(0);
      Dy.col(1) = mat.col(3);
      Dy.col(2) = mat.col(1);
      Dy.col(3) = mat.col(4);
      Dz.col(0) = mat.col(0);
      Dz.col(1) = mat.col(1);
      Dz.col(2) = mat.col(2);
      Dz.col(3) = mat.col(4);
      
      circumcenter = Vec3d(Dx.determinant(), Dy.determinant(), Dz.determinant()) / (2 * alpha);
      
      VertexHandle vh = tomesh->addVertex();
      vdpos[vh] = circumcenter;
      
      dt_tet_2_vd_vert[*tit] = vh;
      vd_vert_2_dt_tet[vh] = *tit;
      dt_vert_2_vd_verts[a].push_back(vh);
      dt_vert_2_vd_verts[b].push_back(vh);
      dt_vert_2_vd_verts[c].push_back(vh);
      dt_vert_2_vd_verts[d].push_back(vh);
    }
    
    // face in dt = edge in vd
    // create a vd edge for each dt face
    FaceProperty<EdgeHandle> dt_face_2_vd_edge(m_mesh);
    EdgeProperty<FaceHandle> vd_edge_2_dt_face(tomesh);
    for (FaceIterator fit = m_mesh->faces_begin(); fit != m_mesh->faces_end(); ++fit)
    {
      std::vector<TetHandle> tets;
      for (FaceTetIterator ftit = m_mesh->ft_iter(*fit); ftit; ++ftit)
        tets.push_back(*ftit);
      
      assert(tets.size() == 1 || tets.size() == 2);
      
      if (tets.size() == 2)
      {
        VertexHandle v0 = dt_tet_2_vd_vert[tets[0]];
        VertexHandle v1 = dt_tet_2_vd_vert[tets[1]];
        
        EdgeHandle eh = tomesh->addEdge(v0, v1);
        
        dt_face_2_vd_edge[*fit] = eh;
        vd_edge_2_dt_face[eh] = *fit;
      }
    }
    
    // edge in dt = face (polygon) in vd
    // triangulate the vd polygonal faces
//    EdgeProperty<std::vector<FaceHandle> > dt_edge_2_vd_faces(m_mesh);
    EdgeProperty<std::vector<VertexHandle> > dt_edge_2_vd_face_verts(m_mesh);
    EdgeProperty<std::vector<EdgeHandle> >   dt_edge_2_vd_face_edges(m_mesh);
//    FaceProperty<EdgeHandle> vd_face_2_dt_edge(tomesh); // this is not injective
    VertexProperty<char> cell_closed(m_mesh);
    cell_closed.assign(true);
    for (EdgeIterator eit = m_mesh->edges_begin(); eit != m_mesh->edges_end(); ++eit)
    {
      std::vector<FaceHandle> dt_faces;
      bool vd_face_incomplete = false;
      for (EdgeFaceIterator efit = m_mesh->ef_iter(*eit); efit; ++efit) 
      {
        dt_faces.push_back(*efit);
        if (!dt_face_2_vd_edge[*efit].isValid())
          vd_face_incomplete = true;
      }
      
      if (vd_face_incomplete)
      {
        cell_closed[m_mesh->fromVertex(*eit)] = false;
        cell_closed[m_mesh->toVertex(*eit)] = false;
        
        continue; // this vd polygon is semi-infinite
      }
      
      assert(dt_faces.size() >= 3);
      
      // reorder edges (vertices) to form a ring
      std::vector<VertexHandle> vd_verts;
      EdgeHandle vd_e0 = dt_face_2_vd_edge[dt_faces[0]];
      VertexHandle vd_v0 = tomesh->fromVertex(vd_e0);
      vd_verts.push_back(vd_v0);
      EdgeHandle vd_e = vd_e0;
      VertexHandle vd_v = tomesh->toVertex(vd_e0);
      while (vd_v != vd_v0)
      {
        vd_verts.push_back(vd_v);
        
        size_t i = 0;
        for (i = 0; i < dt_faces.size(); i++)
        {
          EdgeHandle eh = dt_face_2_vd_edge[dt_faces[i]];
          if (eh == vd_e)
            continue;
          
          if (tomesh->fromVertex(eh) == vd_v)
          {
            vd_v = tomesh->toVertex(eh);
            vd_e = eh;
            break;
          }
          if (tomesh->toVertex(eh) == vd_v)
          {
            vd_v = tomesh->fromVertex(eh);
            vd_e = eh;
            break;
          }
        }
        assert(i < dt_faces.size());
      }
      assert(vd_verts.size() == dt_faces.size());
      
      // find edges from vertices
      std::vector<EdgeHandle> vd_edges;
      for (size_t i = 0; i < vd_verts.size(); i++)
        vd_edges.push_back(findEdge(*tomesh, vd_verts[i], vd_verts[(i + 1) % vd_verts.size()]));
      
      // create faces
      std::vector<FaceHandle> vd_faces;
      VertexHandle v0 = vd_verts[0];
      for (size_t i = 1; i < dt_faces.size() - 1; i++)
      {
        VertexHandle v1 = vd_verts[i];
        VertexHandle v2 = vd_verts[i + 1];
        
//        FaceHandle fh = tomesh->addFace(v0, v1, v2);
//        
//        vd_faces.push_back(fh);
//        vd_face_2_dt_edge[fh] = *eit;
      }
      
//      dt_edge_2_vd_faces[*eit] = vd_faces;
      dt_edge_2_vd_face_verts[*eit] = vd_verts;
      dt_edge_2_vd_face_edges[*eit] = vd_edges;
    }

    /////////////////////////////////////////////////////////////////////////////////////////
    // Clip by bounding box
    
    assert(bbmin.x() < bbmax.x());
    assert(bbmin.y() < bbmax.y());
    assert(bbmin.z() < bbmax.z());
    
    // reshape every voronoi cell
    VertexProperty<std::vector<std::pair<EdgeHandle, EdgeHandle> > > vd_cell_clipping_new_polygon(m_mesh);
    
    // intersection of the polyhedron's interior with the halfspace of one of the walls of the bounding box
    // making use of the convexity of both voronoi cells and the bounding box
    for (int wall = 0; wall < 6; wall++)
    {
      // create a common vertes to shoot edges to from any voronoi cell that gets clipped, and thus needs a new polygon face generated to seal the hole
      VertexHandle common_vh = m_mesh->addVertex();
      
      Vec3d hsnormal = Vec3d((wall % 3 == 0 ? 1.0 : 0.0), (wall % 3 == 1 ? 1.0 : 0.0), (wall % 3 == 2 ? 1.0 : 0.0)) * (wall < 3 ? 1 : -1);
      Scalar hsposition = hsnormal.dot(wall < 3 ? bbmin : bbmax);
      
      // find the intersection points of all vd edges with this bounding box wall
      EdgeProperty<VertexHandle> intersections_with_walls(tomesh);      
      for (EdgeIterator eit = tomesh->edges_begin(); eit != tomesh->edges_end(); ++eit)
      {
        Vec3d x0 = vdpos[tomesh->fromVertex(*eit)];
        Vec3d x1 = vdpos[tomesh->toVertex(*eit)];
        
        if ((x0.dot(hsnormal) - hsposition) * (x1.dot(hsnormal) - hsposition) < 0)
        {
          VertexHandle intersectionvertex = tomesh->addVertex();
          intersections_with_walls[*eit] = intersectionvertex;
          vdpos[intersectionvertex] = x0 + (hsposition - x0.dot(hsnormal)) / (x1 - x0).dot(hsnormal) * (x1 - x0);
        }
      }

      // intersect with one halfspace      
      // iterate through all polygons, clip the ones intersecting the wall, and delete the ones outside the wall
      for (EdgeIterator eit = m_mesh->edges_begin(); eit != m_mesh->edges_end(); ++eit)
      {
        EdgeHandle polygon = *eit;  // edge in dt
        
        if (dt_edge_2_vd_face_edges[polygon].size() == 0)
          continue; // semi-infinite polygon
        
        EdgeHandle edge0;           // edge in vd
        EdgeHandle edge1;           // edge in vd
        bool intersection_found = false;
        for (size_t i = 0; i < dt_edge_2_vd_face_edges[polygon].size(); i++)
        {
          EdgeHandle vd_edge = dt_edge_2_vd_face_edges[polygon][i];
          
          Vec3d x0 = vdpos[tomesh->fromVertex(vd_edge)];
          Vec3d x1 = vdpos[tomesh->toVertex(vd_edge)];
          
          if ((x0.dot(hsnormal) - hsposition) * (x1.dot(hsnormal) - hsposition) < 0)
          {
            // edge vd_edge intersects this wall of the bounding box
            intersection_found = true;
            if (edge0.isValid() && edge1.isValid())
              assert(!"More than two edges of a Voronoi polygon face intersecting this bounding box wall.");
            else if (edge0.isValid())
              edge1 = vd_edge;
            else
              edge0 = vd_edge;
          }
        }
        
        if (intersection_found)
        {
          // clip the current polygon with the wall
          std::vector<EdgeHandle> new_polygon_edges;
          std::vector<EdgeHandle> edges_to_examine = dt_edge_2_vd_face_edges[polygon];
          for (size_t i = 0; i < edges_to_examine.size(); i++)
          {
            EdgeHandle eh = edges_to_examine[i];
            
            VertexHandle v0 = tomesh->fromVertex(eh);
            VertexHandle v1 = tomesh->toVertex(eh);
            
            Vec3d x0 = vdpos[v0];
            Vec3d x1 = vdpos[v1];
            
            if (eh == edge0 || eh == edge1)
            {
              VertexHandle vintersection = intersections_with_walls[eh];
              assert(vintersection.isValid());
              
              assert((x0.dot(hsnormal) - hsposition) * (x1.dot(hsnormal) - hsposition) < 0);

              if (x0.dot(hsnormal) < hsposition)
              {
                // make sure v0 is inside, and v1 is outside
                std::swap(v0, v1);
                std::swap(x0, x1);
              }
              
              // create a new edge to the intersection point
              EdgeHandle newedge = findEdge(*tomesh, v0, vintersection);
              if (!newedge.isValid())
                newedge = tomesh->addEdge(v0, vintersection);
              new_polygon_edges.push_back(newedge);
              vd_edge_2_dt_face[newedge] = FaceHandle();
              
              // remove this edge from the current polygon
              std::vector<EdgeHandle> & edges = dt_edge_2_vd_face_edges[polygon];
              edges.erase(std::remove(edges.begin(), edges.end(), eh), edges.end());
              
              // check to see if the edge is present in any other polygon. if not, delete it, and update primal graph accordingly for the new edge.
              bool edge_in_use = false;
              for (FaceEdgeIterator feit = m_mesh->fe_iter(vd_edge_2_dt_face[eh]); feit; ++feit)
              {
                std::vector<EdgeHandle> & edges = dt_edge_2_vd_face_edges[*feit];
                if (std::find(edges.begin(), edges.end(), eh) != edges.end())
                {
                  edge_in_use = true;
                  break;
                }
              }
              
              if (!edge_in_use)
              {
                // record the old edge's environment - it will become the environment of the new edge
                FaceHandle old_face = vd_edge_2_dt_face[eh];
                std::vector<std::pair<std::vector<FaceHandle>, VertexHandle> > tets_to_create;
                for (FaceTetIterator ftit = m_mesh->ft_iter(old_face); ftit; ++ftit)
                {
                  std::vector<FaceHandle> old_tet_faces;
                  for (TetFaceIterator tfit = m_mesh->tf_iter(*ftit); tfit; ++tfit)
                    old_tet_faces.push_back(*tfit);
                  old_tet_faces.erase(std::remove(old_tet_faces.begin(), old_tet_faces.end(), old_face), old_tet_faces.end());
                  assert(old_tet_faces.size() == 3);
                  
                  tets_to_create.push_back(std::pair<std::vector<FaceHandle>, VertexHandle>(old_tet_faces, dt_tet_2_vd_vert[*ftit]));
                }
                
                std::vector<EdgeHandle> face_to_create;
                for (FaceEdgeIterator feit = m_mesh->fe_iter(old_face); feit; ++feit)
                  face_to_create.push_back(*feit);
                
                // delete old edge from primal (dt)
                std::vector<TetHandle> tets_to_delete;
                for (FaceTetIterator ftit = m_mesh->ft_iter(old_face); ftit; ++ftit)
                  tets_to_delete.push_back(*ftit);
                for (size_t i = 0; i < tets_to_delete.size(); i++)
                  m_mesh->deleteTet(tets_to_delete[i], false);
                
                m_mesh->deleteFace(old_face, false);

                // delete old edge from dual (vd)
                tomesh->deleteEdge(eh, true);
                
                // update primal (dt) for the new edge
                assert(face_to_create.size() == 3);
                FaceHandle new_dt_face = m_mesh->addFace(face_to_create[0], face_to_create[1], face_to_create[2]);
                
                dt_face_2_vd_edge[new_dt_face] = newedge;
                vd_edge_2_dt_face[newedge] = new_dt_face;
                
                for (size_t i = 0; i < tets_to_create.size(); i++)
                {
                  assert(tets_to_create[i].first.size() == 3);
                  TetHandle new_dt_tet = m_mesh->addTet(tets_to_create[i].first[0], tets_to_create[i].first[1], tets_to_create[i].first[2], new_dt_face);
                  
                  dt_tet_2_vd_vert[new_dt_tet] = tets_to_create[i].second;
                  vd_vert_2_dt_tet[tets_to_create[i].second] = new_dt_tet;
                }
                
              }
              
            } else
            {
              assert((x0.dot(hsnormal) - hsposition) * (x1.dot(hsnormal) - hsposition) > 0);
              
              if (x0.dot(hsnormal) > hsposition)
              {
                assert(x1.dot(hsnormal) > hsposition);
                
                // entire edge is inside the halfspace - it's good
                new_polygon_edges.push_back(eh);
              } else
              {
                assert(x1.dot(hsnormal) < hsposition);
                
                // entire edge is outside the halfspace - need to delete it (removing it from its incident polygons too)
                for (FaceEdgeIterator feit = m_mesh->fe_iter(vd_edge_2_dt_face[eh]); feit; ++feit)
                {
                  std::vector<EdgeHandle> & edges = dt_edge_2_vd_face_edges[*feit];
                  edges.erase(std::remove(edges.begin(), edges.end(), eh), edges.end());
                }
                tomesh->deleteEdge(eh, true);                
              }
            }
          }
          
          EdgeHandle newboundaryedge = tomesh->addEdge(intersections_with_walls[edge0], intersections_with_walls[edge1]); // the corresponding primal face will be created later
          EdgeHandle polygon_of_newboundaryedge = polygon;
          
          new_polygon_edges.push_back(newboundaryedge);
          dt_edge_2_vd_face_edges[polygon] = new_polygon_edges;
          vd_edge_2_dt_face[newboundaryedge] = FaceHandle();
          
          vd_cell_clipping_new_polygon[m_mesh->fromVertex(polygon)].push_back(std::pair<EdgeHandle, EdgeHandle>(newboundaryedge, polygon_of_newboundaryedge));
          vd_cell_clipping_new_polygon[m_mesh->toVertex(polygon)  ].push_back(std::pair<EdgeHandle, EdgeHandle>(newboundaryedge, polygon_of_newboundaryedge));
          
        } else
        {
          assert(dt_edge_2_vd_face_edges[polygon].size() > 0);
          VertexHandle v0 = tomesh->fromVertex(dt_edge_2_vd_face_edges[polygon][0]);
          Vec3d x0 = vdpos[v0];
          
          if (x0.dot(hsnormal) < hsposition)
          {
            // this vertex is out of the halfspace. this means the entire polygon must be too, since it doesn't intersect the plane.
            // need to delete this polygon
            
            // delete it from dual graph (vd)
            std::vector<EdgeHandle> edges_to_delete = dt_edge_2_vd_face_edges[polygon]; // make a copy because dt_edge_2_vd_face_edges may be changed during deletion
            for (size_t i = 0; i < edges_to_delete.size(); i++)
            {
              EdgeHandle eh = edges_to_delete[i];  
              for (FaceEdgeIterator feit = m_mesh->fe_iter(vd_edge_2_dt_face[eh]); feit; ++feit)
              {
                std::vector<EdgeHandle> & edges = dt_edge_2_vd_face_edges[*feit];
                edges.erase(std::remove(edges.begin(), edges.end(), eh), edges.end());
              }
              tomesh->deleteEdge(eh, true);
            }
            
            // delete it from primal graph (dt)
            std::vector<TetHandle> tets_to_delete;
            for (EdgeFaceIterator efit = m_mesh->ef_iter(polygon); efit; ++efit)
              for (FaceTetIterator ftit = m_mesh->ft_iter(*efit); ftit; ++ftit)
                tets_to_delete.push_back(*ftit);
            for (size_t i = 0; i < tets_to_delete.size(); i++)
              if (m_mesh->tetExists(tets_to_delete[i]))
                m_mesh->deleteTet(tets_to_delete[i], false);
            
            std::vector<FaceHandle> faces_to_delete;
            for (EdgeFaceIterator efit = m_mesh->ef_iter(polygon); efit; ++efit)
              faces_to_delete.push_back(*efit);
            for (size_t i = 0; i < faces_to_delete.size(); i++)
              m_mesh->deleteFace(faces_to_delete[i], true);
            
          }
          
        }
        
      }
      
      // modify both primal (dt) and dual (vd) to reflect the clipping
      for (VertexIterator vit = m_mesh->vertices_begin(); vit != m_mesh->vertices_end(); ++vit)
      {
        VertexHandle cell = *vit;
        
        if (!cell_closed[cell])
          continue;
        
        if (vd_cell_clipping_new_polygon[cell].size() > 0)
        {
          // this cell has been subject to clipping, i.e. this dt vertex needs a new edge to seal the hole
          EdgeHandle new_polygon = findEdge(*m_mesh, cell, common_vh);
          if (!new_polygon.isValid())
            new_polygon = m_mesh->addEdge(cell, common_vh);  // edge in dt
          dt_edge_2_vd_face_edges[new_polygon].clear();

          std::set<VertexHandle> vd_verts;
          for (size_t i = 0; i < vd_cell_clipping_new_polygon[cell].size(); i++)
          {
            EdgeHandle new_vd_edge =      vd_cell_clipping_new_polygon[cell][i].first;   // edge in vd
            EdgeHandle adjacent_polygon = vd_cell_clipping_new_polygon[cell][i].second;  // edge in dt

            assert(m_mesh->fromVertex(adjacent_polygon) == cell || m_mesh->toVertex(adjacent_polygon) == cell);
            
            VertexHandle adjacent_cell = getEdgesOtherVertex(*m_mesh, adjacent_polygon, cell);
                        
            FaceHandle new_dt_face = findFace(*m_mesh, cell, adjacent_cell, common_vh);
            if (!new_dt_face.isValid())
              new_dt_face = m_mesh->addFace(cell, adjacent_cell, common_vh);
            
            dt_edge_2_vd_face_edges[new_polygon].push_back(new_vd_edge);
            dt_face_2_vd_edge[new_dt_face] = new_vd_edge;
            vd_edge_2_dt_face[new_vd_edge] = new_dt_face;
            
            vd_verts.insert(tomesh->fromVertex(new_vd_edge));
            vd_verts.insert(tomesh->toVertex(new_vd_edge));
          }
          
          // create the tets in dt
          for (std::set<VertexHandle>::iterator i = vd_verts.begin(); i != vd_verts.end(); i++)
          {
            VertexHandle vh = *i;
            
            EdgeHandle existing_edge;
            EdgeHandle new_vd_edge_0;
            EdgeHandle new_vd_edge_1;
            for (VertexEdgeIterator veit = tomesh->ve_iter(vh); veit; ++veit)
            {
              if (!vd_edge_2_dt_face[*veit].isValid())
                continue;
              
              if (faceContainsVertex(*m_mesh, vd_edge_2_dt_face[*veit], common_vh) && faceContainsVertex(*m_mesh, vd_edge_2_dt_face[*veit], cell))
              {
                if (new_vd_edge_0.isValid() && new_vd_edge_1.isValid())
                  assert(!"Error");
                else if (new_vd_edge_0.isValid())
                  new_vd_edge_1 = *veit;
                else
                  new_vd_edge_0 = *veit;
              } else if (faceContainsVertex(*m_mesh, vd_edge_2_dt_face[*veit], cell))
              {
                existing_edge = *veit;
              }
            }
            
            assert(existing_edge.isValid());
            assert(new_vd_edge_0.isValid());
            assert(new_vd_edge_1.isValid());
            
            EdgeHandle edge_between_adjacent_cells = getVertexOppositeEdgeInFace(*m_mesh, vd_edge_2_dt_face[existing_edge], cell);
            FaceHandle new_dt_face_to_common_vh = m_mesh->addFace(common_vh, m_mesh->fromVertex(edge_between_adjacent_cells), m_mesh->toVertex(edge_between_adjacent_cells));
            
            TetHandle new_dt_tet = findTet(*m_mesh, vd_edge_2_dt_face[existing_edge], vd_edge_2_dt_face[new_vd_edge_0], vd_edge_2_dt_face[new_vd_edge_1], new_dt_face_to_common_vh); 
            if (!new_dt_tet.isValid())
              new_dt_tet = m_mesh->addTet(vd_edge_2_dt_face[existing_edge], vd_edge_2_dt_face[new_vd_edge_0], vd_edge_2_dt_face[new_vd_edge_1], new_dt_face_to_common_vh);
            
            dt_face_2_vd_edge[new_dt_face_to_common_vh] = EdgeHandle(); // this edge is semi infinite
            dt_tet_2_vd_vert[new_dt_tet] = vh;
            vd_vert_2_dt_tet[vh] = new_dt_tet;
            
          }
          
          vd_cell_clipping_new_polygon[cell].clear();
          
        }
        
      }      
      
    }

    /////////////////////////////////////////////////////////////////////////////////////////
    // Triangulation of voronoi polygon faces
    
    // triangulate the vd polygonal faces
    EdgeProperty<std::vector<FaceHandle> > dt_edge_2_vd_faces(m_mesh);
    FaceProperty<EdgeHandle> vd_face_2_dt_edge(tomesh); // this is not injective
    for (EdgeIterator eit = m_mesh->edges_begin(); eit != m_mesh->edges_end(); ++eit)
    {
      EdgeHandle polygon = *eit;
      
      if (dt_edge_2_vd_face_edges[polygon].size() == 0)
        continue;
      
      // reorder edges (vertices) to form a ring
      std::vector<VertexHandle> vd_verts;
      EdgeHandle vd_e0 = dt_edge_2_vd_face_edges[polygon][0];
      VertexHandle vd_v0 = tomesh->fromVertex(vd_e0);
      vd_verts.push_back(vd_v0);
      EdgeHandle vd_e = vd_e0;
      VertexHandle vd_v = tomesh->toVertex(vd_e0);
      while (vd_v != vd_v0)
      {
        vd_verts.push_back(vd_v);
        
        size_t i = 0;
        for (i = 0; i < dt_edge_2_vd_face_edges[polygon].size(); i++)
        {
          // search in the polygon for the next edge to grow the ring
          EdgeHandle eh = dt_edge_2_vd_face_edges[polygon][i];
          if (eh == vd_e)
            continue;
          
          if (tomesh->fromVertex(eh) == vd_v)
          {
            vd_v = tomesh->toVertex(eh);
            vd_e = eh;
            break;
          }
          if (tomesh->toVertex(eh) == vd_v)
          {
            vd_v = tomesh->fromVertex(eh);
            vd_e = eh;
            break;
          }
        }
        assert(i < dt_edge_2_vd_face_edges[polygon].size());
      }
      assert(vd_verts.size() == dt_edge_2_vd_face_edges[polygon].size());
            
      // create faces
      std::vector<FaceHandle> vd_faces;
      VertexHandle v0 = vd_verts[0];
      for (size_t i = 1; i < vd_verts.size() - 1; i++)
      {
        VertexHandle v1 = vd_verts[i];
        VertexHandle v2 = vd_verts[i + 1];
        
        FaceHandle fh = tomesh->addFace(v0, v1, v2);
        
        vd_faces.push_back(fh);
        vd_face_2_dt_edge[fh] = *eit;
      }
      
      dt_edge_2_vd_faces[*eit] = vd_faces;
    }

    /////////////////////////////////////////////////////////////////////////////////////////
    // Prune orphan primitives
    std::vector<EdgeHandle> edges_to_delete;
    for (EdgeIterator eit = tomesh->edges_begin(); eit != tomesh->edges_end(); ++eit)
      if (tomesh->edgeIncidentFaces(*eit) == 0)
        edges_to_delete.push_back(*eit);
    
    for (size_t i = 0; i < edges_to_delete.size(); i++)
      tomesh->deleteEdge(edges_to_delete[i], false);
    
    std::vector<VertexHandle> vertices_to_delete;
    for (VertexIterator vit = tomesh->vertices_begin(); vit != tomesh->vertices_end(); ++vit)
      if (tomesh->vertexIncidentEdges(*vit) == 0)
        vertices_to_delete.push_back(*vit);
    
    for (size_t i = 0; i < vertices_to_delete.size(); i++)
      tomesh->deleteVertex(vertices_to_delete[i]);

//    outputDT();
//    outputVD(tomesh, vdpos, dt_edge_2_vd_face_edges);
    

    return true;
  }
  
  Scalar DelaunayTriangulator::predicateOriented(const Vec3d & a, const Vec3d & b, const Vec3d & c, const Vec3d & p)
  {
    typedef Eigen::Matrix<Scalar, 4, 4> Mat4d;
    
    Mat4d mat = Mat4d::Ones();
    mat.block<1, 3>(0, 0) = a;
    mat.block<1, 3>(1, 0) = b;
    mat.block<1, 3>(2, 0) = c;
    mat.block<1, 3>(3, 0) = p;
    
    return -mat.determinant();  // The Ledoux paper uses left-hand rule. Here we use right hand rule.
  }
  
  Scalar DelaunayTriangulator::predicateInSphere(const Vec3d & a, const Vec3d & b, const Vec3d & c, const Vec3d & d, const Vec3d & p)
  {
    Scalar orientation = predicateOriented(a, b, c, d);
    int sign = (orientation > 0 ? 1 : (orientation < 0 ? -1 : 0));
    
    typedef Eigen::Matrix<Scalar, 5, 5> Mat5d;
    
    Mat5d mat = Mat5d::Ones();
    mat.block<1, 3>(0, 0) = a;  mat(0, 3) = a.squaredNorm();
    mat.block<1, 3>(1, 0) = b;  mat(1, 3) = b.squaredNorm();
    mat.block<1, 3>(2, 0) = c;  mat(2, 3) = c.squaredNorm();
    mat.block<1, 3>(3, 0) = d;  mat(3, 3) = d.squaredNorm();
    mat.block<1, 3>(4, 0) = p;  mat(4, 3) = p.squaredNorm();
    
    return sign * mat.determinant();   // return positive if p is outside the out-sphere of tet abcd; negative if inside; zero if on.s. regardless of the orientation of abcd
  }

  bool DelaunayTriangulator::predicateInTetrahedron(const Vec3d & a, const Vec3d & b, const Vec3d & c, const Vec3d & d, const Vec3d & p)
  {
    Scalar abcd = predicateOriented(a, b, c, d);
    Scalar pbcd = predicateOriented(p, b, c, d);
    Scalar apcd = predicateOriented(a, p, c, d);
    Scalar abpd = predicateOriented(a, b, p, d);
    Scalar abcp = predicateOriented(a, b, c, p);
    
    if (abcd >= 0 && pbcd >= 0 && apcd >= 0 && abpd >= 0 && abcp >= 0)
      return true;
    
    if (abcd <= 0 && pbcd <= 0 && apcd <= 0 && abpd <= 0 && abcp <= 0)
      return true;
    
    return false;
  }
  
  Scalar DelaunayTriangulator::predicateOriented(VertexHandle a, VertexHandle b, VertexHandle c, VertexHandle p) const
  {
    return predicateOriented(pos(a), pos(b), pos(c), pos(p));
  }
  
  Scalar DelaunayTriangulator::predicateInSphere(VertexHandle a, VertexHandle b, VertexHandle c, VertexHandle d, VertexHandle p) const
  {
    return predicateInSphere(pos(a), pos(b), pos(c), pos(d), pos(p));
  }
  
  bool DelaunayTriangulator::predicateInTetrahedron(VertexHandle a, VertexHandle b, VertexHandle c, VertexHandle d, VertexHandle p) const
  {
    return predicateInTetrahedron(pos(a), pos(b), pos(c), pos(d), pos(p));
  }
 
  void DelaunayTriangulator::outputVD(Mesh * tomesh, VertexProperty<Vec3d> & vdpos, EdgeProperty<std::vector<EdgeHandle> > & dt_edge_2_vd_face_edges)
  {
    std::cout << "VD mesh" << std::endl;
    for (VertexIterator it = tomesh->vertices_begin(); it != tomesh->vertices_end(); ++it)
    {
      if (tomesh->vertexIncidentEdges(*it) > 0)
        std::cout << "vertex " << (*it).idx() << ": " << vdpos[*it] << std::endl;
    }
    for (EdgeIterator it = tomesh->edges_begin(); it != tomesh->edges_end(); ++it)
    {
      std::cout << "edge " << (*it).idx() << ": vertices: ";
      for (EdgeVertexIterator evit = tomesh->ev_iter(*it); evit; ++evit)
        std::cout << (*evit).idx() << " ";
      std::cout << std::endl;
    }
    for (EdgeIterator it = m_mesh->edges_begin(); it != m_mesh->edges_end(); ++it)
    {
      std::cout << "polygon " << (*it).idx() << ": edges: ";
      for (size_t i = 0; i < dt_edge_2_vd_face_edges[*it].size(); i++)
        std::cout << dt_edge_2_vd_face_edges[*it][i].idx() << " ";
      std::cout << std::endl;
    }
  }
  
  void DelaunayTriangulator::outputDT()
  {
    std::cout << "DT mesh" << std::endl;
    for (EdgeIterator it = m_mesh->edges_begin(); it != m_mesh->edges_end(); ++it)
    {
      std::cout << "edge " << (*it).idx() << ": " << m_mesh->fromVertex(*it).idx() << " " << m_mesh->toVertex(*it).idx() << " ";
      int count = 0;
      for (EdgeFaceIterator efit = m_mesh->ef_iter(*it); efit; ++efit)
        count++;
      std::cout << "face count = " << count << std::endl;
    }
    for (FaceIterator it = m_mesh->faces_begin(); it != m_mesh->faces_end(); ++it)
    {
      std::cout << "face " << (*it).idx() << ": vertices: ";
      for (FaceVertexIterator fvit = m_mesh->fv_iter(*it); fvit; ++fvit)
        std::cout << (*fvit).idx() << " ";
      std::cout << " edges: ";
      for (FaceEdgeIterator feit = m_mesh->fe_iter(*it); feit; ++feit)
        std::cout << (*feit).idx() << " ";
      std::cout << std::endl;
    }
    for (TetIterator it = m_mesh->tets_begin(); it != m_mesh->tets_end(); ++it)
    {
      std::cout << "tet " << (*it).idx() << ": vertices: ";
      for (TetVertexIterator tvit = m_mesh->tv_iter(*it); tvit; ++tvit)
        std::cout << (*tvit).idx() << " ";
      std::cout << " faces: ";
      for (TetFaceIterator tfit = m_mesh->tf_iter(*it); tfit; ++tfit)
        std::cout << (*tfit).idx() << " ";
      std::cout << std::endl;
    }
  }

}