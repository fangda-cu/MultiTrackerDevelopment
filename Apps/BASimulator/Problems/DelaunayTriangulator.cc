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
    // vertex in dt = region in vd
    if (false)
    {
      for (VertexIterator vit = m_mesh->vertices_begin(); vit != m_mesh->vertices_end(); ++vit)
      {
        VertexHandle vh = tomesh->addVertex();
        vdpos[vh] = pos(*vit);
      }
      
      for (EdgeIterator eit = m_mesh->edges_begin(); eit != m_mesh->edges_end(); ++eit)
      {
        EdgeHandle eh = tomesh->addEdge(m_mesh->fromVertex(*eit), m_mesh->toVertex(*eit));
      }
    }
    
    // debugging output
    if (true)
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
        continue;
      
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
    }

    // debugging output
    if (true)
    {
      std::cout << "VD mesh" << std::endl;
      for (VertexIterator it = tomesh->vertices_begin(); it != tomesh->vertices_end(); ++it)
      {
        std::cout << "vertex " << (*it).idx() << ": " << vdpos[*it] << std::endl;
      }
      for (FaceIterator it = tomesh->faces_begin(); it != tomesh->faces_end(); ++it)
      {
        std::cout << "face " << (*it).idx() << ": vertices: ";
        for (FaceVertexIterator fvit = tomesh->fv_iter(*it); fvit; ++fvit)
          std::cout << (*fvit).idx() << " ";
        std::cout << std::endl;
      }
//      for (EdgeIterator it = m_mesh->edges_begin(); it != m_mesh->edges_end(); ++it)
//      {
//        std::cout << "polygon " << (*it).idx() << ": faces: ";
//        for (size_t i = 0; i < dt_edge_2_vd_faces[*it].size(); i++)
//          std::cout << dt_edge_2_vd_faces[*it][i].idx() << " ";
//        std::cout << std::endl;
//      }
    }
    
    /////////////////////////////////////////////////////////////////////////////////////////
    // Clip by bounding box
    
    assert(bbmin.x() < bbmax.x());
    assert(bbmin.y() < bbmax.y());
    assert(bbmin.z() < bbmax.z());
    
    // find the intersection points of all vd edges with bounding box walls
    EdgeProperty<VertexHandle> * intersections_with_walls[6];
    for (int wall = 0; wall < 6; wall++)
    {
      intersections_with_walls[wall] = new EdgeProperty<VertexHandle>(tomesh);
      Vec3d hsnormal = Vec3d((wall % 3 == 0 ? 1.0 : 0.0), (wall % 3 == 1 ? 1.0 : 0.0), (wall % 3 == 2 ? 1.0 : 0.0)) * (wall < 3 ? 1 : -1);
      Scalar hsposition = hsnormal.dot(wall < 3 ? bbmin : bbmax);
      
      for (EdgeIterator eit = tomesh->edges_begin(); eit != tomesh->edges_end(); ++eit)
      {
        Vec3d x0 = vdpos[tomesh->fromVertex(*eit)];
        Vec3d x1 = vdpos[tomesh->toVertex(*eit)];
        
        if ((x0.dot(hsnormal) - hsposition) * (x1.dot(hsnormal) - hsposition) < 0)
        {
          VertexHandle intersectionvertex = tomesh->addVertex();
          (*intersections_with_walls[wall])[*eit] = intersectionvertex;
          vdpos[intersectionvertex] = x0 + (hsposition - x0.dot(hsnormal)) / (x1 - x0).dot(hsnormal) * (x1 - x0);
        }
      }
    }
    
//    EdgeProperty<int> first_intersection(tomesh);
//    for (EdgeIterator eit = tomesh->edges_begin(); eit != tomesh->edges_end(); ++eit)
//    {
//      VertexHandle v0 = tomesh->fromVertex(*eit);
//      VertexHandle v1 = tomesh->toVertex(*eit);
//      
//      Vec3d x0 = vdpos[v0];
//      Vec3d x1 = vdpos[v1];
//      
//      bool inside0 = (x0.x() > bbmin.x() && x0.x() < bbmax.x() && x0.y() > bbmin.y() && x0.y() < bbmax.y() && x0.z() > bbmin.z() && x0.z() < bbmax.z());
//      bool inside1 = (x1.x() > bbmin.x() && x1.x() < bbmax.x() && x1.y() > bbmin.y() && x1.y() < bbmax.y() && x1.z() > bbmin.z() && x1.z() < bbmax.z());
//      if (inside0 != inside1)
//      {
//        // this edge intersects the bounding box walls - find the intersection
//        if (inside1)
//        {
//          std::swap(v0, v1);
//          std::swap(x0, x1);
//          std::swap(inside0, inside1);
//        }
//        assert(inside0);
//        assert(!inside1);
//        
//        Scalar d[6];
//        d[0] = (bbmin.x() - x0.x()) / (x1.x() - x0.x());
//        d[1] = (bbmin.y() - x0.y()) / (x1.y() - x0.y());
//        d[2] = (bbmin.z() - x0.z()) / (x1.z() - x0.z());
//        d[3] = (bbmax.x() - x0.x()) / (x1.x() - x0.x());
//        d[4] = (bbmax.y() - x0.y()) / (x1.y() - x0.y());
//        d[5] = (bbmax.z() - x0.z()) / (x1.z() - x0.z());
//        
//        // add an intersection vertex for the intersection between this edge and each bounding box wall (if they intersect).
//        // also find the first intersection going from the interior endpoint x0
//        Scalar firstd = std::numeric_limits<Scalar>::infinity();
//        int firstdi = -1;
//        for (int i = 0; i < 6; i++)
//        {
//          if (d[i] >= 0 && d[i] <= 1 && d[i] < firstd)
//          {
//            firstd = d[i];
//            firstdi = i;
//          }
//        }
//        assert(firstd >= 0 && firstd <= 1);
//        assert(firstdi >= 0);
//        
//        first_intersection[*eit] = firstdi;
//      }
//    }
    
    // reshape every voronoi cell
    EdgeProperty<std::vector<VertexHandle> > dt_edge_2_vd_face_verts_new(m_mesh);
    dt_edge_2_vd_face_verts_new = dt_edge_2_vd_face_verts;
    EdgeProperty<std::vector<EdgeHandle> >   dt_edge_2_vd_face_edges_new(m_mesh);
    dt_edge_2_vd_face_edges_new = dt_edge_2_vd_face_edges;
    
    // intersection of the polyhedron's interior with the halfspace of one of the walls of the bounding box
    // making use of the convexity of both voronoi cells and the bounding box
    for (int wall = 0; wall < 6; wall++)
    {
      Vec3d hsnormal = Vec3d((wall % 3 == 0 ? 1.0 : 0.0), (wall % 3 == 1 ? 1.0 : 0.0), (wall % 3 == 2 ? 1.0 : 0.0)) * (wall < 3 ? 1 : -1);
      Scalar hsposition = hsnormal.dot(wall < 3 ? bbmin : bbmax);
      
      // intersect with one halfspace
      
      // this type represents a vd edge that is cut by a bounding box wall.
      // the first edge is an edge in dt representing the vd polygon
      // the other two edges are the two vd edges in the vd polygon that are cut by the wall
      typedef Eigen::Matrix<EdgeHandle, 3, 1> EdgeHandle3;
      std::vector<EdgeHandle3> open_intersected_edges;

      // find all the polygons that intersect the wall
      EdgeHandle intersected_polygon;         // edge in dt
      EdgeHandle intersected_polygon_edge_0;  // edges in vd
      EdgeHandle intersected_polygon_edge_1;
      for (EdgeIterator eit = m_mesh->edges_begin(); eit != m_mesh->edges_end(); ++eit)
      {
        EdgeHandle polygon = *eit;
        
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
            intersected_polygon = polygon;
            
            if (intersected_polygon_edge_0.isValid() && intersected_polygon_edge_1.isValid())
              assert(!"More than two edges of a Voronoi polygon face intersecting this bounding box wall.");
            else if (intersected_polygon_edge_0.isValid())
              intersected_polygon_edge_1 = vd_edge;
            else
              intersected_polygon_edge_0 = vd_edge;
          }
        }
        
        if (intersection_found)
          open_intersected_edges.push_back(EdgeHandle3(intersected_polygon, intersected_polygon_edge_0, intersected_polygon_edge_1));
      }
      
      // a search that traces the loop of intersection between voronoi cell surface and bounding box walls, processing the clipping of polygons along the way
      EdgeProperty<bool> visited(m_mesh);
      visited.assign(false);
      
      while (open_intersected_edges.size() > 0)
      {
        EdgeHandle3 current = open_intersected_edges.back();
        open_intersected_edges.pop_back();
        
        EdgeHandle polygon = current.x();
        EdgeHandle edge0 = current.y();
        EdgeHandle edge1 = current.y();
        
        if (visited[polygon])
          continue;
        
        // clip the current polygon with the wall
        std::vector<EdgeHandle> new_polygon_edges;
        for (size_t i = 0; i < dt_edge_2_vd_face_edges[polygon].size(); i++)
        {
          EdgeHandle eh = dt_edge_2_vd_face_edges[polygon][i];
          
          VertexHandle v0 = tomesh->fromVertex(eh);
          VertexHandle v1 = tomesh->toVertex(eh);
          
          Vec3d x0 = vdpos[v0];
          Vec3d x1 = vdpos[v1];
          
          if (eh == edge0 || eh == edge1)
          {
            VertexHandle vintersection = (*intersections_with_walls[wall])[eh];
            assert(vintersection.isValid());
            
            assert((x0.dot(hsnormal) - hsposition) * (x1.dot(hsnormal) - hsposition) < 0);

            if (x0.dot(hsnormal) < hsposition)
            {
              // make sure v0 is inside, and v1 is outside
              std::swap(v0, v1);
              std::swap(x0, x1);
            }
            
            EdgeHandle newedge = tomesh->addEdge(v0, vintersection);
            new_polygon_edges.push_back(newedge);
            
          } else
          {
            assert((x0.dot(hsnormal) - hsposition) * (x1.dot(hsnormal) - hsposition) > 0);
            
            if (x0.dot(hsnormal) > hsposition)
            {
              assert(x1.dot(hsnormal) > hsposition);
              
              // entire edge is inside the halfspace - it's good
              new_polygon_edges.push_back(eh);
            }
          }
        }
        
        EdgeHandle newboundaryedge = tomesh->addEdge((*intersections_with_walls[wall])[edge0], (*intersections_with_walls[wall])[edge1]);
        new_polygon_edges.push_back(newboundaryedge);
        
        dt_edge_2_vd_face_edges_new[polygon] = new_polygon_edges;
        visited[polygon] = true;
        
        // find the polygons incident to these cut edges, because they are also cut by this wall
        std::vector<EdgeHandle> next_polygons;
        for (FaceEdgeIterator feit = m_mesh->fe_iter(vd_edge_2_dt_face[edge0]); feit; ++feit)
          if (!visited[*feit])
            next_polygons.push_back(*feit);
        for (FaceEdgeIterator feit = m_mesh->fe_iter(vd_edge_2_dt_face[edge1]); feit; ++feit)
          if (!visited[*feit])
            next_polygons.push_back(*feit);
        
        for (size_t i = 0; i < next_polygons.size(); i++)
        {
          EdgeHandle next_polygon = next_polygons[i];
          EdgeHandle next_edge_0;
          EdgeHandle next_edge_1;
          
          for (size_t j = 0; j < dt_edge_2_vd_face_edges[next_polygon].size(); j++)
          {
            EdgeHandle vd_edge = dt_edge_2_vd_face_edges[next_polygon][j];
            
            Vec3d x0 = vdpos[tomesh->fromVertex(vd_edge)];
            Vec3d x1 = vdpos[tomesh->toVertex(vd_edge)];
            
            if ((x0.dot(hsnormal) - hsposition) * (x1.dot(hsnormal) - hsposition) < 0)
            {
              // edge vd_edge intersects this wall of the bounding box
              if (next_edge_0.isValid() && next_edge_1.isValid())
                assert(!"More than two edges of a Voronoi polygon face intersecting this bounding box wall.");
              else if (next_edge_0.isValid())
                next_edge_1 = vd_edge;
              else
                next_edge_0 = vd_edge;
            }
          }
          assert(next_edge_0.isValid());
          assert(next_edge_1.isValid());
          
          open_intersected_edges.push_back(EdgeHandle3(next_polygon, next_edge_0, next_edge_1));
        }
      }
      
    }
    
    

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
  
}