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
      
      if (predicateInSphere(a, b, c, p, d))
      {
        // local non-delaunay, need flpping
        // algorithm in figure 12 of Ledoux paper        
        
        // determine the intersection of line pd with facd abc
        Vec3d nabc = (pos(b) - pos(a)).cross(pos(c) - pos(a));
        Vec3d cabc = (pos(a) + pos(b) + pos(c)) / 3;
        Scalar t = (cabc - pos(p)).dot(nabc) / (pos(d) - pos(p)).dot(nabc);
        assert(t == t);
        
        bool inside = predicateInTetrahedron(Vec3d(0, 0, 0), pos(a), pos(b), pos(c), pos(p) + (pos(d) - pos(p)) * t);
        
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
          TetHandle pabd = findTet(*m_mesh, p, a, b, d);
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
    
    m_mesh->deleteTet(pabc, true);
    m_mesh->deleteTet(abcd, true);
    
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
    
    m_mesh->deleteTet(pabc, true);
    m_mesh->deleteTet(abcd, true);
    m_mesh->deleteTet(pabd, true);
    
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

  bool DelaunayTriangulator::extractVoronoiDiagram(Mesh * tomesh)
  {
    
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
    
    Mat4d mat = Mat4d::Ones();
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