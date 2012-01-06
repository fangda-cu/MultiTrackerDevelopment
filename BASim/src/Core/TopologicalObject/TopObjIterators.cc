#include "BASim/src/Core/TopologicalObject/TopObjIterators.hh"
#include "BASim/src/Core/TopologicalObject/TopologicalObject.hh"

namespace BASim {


  //* Vertex Iterator routines */

  VertexIterator& VertexIterator::operator++ ()
  {
    do {
      m_hnd.__increment();
    } while(m_hnd.idx() < (int)m_obj->numVertexSlots() && !m_obj->vertexExists(m_hnd));
    return *this;
  }

  VertexIterator& VertexIterator::operator-- ()
  {
    do {
      m_hnd.__decrement();
    } while(m_hnd.idx() >= 0 && !m_obj->vertexExists(m_hnd));
    return *this;
  }

  //* Edge Iterator routines */

  EdgeIterator& EdgeIterator::operator++ ()
  {
    do {
      m_hnd.__increment();
    } while(m_hnd.idx() < (int)m_obj->numEdgeSlots() && !m_obj->edgeExists(m_hnd));
    return *this;
  }

  EdgeIterator& EdgeIterator::operator-- ()
  {
    do {
      m_hnd.__decrement();
    } while (m_hnd.idx() >= 0 && !m_obj->edgeExists(m_hnd));
    return *this;
  }

  //* Face Iterator routines */

  FaceIterator& FaceIterator::operator++ ()
  {
    do {
      m_hnd.__increment();
    } while(m_hnd.idx() < (int)m_obj->numFaceSlots() && !m_obj->faceExists(m_hnd));
    return *this;
  }

  FaceIterator& FaceIterator::operator-- ()
  {
    do {
      m_hnd.__decrement();
    } while (m_hnd.idx() >= 0 && m_obj->faceExists(m_hnd));
    return *this;
  }

  //* Tet Iterator routines */

  TetIterator& TetIterator::operator++ ()
  {
    do {
      m_hnd.__increment();
    } while(m_hnd.idx() < (int)m_obj->numTetSlots() && !m_obj->tetExists(m_hnd));
    return *this;
  }

  TetIterator& TetIterator::operator-- ()
  {
    do {
      m_hnd.__decrement();
    } while (m_hnd.idx() >= 0 && m_obj->tetExists(m_hnd));
    return *this;
  }

  //* VertexEdge iterator routines */

  VertexEdgeIterator::VertexEdgeIterator() : m_obj(NULL), m_hnd(-1), m_idx(0) {}

  VertexEdgeIterator::VertexEdgeIterator(const VertexEdgeIterator& veit)
    : m_obj(veit.m_obj)
    , m_hnd(veit.m_hnd)
    , m_idx(0)
  {}

  VertexEdgeIterator::VertexEdgeIterator(const TopologicalObject* t, const VertexHandle& vh)
    : m_obj(t)
    , m_hnd(vh)
    , m_idx(0)
  {}

  VertexEdgeIterator& VertexEdgeIterator::operator= (const VertexEdgeIterator& veit)
  {
    m_obj = veit.m_obj;
    m_hnd = veit.m_hnd;
    m_idx = veit.m_idx;
    return *this;
  }


  bool VertexEdgeIterator::operator== (const VertexEdgeIterator& veit) const
  {
    assert(m_obj);
    return m_obj == veit.m_obj
      && m_hnd == veit.m_hnd 
      && m_idx == veit.m_idx; 
  }

  bool VertexEdgeIterator::operator!= (const VertexEdgeIterator& veit) const
  {
    return !(operator==(veit));
  }

  VertexEdgeIterator& VertexEdgeIterator::operator++ ()
  {
    ++m_idx;
    return *this;
  }

  VertexEdgeIterator& VertexEdgeIterator::operator-- ()
  {
    --m_idx;
    return *this;
  }

  VertexEdgeIterator::value_type VertexEdgeIterator::operator* ()
  {
    assert(m_obj);
    assert(m_obj->vertexExists(m_hnd));

    if(m_idx >= m_obj->m_VE.getNumEntriesInRow(m_hnd.idx()))
      return EdgeHandle(-1);
    return EdgeHandle(m_obj->m_VE.getColByIndex(m_hnd.idx(), m_idx));
  }


  VertexEdgeIterator::operator bool() const
  {
    assert(m_obj);
    assert(m_obj->vertexExists(m_hnd));
    return m_idx < m_obj->m_VE.getNumEntriesInRow(m_hnd.idx());
  }

  //* EdgeVertex circulator routines */

  EdgeVertexIterator::EdgeVertexIterator() : m_obj(NULL), m_hnd(-1), m_idx(0) {}

  EdgeVertexIterator::EdgeVertexIterator(const EdgeVertexIterator& evit)
    : m_obj(evit.m_obj)
    , m_hnd(evit.m_hnd)
    , m_idx(0)
  {}

  EdgeVertexIterator::EdgeVertexIterator(const TopologicalObject* t, const EdgeHandle& eh)
    : m_obj(t)
    , m_hnd(eh)
    , m_idx(0)
  {
    assert(m_obj);    
    assert(m_obj->edgeExists(eh));
  }
  
  EdgeVertexIterator& EdgeVertexIterator::operator= (const EdgeVertexIterator& evit)
  {
    m_obj = evit.m_obj;
    m_hnd = evit.m_hnd;
    m_idx = evit.m_idx;
    return *this;
  }

  bool EdgeVertexIterator::operator== (const EdgeVertexIterator& evit) const
  {
    assert(m_obj);
    assert(evit.m_obj);

    return ((m_obj == evit.m_obj)
      && (m_hnd == evit.m_hnd) 
      && (m_idx == evit.m_idx));
  }

  bool  EdgeVertexIterator::operator!= (const EdgeVertexIterator& evit) const
  {
    return !(operator==(evit));
  }

  EdgeVertexIterator&  EdgeVertexIterator::operator++ ()
  {
    ++m_idx;
    return *this;
  }

  EdgeVertexIterator&  EdgeVertexIterator::operator-- ()
  {
    --m_idx;
    return *this;
  }

  EdgeVertexIterator::value_type EdgeVertexIterator::operator* ()
  {
    assert(m_obj);
    assert(m_obj->edgeExists(m_hnd));

    if(m_idx < 0 || m_idx >= m_obj->m_EV.getNumEntriesInRow(m_hnd.idx()))
      return VertexHandle(-1);

    //Return in directional order: 'from' vertex = 0, 'to' vertex = 1
    //There are only two vertices to check.
    int targetSign = m_idx == 0 ? -1 : 1; 
    int vert_index = (m_obj->m_EV.getValueByIndex(m_hnd.idx(), 0) == targetSign)?
      m_obj->m_EV.getColByIndex(m_hnd.idx(), 0) :
      m_obj->m_EV.getColByIndex(m_hnd.idx(), 1);

    return VertexHandle(vert_index);
          
    //Return in lexicographic order
    //return VertexHandle(m_obj->m_EV.getColByIndex(m_hnd.idx(), m_idx));
  }

 
  EdgeVertexIterator::operator bool() const
  {
    assert(m_obj);
    
    return m_idx < m_obj->m_EV.getNumEntriesInRow(m_hnd.idx());
  }

  //* FaceEdge circulator routines */

  FaceEdgeIterator::FaceEdgeIterator() : m_obj(NULL), m_hnd(-1), m_curEdge(-1), m_idx(0)  {}

  FaceEdgeIterator::FaceEdgeIterator(const FaceEdgeIterator& fvit)
    : m_obj(fvit.m_obj)
    , m_hnd(fvit.m_hnd)
    , m_curEdge(fvit.m_curEdge)
    , m_idx(fvit.m_idx)
  {}

  FaceEdgeIterator::FaceEdgeIterator(const TopologicalObject* t, const FaceHandle& fh)
    : m_obj(t)
    , m_hnd(fh)
    , m_curEdge(EdgeHandle(t->m_FE.getColByIndex(fh.idx(), 0)))
    , m_idx(0)
  {
    assert(m_obj);    
    assert(m_obj->faceExists(fh));
  }

  FaceEdgeIterator& FaceEdgeIterator::operator= (const FaceEdgeIterator& fvit)
  {
    m_obj = fvit.m_obj;
    m_hnd = fvit.m_hnd;
    m_curEdge = fvit.m_curEdge;
    m_idx = fvit.m_idx;
    return *this;
  }

  bool FaceEdgeIterator::operator== (const FaceEdgeIterator& fvit) const
  {
    assert(m_obj);
    assert(fvit.m_obj);

    return ((m_obj == fvit.m_obj)
      && (m_hnd == fvit.m_hnd) 
      && (m_curEdge == fvit.m_curEdge)
      && (m_idx == fvit.m_idx));
  }

  bool  FaceEdgeIterator::operator!= (const FaceEdgeIterator& fvit) const
  {
    return !(operator==(fvit));
  }

  FaceEdgeIterator&  FaceEdgeIterator::operator++ ()
  {
    ++m_idx;
    m_curEdge = m_obj->nextEdge(m_hnd, m_curEdge, 1);
    return *this;
  }

  FaceEdgeIterator&  FaceEdgeIterator::operator-- ()
  {
    --m_idx;
    m_curEdge = m_obj->nextEdge(m_hnd, m_curEdge, -1);
    return *this;
  }

  FaceEdgeIterator::value_type FaceEdgeIterator::operator* ()
  {
    assert(m_obj);
    assert(m_obj->faceExists(m_hnd));

    if(m_idx < 0 || m_idx >= (int)m_obj->m_FE.getNumEntriesInRow(m_hnd.idx()))
      return EdgeHandle(-1);
    
    return m_curEdge;
  }


  FaceEdgeIterator::operator bool() const
  {
    assert(m_obj);
    return m_idx < (int)m_obj->m_FE.getNumEntriesInRow(m_hnd.idx());
  }

  //* VertexVertex circulator routines */

  VertexVertexIterator::VertexVertexIterator() : m_obj(NULL), m_hnd(-1), m_idx(0) {}

  VertexVertexIterator::VertexVertexIterator(const VertexVertexIterator& vvit)
    : m_obj(vvit.m_obj)
    , m_hnd(vvit.m_hnd)
    , m_idx(0)
  {}

  VertexVertexIterator::VertexVertexIterator(const TopologicalObject* t, const VertexHandle& vh)
    : m_obj(t)
    , m_hnd(vh)
    , m_idx(0)
  {}
  
  VertexVertexIterator& VertexVertexIterator::operator= (const VertexVertexIterator& vvit)
  {
    m_obj = vvit.m_obj;
    m_hnd = vvit.m_hnd;
    m_idx = vvit.m_idx;
    return *this;
  }

  bool VertexVertexIterator::operator== (const VertexVertexIterator& vvit) const
  {
    assert(m_obj);
    assert(m_obj->vertexExists(m_hnd));

    return m_obj == vvit.m_obj
      && m_hnd == vvit.m_hnd 
      && m_idx == vvit.m_idx;
  }

  bool VertexVertexIterator::operator!= (const VertexVertexIterator& vvit) const
  {
    return !(operator==(vvit));
  }

  VertexVertexIterator& VertexVertexIterator::operator++ ()
  {
    ++m_idx;
    return *this;
  }

  VertexVertexIterator& VertexVertexIterator::operator-- ()
  {
    --m_idx;
    return *this;
  }

  VertexVertexIterator::value_type VertexVertexIterator::operator* ()
  {
    assert(m_obj);
    assert(m_obj->vertexExists(m_hnd));

    if(m_idx >=  m_obj->m_VE.getNumEntriesInRow(m_hnd.idx()))
      return VertexHandle(-1);

    //get nbr edge
    unsigned int e_idx = m_obj->m_VE.getColByIndex(m_hnd.idx(), m_idx);
    
    //get vertices that constitute it
    int v0 = m_obj->m_EV.getColByIndex(e_idx, 0);
    int v1 = m_obj->m_EV.getColByIndex(e_idx, 1);
    
    assert(v0 == m_hnd.idx() || v1 == m_hnd.idx());

    //pick the one that's not the central vertex
    return VertexHandle(m_hnd.idx() == v0? v1 : v0);
  }


  VertexVertexIterator::operator bool() const
  {
    assert(m_obj);
    assert(m_obj->vertexExists(m_hnd));
    return m_idx < m_obj->m_VE.getNumEntriesInRow(m_hnd.idx());
  }



  //* FaceVertex circulator routines */

  FaceVertexIterator::FaceVertexIterator() : m_feit() {}

  FaceVertexIterator::FaceVertexIterator(const FaceVertexIterator& fvit)
    : m_feit(fvit.m_feit)
    , m_obj(fvit.m_obj)
    , m_hnd(fvit.m_hnd)
  {}

  FaceVertexIterator::FaceVertexIterator(const TopologicalObject* t, const FaceHandle& fh)
    : m_feit(t, fh)
    , m_obj(t)
    , m_hnd(fh)
  {
 
  }

  FaceVertexIterator& FaceVertexIterator::operator= (const FaceVertexIterator& fvit)
  {
    m_feit = fvit.m_feit;
    m_obj = fvit.m_obj;
    m_hnd = fvit.m_hnd;
    return *this;
  }

  bool FaceVertexIterator::operator== (const FaceVertexIterator& fvit) const
  {
    return m_feit == fvit.m_feit && m_obj == fvit.m_obj && m_hnd == fvit.m_hnd;
  }

  bool FaceVertexIterator::operator!= (const FaceVertexIterator& fvit) const
  {
    return !(operator==(fvit));
  }

  FaceVertexIterator& FaceVertexIterator::operator++ ()
  {
    ++m_feit;
    return *this;
  }

  FaceVertexIterator& FaceVertexIterator::operator-- ()
  {
    --m_feit;
    return *this;
  }

  FaceVertexIterator::value_type FaceVertexIterator::operator* ()
  {
    EdgeHandle edge = *m_feit;
    assert(m_obj->edgeExists(edge));

    int sign = m_obj->m_FE.get(m_hnd.idx(), edge.idx());

    return VertexHandle(sign == 1? m_obj->fromVertex(edge) : m_obj->toVertex(edge));
  }

 
  FaceVertexIterator::operator bool() const
  {
    assert(m_obj);
    return m_feit;
  }


  //* VertexFace circulator routines */

  VertexFaceIterator::VertexFaceIterator() : m_obj(NULL), m_hnd(-1), m_idx(0) {}

  VertexFaceIterator::VertexFaceIterator(const VertexFaceIterator& vfit)
    : m_obj(vfit.m_obj)
    , m_hnd(vfit.m_hnd)
    , m_idx(vfit.m_idx)
  {}

  VertexFaceIterator::VertexFaceIterator(const TopologicalObject* t, const VertexHandle& vh)
    : m_obj(t)
    , m_hnd(vh)
    , m_idx(0)
  {
    if(!m_obj->m_validVF)
      m_obj->compute_nbrsVF();

  }

  VertexFaceIterator& VertexFaceIterator::operator= (const VertexFaceIterator& vfit)
  {
    m_obj = vfit.m_obj;
    m_hnd = vfit.m_hnd;
    m_idx = vfit.m_idx;
    return *this;
  }

  bool VertexFaceIterator::operator== (const VertexFaceIterator& vfit) const
  {
    assert(m_obj);
    assert(vfit.m_obj);
    return m_obj == vfit.m_obj
      && m_hnd == vfit.m_hnd
      && m_idx == vfit.m_idx;
  }

  bool VertexFaceIterator::operator!= (const VertexFaceIterator& vfit) const
  {
    return !(operator==(vfit));
  }

  VertexFaceIterator& VertexFaceIterator::operator++ ()
  {
    ++m_idx;
    return *this;
  }

  VertexFaceIterator& VertexFaceIterator::operator-- ()
  {
    --m_idx;
    return *this;
  }

  VertexFaceIterator::value_type VertexFaceIterator::operator* ()
  {
    assert(m_obj);
    if(m_idx >= m_obj->m_nbrsVF.getNumEntriesInRow(m_hnd.idx()))
      return FaceHandle(-1);
    return FaceHandle(m_obj->m_nbrsVF.getColByIndex(m_hnd.idx(), m_idx));
  }


  VertexFaceIterator::operator bool() const
  {
    assert(m_obj);
    return m_idx < m_obj->m_nbrsVF.getNumEntriesInRow(m_hnd.idx());
  }

  //* EdgeFace iterator routines */

  EdgeFaceIterator::EdgeFaceIterator() : m_obj(NULL), m_hnd(-1), m_idx(0) {}

  EdgeFaceIterator::EdgeFaceIterator(const EdgeFaceIterator& efit)
    : m_obj(efit.m_obj)
    , m_hnd(efit.m_hnd)
    , m_idx(efit.m_idx)
  {}

  EdgeFaceIterator::EdgeFaceIterator(const TopologicalObject* t, const EdgeHandle& eh)
    : m_obj(t)
    , m_hnd(eh)
    , m_idx(0)
  {}

  EdgeFaceIterator& EdgeFaceIterator::operator= (const EdgeFaceIterator& efit)
  {
    m_obj = efit.m_obj;
    m_hnd = efit.m_hnd;
    m_idx = efit.m_idx;
    return *this;
  }


  bool EdgeFaceIterator::operator== (const EdgeFaceIterator& efit) const
  {
    assert(m_obj);
    return m_obj == efit.m_obj
      && m_hnd == efit.m_hnd 
      && m_idx == efit.m_idx; 
  }

  bool EdgeFaceIterator::operator!= (const EdgeFaceIterator& efit) const
  {
    return !(operator==(efit));
  }

  EdgeFaceIterator& EdgeFaceIterator::operator++ ()
  {
    ++m_idx;
    return *this;
  }

  EdgeFaceIterator& EdgeFaceIterator::operator-- ()
  {
    --m_idx;
    return *this;
  }

  EdgeFaceIterator::value_type EdgeFaceIterator::operator* ()
  {
    assert(m_obj);
    assert(m_obj->edgeExists(m_hnd));

    if(m_idx >= m_obj->m_EF.getNumEntriesInRow(m_hnd.idx()) || m_idx < 0)
      return FaceHandle(-1);
    return FaceHandle(m_obj->m_EF.getColByIndex(m_hnd.idx(), m_idx));
  }


  EdgeFaceIterator::operator bool() const
  {
    assert(m_obj);
    assert(m_obj->edgeExists(m_hnd));
    return m_idx < m_obj->m_EF.getNumEntriesInRow(m_hnd.idx());
  }

}