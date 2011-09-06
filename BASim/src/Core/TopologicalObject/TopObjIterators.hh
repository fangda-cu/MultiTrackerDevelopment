/**
 * \file TopObjIterators.hh
 *
 * \author miklos@cs.columbia.edu
 * \date 09/13/2009
 * \author smith@cs.columbia.edu
 * \date 03/14/2010
 * \author batty@cs.columbia.edu
 * \date 04/10/2011
 */

#ifndef TOPOBJITERATORS_HH
#define TOPOBJITERATORS_HH

#include "BASim/src/Core/Definitions.hh"
#include "BASim/src/Core/TopologicalObject/TopObjHandles.hh"
#include "BASim/src/Core/TopologicalObject/Topology.hh"

namespace BASim {


class TopologicalObject;

/** Iterator over vertices */
class VertexIterator
{
public:

  typedef VertexHandle                    value_type;
  typedef value_type&                     reference;
  typedef const value_type&               const_reference;
  typedef value_type*                     pointer;
  typedef const value_type*               const_pointer;

  VertexIterator() : m_obj(NULL), m_hnd(-1) {}

  VertexIterator(const VertexIterator& vit)
    : m_obj(vit.m_obj)
    , m_hnd(vit.m_hnd)
  {}

  VertexIterator(const TopologicalObject* t, value_type val)
    : m_obj(t)
    , m_hnd(val)
  {}

  VertexIterator& operator= (const VertexIterator& vit)
  {
    m_obj = vit.m_obj;
    m_hnd = vit.m_hnd;
    return *this;
  }

  bool operator== (const VertexIterator& vit) const
  {
    return ((m_obj == vit.m_obj) && (m_hnd == vit.m_hnd));
  }

  bool operator!= (const VertexIterator& vit) const
  {
    return !(operator==(vit));
  }

  VertexIterator& operator++ ();

  VertexIterator& operator-- ();


  reference operator* ()
  {
    return m_hnd;
  }

  const_reference operator* () const
  {
    return m_hnd;
  }

  pointer operator-> ()
  {
    return &m_hnd;
  }

  const_pointer operator-> () const
  {
    return &m_hnd;
  }

protected:

  const TopologicalObject* m_obj;
  value_type m_hnd;
};

/** Iterator over edges */

class EdgeIterator
{
public:

  typedef EdgeHandle                value_type;
  typedef value_type&               reference;
  typedef const value_type&         const_reference;
  typedef value_type*               pointer;
  typedef const value_type*         const_pointer;

  EdgeIterator() : m_obj(NULL), m_hnd(-1) {}

  EdgeIterator(const EdgeIterator& eit)
    : m_obj(eit.m_obj)
    , m_hnd(eit.m_hnd)
  {}

  EdgeIterator(const TopologicalObject* t, value_type val)
    : m_obj(t)
    , m_hnd(val)
  {}

  EdgeIterator& operator= (const EdgeIterator& eit)
  {
    m_obj = eit.m_obj;
    m_hnd = eit.m_hnd;
    return *this;
  }

  bool operator== (const EdgeIterator& eit) const
  {
    return ((m_obj == eit.m_obj) && (m_hnd == eit.m_hnd));
  }

  bool operator!= (const EdgeIterator& eit) const
  {
    return !(operator==(eit));
  }

  EdgeIterator& operator++ ();

  EdgeIterator& operator-- ();

  reference operator* ()
  {
    return m_hnd;
  }

  const_reference operator* () const
  {
    return m_hnd;
  }

  pointer operator-> ()
  {
    return &m_hnd;
  }

  const_pointer operator-> () const
  {
    return &m_hnd;
  }

protected:

  const TopologicalObject* m_obj;
  value_type m_hnd;
};

/** Iterator over faces */
class FaceIterator
{
public:

  typedef FaceHandle                value_type;
  typedef value_type&               reference;
  typedef const value_type&         const_reference;
  typedef value_type*               pointer;
  typedef const value_type*         const_pointer;

  FaceIterator() : m_obj(NULL), m_hnd(-1) {}

  FaceIterator(const FaceIterator& fit)
    : m_obj(fit.m_obj)
    , m_hnd(fit.m_hnd)
  {}

  FaceIterator(const TopologicalObject* t, value_type val)
    : m_obj(t)
    , m_hnd(val)
  {}

  FaceIterator& operator= (const FaceIterator& fit)
  {
    m_obj = fit.m_obj;
    m_hnd = fit.m_hnd;
    return *this;
  }

  bool operator== (const FaceIterator& fit) const
  {
    return ((m_obj == fit.m_obj) && (m_hnd == fit.m_hnd));
  }

  bool operator!= (const FaceIterator& fit) const
  {
    return !(operator==(fit));
  }

  FaceIterator& operator++ ();
 
  FaceIterator& operator-- ();
 
  reference operator* ()
  {
    return m_hnd;
  }

  const_reference operator* () const
  {
    return m_hnd;
  }

  pointer operator-> ()
  {
    return &m_hnd;
  }

  const_pointer operator-> () const
  {
    return &m_hnd;
  }

protected:

  const TopologicalObject* m_obj;
  value_type m_hnd;
};

/** Iterator over tetrahedra */
class TetIterator
{
public:

   typedef TetHandle                 value_type;
   typedef value_type&               reference;
   typedef const value_type&         const_reference;
   typedef value_type*               pointer;
   typedef const value_type*         const_pointer;

   TetIterator() : m_obj(NULL), m_hnd(-1) {}

   TetIterator(const TetIterator& tit)
      : m_obj(tit.m_obj)
      , m_hnd(tit.m_hnd)
   {}

   TetIterator(const TopologicalObject* t, value_type val)
      : m_obj(t)
      , m_hnd(val)
   {}

   TetIterator& operator= (const TetIterator& tit)
   {
      m_obj = tit.m_obj;
      m_hnd = tit.m_hnd;
      return *this;
   }

   bool operator== (const TetIterator& tit) const
   {
      return ((m_obj == tit.m_obj) && (m_hnd == tit.m_hnd));
   }

   bool operator!= (const TetIterator& tit) const
   {
      return !(operator==(tit));
   }

   TetIterator& operator++ ();
   

   TetIterator& operator-- ();
  

   reference operator* ()
   {
      return m_hnd;
   }

   const_reference operator* () const
   {
      return m_hnd;
   }

   pointer operator-> ()
   {
      return &m_hnd;
   }

   const_pointer operator-> () const
   {
      return &m_hnd;
   }

protected:

   const TopologicalObject* m_obj;
   value_type m_hnd;
};


/** Iterator over edges adjacent to a vertex */
class VertexEdgeIterator
{
public:

  typedef EdgeHandle                  value_type;
  typedef value_type&                 reference;
  typedef const value_type&           const_reference;
  typedef value_type*                 pointer;
  typedef const value_type*           const_pointer;

  VertexEdgeIterator();
  VertexEdgeIterator(const VertexEdgeIterator& veit);
  VertexEdgeIterator(const TopologicalObject* t, const VertexHandle& vh);

  VertexEdgeIterator& operator= (const VertexEdgeIterator& veit);

  bool operator== (const VertexEdgeIterator& veit) const;
  bool operator!= (const VertexEdgeIterator& veit) const;

  VertexEdgeIterator& operator++ ();
  VertexEdgeIterator& operator-- ();

  value_type operator* ();

  operator bool() const;

protected:

  const TopologicalObject* m_obj; ///< Pointer to the topological object
  VertexHandle m_hnd;            ///< Handle to the vertex whose topology is being considered
  size_t m_idx;                   ///< Current index into topological information
};

/** Iterator over vertices adjacent to an edge */
class EdgeVertexIterator
{
public:

  typedef VertexHandle               value_type;
  typedef value_type&                 reference;
  typedef const value_type&           const_reference;
  typedef value_type*                 pointer;
  typedef const value_type*           const_pointer;

  EdgeVertexIterator();
  EdgeVertexIterator(const EdgeVertexIterator& evit);
  EdgeVertexIterator(const TopologicalObject* t, const EdgeHandle& eh);
 
  EdgeVertexIterator& operator= (const EdgeVertexIterator& evit);
 
  bool operator== (const EdgeVertexIterator& evit) const;
  bool operator!= (const EdgeVertexIterator& evit) const;
 
  EdgeVertexIterator& operator++ ();
  EdgeVertexIterator& operator-- ();
 
  value_type operator* ();
 
  operator bool() const;
 
protected:

  const TopologicalObject* m_obj; ///< Pointer to topological object
  EdgeHandle m_hnd;               ///< Handle to edge whose connectivity is being considered
  size_t m_idx;                   ///< Current index into topology information
  };

//** Iterator over edges adjacent to a face */
class FaceEdgeIterator
{
public:

  typedef EdgeHandle                  value_type;
  typedef value_type&                 reference;
  typedef const value_type&           const_reference;
  typedef value_type*                 pointer;
  typedef const value_type*           const_pointer;

  FaceEdgeIterator();
  FaceEdgeIterator(const FaceEdgeIterator& evit);
  FaceEdgeIterator(const TopologicalObject* t, const FaceHandle& eh);

  FaceEdgeIterator& operator= (const FaceEdgeIterator& evit);

  bool operator== (const FaceEdgeIterator& evit) const;
  bool operator!= (const FaceEdgeIterator& evit) const;

  FaceEdgeIterator& operator++ ();
  FaceEdgeIterator& operator-- ();

  value_type operator* ();

  operator bool() const;

protected:

  const TopologicalObject* m_obj; ///< Pointer to topological object
  FaceHandle m_hnd;               ///< Handle to face whose connectivity is being considered
  EdgeHandle m_curEdge;         ///< Handle to current edge
  int m_idx;                    ///< Index in the iterator, to prevent infinite cycling
};

/** Iterator over vertices in the one ring of a vertex */
class VertexVertexIterator
{
public:

  typedef VertexHandle               value_type;
  typedef value_type&                 reference;
  typedef const value_type&           const_reference;
  typedef value_type*                 pointer;
  typedef const value_type*           const_pointer;

  VertexVertexIterator();
  VertexVertexIterator(const VertexVertexIterator& vvit);
  VertexVertexIterator(const TopologicalObject* t, const VertexHandle& vh);

  VertexVertexIterator& operator= (const VertexVertexIterator& vvit);

  bool operator== (const VertexVertexIterator& vvit) const;
  bool operator!= (const VertexVertexIterator& vvit) const;

  VertexVertexIterator& operator++ ();
  VertexVertexIterator& operator-- ();
 
  value_type operator* ();

  operator bool() const;

protected:

  const TopologicalObject* m_obj;   ///< Pointer to topological object
  VertexHandle m_hnd;               ///< Handle to vertex whose one-ring is being considered
  size_t m_idx;                     ///< Current index into topology information
};

/** Iterator over vertices adjacent to a face */
class FaceVertexIterator
{
public:

  typedef VertexHandle               value_type;
  typedef value_type&                 reference;
  typedef const value_type&           const_reference;
  typedef value_type*                 pointer;
  typedef const value_type*           const_pointer;

  FaceVertexIterator();
  FaceVertexIterator(const FaceVertexIterator& vfit);
  FaceVertexIterator(const TopologicalObject* t, const FaceHandle& fh);

  FaceVertexIterator& operator= (const FaceVertexIterator& vfit);

  bool operator== (const FaceVertexIterator& vfit) const;
  bool operator!= (const FaceVertexIterator& vfit) const;
 
  FaceVertexIterator& operator++ ();
  FaceVertexIterator& operator-- ();

  value_type operator* ();

  operator bool() const;
  
protected:
  FaceEdgeIterator m_feit;
  const TopologicalObject* m_obj;
  FaceHandle m_hnd;
  /*
  const TopologicalObject* m_obj;        ///< Pointer to topological object
  FaceHandle m_hnd;     ///< Handle to face whose connectivity is being considered
  size_t m_idx;          ///< Current index into topology information
  */
};

/** Iterator over faces adjacent to a vertex */
class VertexFaceIterator
{
public:

  typedef FaceHandle                  value_type;
  typedef value_type&                 reference;
  typedef const value_type&           const_reference;
  typedef value_type*                 pointer;
  typedef const value_type*           const_pointer;

  VertexFaceIterator();
  VertexFaceIterator(const VertexFaceIterator& vfit);
  VertexFaceIterator(const TopologicalObject* t, const VertexHandle& fh);


  VertexFaceIterator& operator= (const VertexFaceIterator& vfit);

  bool operator== (const VertexFaceIterator& vfit) const;
  bool operator!= (const VertexFaceIterator& vfit) const;

  VertexFaceIterator& operator++ ();
  VertexFaceIterator& operator-- ();

  value_type operator* ();

  operator bool() const;

protected:

  const TopologicalObject* m_obj;        ///< Pointer to topological object
  VertexHandle m_hnd;     ///< Handle to vertex whose connectivity is being considered
  size_t m_idx;          ///< Current index into topology information
};

/** Iterator over edges adjacent to a vertex */
class EdgeFaceIterator
{
public:

  typedef FaceHandle                  value_type;
  typedef value_type&                 reference;
  typedef const value_type&           const_reference;
  typedef value_type*                 pointer;
  typedef const value_type*           const_pointer;

  EdgeFaceIterator();
  EdgeFaceIterator(const EdgeFaceIterator& efit);
  EdgeFaceIterator(const TopologicalObject* t, const EdgeHandle& eh);

  EdgeFaceIterator& operator= (const EdgeFaceIterator& veit);

  bool operator== (const EdgeFaceIterator& veit) const;
  bool operator!= (const EdgeFaceIterator& veit) const;

  EdgeFaceIterator& operator++ ();
  EdgeFaceIterator& operator-- ();

  value_type operator* ();

  operator bool() const;

protected:

  const TopologicalObject* m_obj; ///< Pointer to the topological object
  EdgeHandle m_hnd;               ///< Handle to the edge whose topology is being considered
  size_t m_idx;                   ///< Current index into topological information
};


} // namespace BASim

#endif // TOPOBJITERATORS_HH
