/**
 * \file DegreeOfFreedom.hh
 *
 * \author miklos@cs.columbia.edu
 * \date 09/17/2009
 */

#ifndef DEGREEOFFREEDOM_HH
#define DEGREEOFFREEDOM_HH

namespace BASim {
// TODO: Why does DofHandle have a HandleBase and extend form HandleBase?
//       It is not like its decorating methods or anything... look into this further

/** Handle for referring to a degree of freedom.  */
class DofHandle : public HandleBase
{
public:

  /// Degrees of freedom can be associated to vertices, edges
  enum Type { VERTEX_DOF, EDGE_DOF };

  explicit DofHandle(int idx = -1)
    : HandleBase(idx)
    , m_dofNum(-1)
  {}

  DofHandle(const DofHandle& handle)
    : HandleBase(handle.m_idx)
    , m_type(handle.m_type)
    , m_handle(handle.m_handle)
    , m_dofNum(handle.m_dofNum)
  {}

  DofHandle& operator= (const DofHandle& handle)
  {
    m_idx = handle.m_idx;
    m_type = handle.m_type;
    m_handle = handle.m_handle;
    m_dofNum = handle.m_dofNum;
    return *this;
  }

  bool operator<( const DofHandle& rhs ) const
  {
    if( this->getType() != rhs.getType() ) return this->getType() < rhs.getType();
    if( this->getNum() != rhs.getNum() ) return this->getNum() < rhs.getNum();
    if( this->getHandle() != rhs.getHandle() ) return this->getHandle() < rhs.getHandle();
    // Should this be here???
    if( this->m_idx != rhs.m_idx ) return this->m_idx < rhs.m_idx;
    return false;
  }

  bool operator==( const DofHandle& rhs ) const
  {
    if( this->getType() != rhs.getType() ) return false;
    if( this->getNum() != rhs.getNum() ) return false;
    if( this->getHandle() != rhs.getHandle() ) return false;
    // Should this be here???
    if( this->m_idx != rhs.m_idx ) return false;
    return true;
  }
  
  bool operator!=( const DofHandle& rhs ) const
  {
    return !((*this)==rhs);
  }

  Type getType() const { return m_type; }
  void setType(const Type& type) { m_type = type; }

  int getNum() const { return m_dofNum; }
  void setNum(int num) { m_dofNum = num; }

  const HandleBase& getHandle() const { return m_handle; }
  void setHandle(const HandleBase& handle) { m_handle = handle; }

protected:

  Type m_type;
  HandleBase m_handle;
  int m_dofNum; ///< which degree of freedom
};

/** Mapping from degrees of freedom to indices. */
class DOFMap
{
public:

  DOFMap() { assert( m_dofToIndex.size() == m_indexToDof.size() ); }
  ~DOFMap() {}

  void addMapping(const DofHandle& dof, int index)
  {
    assert( m_dofToIndex.size() == m_indexToDof.size() );
    m_dofToIndex.insert(std::make_pair(dof, index));
    m_indexToDof.insert(std::make_pair(index, dof));
    assert( m_dofToIndex.size() == m_indexToDof.size() );
  }

  const DofHandle& getDof(int index) const
  {
    std::map<int, DofHandle>::const_iterator it = m_indexToDof.find(index);
    assert(it != m_indexToDof.end());
    return it->second;
  }

  int getIndex(const DofHandle& handle) const
  {
    std::map<DofHandle, int>::const_iterator it = m_dofToIndex.find(handle);
    assert(it != m_dofToIndex.end());
    return it->second;
  }

  const std::map<DofHandle,int>& getDofToIndexMap() const
  {
    return m_dofToIndex;
  }

  const std::map<int,DofHandle>& getIndexToDofMap() const
  {
    return m_indexToDof;
  }

  void clearMappings()
  {
    m_dofToIndex.clear();
    m_indexToDof.clear();
    assert( m_dofToIndex.size() == 0 );
    assert( m_indexToDof.size() == 0 );
  }
  
  size_t size() const
  {
    assert( m_dofToIndex.size() == m_indexToDof.size() );
    return m_dofToIndex.size();
  }

protected:

  std::map<DofHandle, int> m_dofToIndex;
  std::map<int, DofHandle> m_indexToDof;
};

} // namespace BASim

#endif // DEGREEOFFREEDOM_HH
