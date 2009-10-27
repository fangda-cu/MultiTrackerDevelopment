/**
 * \file Property.hh
 *
 * \author miklos@cs.columbia.edu
 * \date 08/29/2009
 */

#ifndef PROPERTY_HH
#define PROPERTY_HH

namespace BASim {

/** Base class for storing properties. A property essentially just
    represents any data associated to a certain object. \sa
    Property */
class PropertyBase
{
public:

  PropertyBase(const std::string& name) : m_name(name) {}

  virtual ~PropertyBase() {}

  virtual void reserve(size_t n) = 0;

  virtual void resize(size_t n) = 0;

  virtual void push_back() = 0;

  const std::string& name() const { return m_name; }

protected:

  std::string m_name;
};

/** Handle for referring to properties. */
template <typename T>
class PropertyHandleBase : public HandleBase
{
public:

  explicit PropertyHandleBase(int idx = -1) : HandleBase(idx) {}
};

/** A generic property class using an array-like interface for storing
    multiple properties of the same type. A property can have any type
    and a name associated to it. */
template <typename T>
class Property : public PropertyBase
{
public:

  typedef T                                       Value;
  typedef std::vector<T>                          vector_type;
  typedef T                                       value_type;
  typedef typename vector_type::reference         reference;
  typedef typename vector_type::const_reference   const_reference;
  typedef typename vector_type::pointer           pointer;
  typedef typename vector_type::const_pointer     const_pointer;
  typedef typename vector_type::iterator          iterator;

  Property(const std::string& name)
    : PropertyBase(name)
  {}

  virtual void reserve(size_t n)
  {
    m_data.reserve(n);
  }

  virtual void resize(size_t n)
  {
    m_data.resize(n);
  }

  virtual void push_back()
  {
    m_data.push_back(T());
  }

  const_pointer data() const
  {
    if (m_data.empty()) return NULL;
    return &m_data[0];
  }

  reference operator[] (size_t idx)
  {
    assert(idx < m_data.size());
    return m_data[idx];
  }

  const_reference operator[] (size_t idx) const
  {
    assert(idx < m_data.size());
    return m_data[idx];
  }

  reference operator[] (const HandleBase& h)
  {
    assert(h.isValid());
    assert(h.idx() >= 0);
    assert((size_t) h.idx() < m_data.size());
    return m_data[h.idx()];
  }

  const_reference operator[] (const HandleBase& h) const
  {
    assert(h.isValid());
    assert(h.idx() >= 0);
    assert((size_t)  h.idx() < m_data.size());
    return m_data[h.idx()];
  }

  void set_default(const T& def)
  {
    for (iterator it = m_data.begin(); it != m_data.end(); ++it) *it = def;
  }

protected:

  vector_type m_data;
};

/** Container for properties. */
class PropertyContainer
{
public:

  typedef std::vector<PropertyBase*>   Properties;
  typedef Properties::iterator         iterator;

  PropertyContainer() {}

  ~PropertyContainer()
  {
    iterator pIt;
    for (pIt = m_properties.begin(); pIt != m_properties.end(); ++pIt) {
      delete (*pIt);
    }
  }

  template <class T> PropertyHandleBase<T>
  add(const T&, const std::string& name = "<unknown>")
  {
    int idx = m_properties.size();
    m_properties.push_back(new Property<T>(name));

    return PropertyHandleBase<T>(idx);
  }

  template <class T> Property<T>&
  property(const PropertyHandleBase<T>& h)
  {
    assert(h.idx() >= 0);
    assert((size_t) h.idx() < m_properties.size());
    assert(m_properties[h.idx()] != NULL);

    return *smart_cast<Property<T>*>(m_properties[h.idx()]);
  }

  template <class T> const Property<T>&
  property(const PropertyHandleBase<T>& h) const
  {
    assert(h.idx() >= 0);
    assert((size_t) h.idx() < m_properties.size());
    assert(m_properties[h.idx()] != NULL);

    return *smart_cast<const Property<T>*>(m_properties[h.idx()]);
  }

  template <class T>
  bool exists(const T&, const std::string& name) const
  {
    for (size_t i = 0; i < m_properties.size(); ++i) {
      if (m_properties[i]->name() == name) {
        Property<T>* test = dynamic_cast<Property<T>*>(m_properties[i]);
        if (test != NULL) return true;
      }
    }
    return false;
  }

  template <class T> PropertyHandleBase<T>
  handle(const T&, const std::string& name) const
  {
    assert(exists(T(), name));

    for (size_t i = 0; i < m_properties.size(); ++i) {
      if (m_properties[i]->name() == name) {
        Property<T>* test = dynamic_cast<Property<T>*>(m_properties[i]);
        if (test != NULL) return PropertyHandleBase<T>((int) i);
      }
    }

    std::cerr << "Property " << name << " not found in container"
              << std::endl;
    return PropertyHandleBase<T>(m_properties.size());
  }

  void resize(size_t size)
  {
    iterator it;
    for (it = m_properties.begin(); it != m_properties.end(); ++it) {
      (*it)->resize(size);
    }
  }

protected:

  Properties m_properties;
};

} // namespace BASim

#endif // PROPERTY_HH
