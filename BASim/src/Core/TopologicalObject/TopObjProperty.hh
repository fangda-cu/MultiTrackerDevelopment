#ifndef TOPOBJECTPROPERTY_H
#define TOPOBJECTPROPERTY_H

#include "BASim/src/Core/Definitions.hh"
#include "BASim/src/Core/Handle.hh"
#include "BASim/src/Core/TopologicalObject/TopologicalObject.hh"

namespace BASim {

//Base class, so TopologicalObject can store pointers to TopObjProperties of different types in a single list.
class TopObjPropertyBase {

public:

  TopObjPropertyBase(TopologicalObject* obj) : m_obj(obj) {}
  virtual ~TopObjPropertyBase() {}

protected:

  //Only the TopologicalObject that owns this object is allowed to manipulate its size, 
  //in order to match the number of the given simplex type. (eg. #of edge properties should equal #of edges)
  virtual size_t size() const = 0;
  virtual void resize(size_t n) = 0;

  //The simplex mesh this property is associated with.
  TopologicalObject* m_obj;

  friend class TopologicalObject;
};


//The templated object property that stores the data. This is a base class that should not be used, since it doesn't get registered with
//a particular simplex type (edge, face, etc.), and cannot be resized.
template <class T>
class TopObjProperty : public TopObjPropertyBase {

public:
  TopObjProperty(TopologicalObject* obj, size_t n) : TopObjPropertyBase(obj), m_data(n) {}

  //TopObjProperty(const TopObjProperty& ) {}

  virtual ~TopObjProperty() {}

  void assign(const T& data_value) { for(unsigned int i = 0; i < m_data.size(); ++i) m_data[i] = data_value; }
  
  T& operator[] (const HandleBase& h) { 
    assert(h.idx() >= 0 && h.idx() < (int)m_data.size());
    return m_data[h.idx()]; 
  }
  
  T const& operator[] (const HandleBase& h) const { 
    assert(h.idx() >= 0 && h.idx() < (int)m_data.size()); 
    return m_data[h.idx()]; 
  }

protected: 
  
  size_t size() const { return m_data.size(); }
  void resize(size_t n) { m_data.resize(n); }

  std::vector<T> m_data;
  
};

//The properties associated to simplices of different dimensions. They are instantiated with a pointer to the the TopologicalObject that
//"owns" them, and registered with that object, so they can be automatically resized to match the mesh when necessary.
template <class T>
class VertexProperty : public TopObjProperty<T> {

public:
  VertexProperty(TopologicalObject* obj) : TopObjProperty<T>(obj,obj->numVertexSlots()) {
    this->m_obj->registerVertexProperty(this);
  }

  explicit VertexProperty(const VertexProperty& prop) : TopObjProperty<T>(prop.m_obj,prop.m_obj->numVertexSlots()) {
    this->m_obj->registerVertexProperty(this);
    this->m_data = prop.m_data;
  }

  ~VertexProperty() {
    this->m_obj->removeVertexProperty(this);
  }

  VertexProperty& operator=(const VertexProperty& other) {

    if(this != &other) {
        this->m_obj->removeVertexProperty(this);

        this->m_obj = other.m_obj;
        this->m_data = other.m_data;

        this->m_obj->registerVertexProperty(this);
    }

    return *this;
  }

};


template <class T>
class EdgeProperty : public TopObjProperty<T> {
public:
  EdgeProperty(TopologicalObject* obj) : TopObjProperty<T>(obj,obj->numEdgeSlots()) {
    this->m_obj->registerEdgeProperty(this);
  }

  explicit EdgeProperty(const EdgeProperty& prop) : TopObjProperty<T>(prop.m_obj,prop.m_obj->numEdgeSlots()) {
    this->m_obj->registerEdgeProperty(this);
    this->m_data = prop.m_data;
  }

  ~EdgeProperty() {
    this->m_obj->removeEdgeProperty(this);
  }

  EdgeProperty& operator=(const EdgeProperty& other) {

    if(this != &other) {
      this->m_obj->removeEdgeProperty(this);

      this->m_obj = other.m_obj;
      this->m_data = other.m_data;

      this->m_obj->registerEdgeProperty(this);
    }

    return *this;
  }
};

template <class T>
class FaceProperty : public TopObjProperty<T> {
public:
  FaceProperty(TopologicalObject* obj) : TopObjProperty<T>(obj,obj->numFaceSlots()) {
    this->m_obj->registerFaceProperty(this);
  }

  explicit FaceProperty(const FaceProperty& prop) : TopObjProperty<T>(prop.m_obj,prop.m_obj->numFaceSlots()) {
    this->m_obj->registerFaceProperty(this);
    this->m_data = prop.m_data;
  }

  ~FaceProperty() {
    this->m_obj->removeFaceProperty(this);
  }

  FaceProperty& operator=(const FaceProperty& other) {

    if(this != &other) {
      this->m_obj->removeFaceProperty(this);

      this->m_obj = other.m_obj;
      this->m_data = other.m_data;

      this->m_obj->registerFaceProperty(this);
    }

    return *this;
  }
};


template <class T>
class TetProperty : public TopObjProperty<T> {
public:
  TetProperty(TopologicalObject* obj) : TopObjProperty<T>(obj, obj->numTetSlots()) {
    this->m_obj->registerTetProperty(this);
  }

  explicit TetProperty(const TetProperty& prop) : TopObjProperty<T>(prop.m_obj,prop.m_obj->numTetSlots()) {
    this->m_obj->registerTetProperty(this);
    this->m_data = prop.m_data;
  }

  ~TetProperty() {
    this->m_obj->removeTetProperty(this);
  }

  TetProperty& operator=(const TetProperty& other) {

    if(this != &other) {
      this->m_obj->removeTetProperty(this);

      this->m_obj = other.m_obj;
      this->m_data = other.m_data;

      this->m_obj->registerTetProperty(this);
    }

    return *this;
  }
  
};

}

#endif
