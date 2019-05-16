// Copyright (c) 2012-2019 University of Lyon and CNRS (France).
// All rights reserved.
//
// This file is part of MEPP2; you can redistribute it and/or modify
// it under the terms of the GNU Lesser General Public License as
// published by the Free Software Foundation; either version 3 of
// the License, or (at your option) any later version.
//
// This file is provided AS IS with NO WARRANTY OF ANY KIND, INCLUDING THE
// WARRANTY OF DESIGN, MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE.
#pragma once

#include <vector>
#include <limits>
#include <stdexcept>
#include <cmath>

namespace FEVV {
namespace DataStructures {
namespace AIF {
template< typename CoordinateT = double, unsigned short DIM = 3 >
class AIFVector;
/*
 * The class to store vertex coordinates.
 * It extends std::array<...,3>
 * with a constructor that takes 3 coordinates.
 * See FEVV/Filters/Generic/scaling.hpp usecase.
 */
template< typename CoordinateT = double, unsigned short DIM = 3 >
class AIFPoint : public std::array< CoordinateT, DIM >
{
public:
  typedef CoordinateT                       CoordinateType;
  typedef std::array< CoordinateType, DIM > SuperClass;
  // TODO-elo refactor all function to handle DIM!=3 case !

  /*!
   *			AIFPoint default constructor.
   */
  AIFPoint(void)
  {
    if(DIM != 3)
      throw std::runtime_error("AIFPoint only support dim=3 for now.");

    (*this)[0] = 0;
    (*this)[1] = 0;
    (*this)[2] = 0;
  }

  /*!
   *			AIFPoint constructor with 3 coordinates.
   */
  AIFPoint(CoordinateType x, CoordinateType y, CoordinateType z)
  {
    if(DIM != 3)
      throw std::runtime_error("AIFPoint only support dim=3 for now.");

    (*this)[0] = x;
    (*this)[1] = y;
    (*this)[2] = z;
  }

  /*!
   * 			+ operator
   * \param	other	The second AIFPoint to add to the current one.
   * \return The resulting AIFPoint object.
   */
  const AIFPoint< CoordinateType, DIM >
  operator+(const AIFPoint< CoordinateType, DIM > &other) const
  {
    return AIFPoint< CoordinateType, DIM >(
        (*this)[0] + other[0], (*this)[1] + other[1], (*this)[2] + other[2]);
  }

  /*!
   * 			<< operator
   * \param	stream	The output stream to update.
   * \param	p	The AIFPoint to add to the output stream.
   * \return The resulting output stream object.
   */
  friend std::ostream &operator<<(std::ostream &stream,
                                  const AIFPoint< CoordinateType, DIM > &p)
  {
    // stream << "("; // do not add this to get the same output as CGAL or
    for(auto itVal = p.begin(); itVal != p.end(); ++itVal)
    {
      stream << (*itVal);
      if(itVal != (p.end() - 1))
        stream << " ";
    }
    // stream << ")";

    return stream;
  }

  /*!
   * 			- operator
   * \param	other	The second AIFPoint to substract to the current one.
   * \return The resulting AIFVector object.
   */
  AIFVector< CoordinateType, DIM >
  operator-(const AIFPoint< CoordinateType, DIM > &other) const
  {
    return AIFVector< CoordinateType, DIM >(
        (*this)[0] - other[0], (*this)[1] - other[1], (*this)[2] - other[2]);
  }
  /*!
   * 			+ operator
   * \param	v	The AIFVector to add to the current AIFPoint.
   * \return The resulting AIFPoint object.
   */
  AIFPoint< CoordinateType, DIM >
  operator+(const AIFVector< CoordinateType, DIM > &v) const
  {
    return AIFPoint< CoordinateType, DIM >(
        (*this)[0] + v[0], (*this)[1] + v[1], (*this)[2] + v[2]);
  }

  /*!
   * 			< operator
   * \param	other	The second AIFPoint to compare to the current AIFPoint.
   * \return The lexicographical compare.
   */
  bool operator<(const AIFPoint< CoordinateType, DIM > &p) const
  {
    if((*this)[0] < p[0])
      return true;
    else if((*this)[0] > p[0])
      return false;

    if((*this)[1] < p[1])
      return true;
    else if((*this)[1] > p[1])
      return false;

    if((*this)[2] < p[2])
      return true;
    else
      return false;
  }
  
  /*!
   * 			== operator
   * \param	other	The second AIFPoint to compare to the current AIFPoint.
   * \return The equality compare.
   */
  bool operator==(const AIFPoint< CoordinateType, DIM > &p) const
  {
    if( this == &p )
      return true;
  
    if( fabs((*this)[0] - p[0]) > std::numeric_limits<CoordinateType>::epsilon() )
      return false;
    
	if( fabs((*this)[1] - p[1]) > std::numeric_limits<CoordinateType>::epsilon() )
      return false;

	if( fabs((*this)[2] - p[2]) > std::numeric_limits<CoordinateType>::epsilon() )
      return false;
  
    return true;
  }  
};


/*
 * A class to store a vector (usefull for normals).
 * It extends std::array<...,3>
 * with a constructor that takes 3 coordinates.
 */
template< typename CoordinateT /*=double*/, unsigned short DIM /*=3*/ >
class AIFVector : public std::array< CoordinateT, DIM >
{
public:
  typedef CoordinateT                       CoordinateType;
  // TODO-elo refactor all function to handle DIM!=3 case !

  /*!
   *			AIFVector default constructor.
   */
  AIFVector(void)
  {
    if(DIM != 3)
      throw std::runtime_error("AIFVector only support dim=3 for now.");

    (*this)[0] = 0;
    (*this)[1] = 0;
    (*this)[2] = 0;
  }

  /*!
   *			AIFVector constructor with 3 coordinates.
   */
  AIFVector(CoordinateType x, CoordinateType y, CoordinateType z)
  {
    if(DIM != 3)
      throw std::runtime_error("AIFVector only support dim=3 for now.");

    (*this)[0] = x;
    (*this)[1] = y;
    (*this)[2] = z;
  }

  /*!
   * 			The Euclidean length.
   * \return The current AIFVector object Euclidean length.
   */
  CoordinateType length(void) const
  {
    CoordinateType tmp = 0;
    for(auto v : *this)
      tmp += v * v;
    return sqrt(tmp);
  }

  /*!
   * 			+ operator
   * \param	other	The second AIFVector to add to the current one.
   * \return The resulting AIFVector object.
   */
  AIFVector< CoordinateType, DIM >
  operator+(const AIFVector< CoordinateType, DIM > &other) const
  {
    return AIFVector(
        (*this)[0] + other[0], (*this)[1] + other[1], (*this)[2] + other[2]);
  }
  /*!
   * 			- operator
   * \param	other	The second AIFVector to substract to the current one.
   * \return The resulting AIFVector object.
   */
  AIFVector< CoordinateType, DIM >
  operator-(const AIFVector< CoordinateType, DIM > &other) const
  {
    return AIFVector(
        (*this)[0] - other[0], (*this)[1] - other[1], (*this)[2] - other[2]);
  }
  /*!
   * 			/ operator
   * \param	k	The scalar to divide each coordinate.
   * \return The resulting AIFVector object.
   */
  AIFVector< CoordinateType, DIM > operator/(double k) const
  {
    return AIFVector< CoordinateType, DIM >(
        (*this)[0] / k, (*this)[1] / k, (*this)[2] / k);
  }
  /*!
   * 			* operator
   * \param	other	The second AIFVector to do a dot product with the current
   * one. 
   * \return The resulting dot product (scalar) value.
   */
  double operator*(const AIFVector< CoordinateType, DIM > &other) const
  {
    return ((*this)[0] * other[0] + (*this)[1] * other[1] +
            (*this)[2] * other[2]);
  }
  /*!
   * 			* operator
   * \param	k	The scalar to multiply each coordinate.
   * \return The resulting AIFVector object.
   */
  AIFVector< CoordinateType, DIM > operator*(double k) const
  {
    return AIFVector< CoordinateType, DIM >(
        (*this)[0] * k, (*this)[1] * k, (*this)[2] * k);
  }
  /*!
   * 			+ operator
   * \param	other	The second AIFPoint to add to the current AIFVector.
   * \return The resulting AIFPoint object.
   */
  AIFPoint< CoordinateType, DIM >
  operator+(const AIFPoint< CoordinateType, DIM > &p) const
  {
    return p + *this;
  }

  /*!
   * 			< operator
   * \param	other	The second AIFVector to compare to the current
   * AIFVector. 
   * \return The lexicographical compare.
   */
  bool operator<(const AIFVector< CoordinateType, DIM > &v) const
  {
    if((*this)[0] < v[0])
      return true;
    else if((*this)[0] > v[0])
      return false;

    if((*this)[1] < v[1])
      return true;
    else if((*this)[1] > v[1])
      return false;

    if((*this)[2] < v[2])
      return true;
    else
      return false;
  }
  
  /*!
   * 			== operator
   * \param	other	The second AIFVector to compare to the current AIFVector.
   * \return The equality compare.
   */
  bool operator==(const AIFVector< CoordinateType, DIM > &v) const
  {
    if( this == &v )
      return true;
  
    if( fabs((*this)[0] - v[0]) > std::numeric_limits<CoordinateType>::epsilon() )
      return false;
    
	if( fabs((*this)[1] - v[1]) > std::numeric_limits<CoordinateType>::epsilon() )
      return false;

	if( fabs((*this)[2] - v[2]) > std::numeric_limits<CoordinateType>::epsilon() )
      return false;
  
    return true;
  }    

  /*!
   * 			* operator
   * \param	s	The scalar to multiply each coordinate.
   * \param	v	The AIFVector object whose coordinates are multiplied by
   * s. 
   * \return The resulting AIFVector object.
   */
  friend AIFVector< CoordinateType, DIM >
  operator*(double s, const AIFVector< CoordinateType, DIM > &v)
  {
    return v * s;
  }
  
  /*!
  * 			<< operator
  * \param	stream	The output stream to update.
  * \param	v	The AIFVector to add to the output stream.
  * \return The resulting output stream object.
  */
  friend std::ostream& operator<<(std::ostream& stream, const AIFVector<CoordinateType, DIM>& v)
  {
    //stream << "("; // do not add this to get the same output as CGAL or OpenMesh
    for(auto itVal = v.begin(); itVal != v.end(); ++itVal)
    {
      stream << (*itVal);
      if( itVal != (v.end() - 1) )
        stream << " ";
    }
    //stream << ")";

    return stream;
  }  
};

/*
 * Base class for property maps.
 * To be able to store pointers to property maps of different types
 * in the same container.
 */
class BasePropertyMap
{
public:
  virtual ~BasePropertyMap() {}
  virtual void forcePolymorphic(void) {}
  virtual void remove(std::size_t idx, std::size_t cLastIdx) {}
};

/*
 * The purpose of this class is just to boost::vector_property_map
 * with an operator[] that takes a vertex_descriptor as parameter.
 * See FEVV/Filters/Generic/scaling.hpp usecase.
 */
template< typename T >
class PropertyMap : public boost::vector_property_map< T >,
                    public BasePropertyMap
{
public:
  /// make base class operator[] visible
  using boost::vector_property_map< T >::operator[];

  /// add operator[] for vertex_descriptor
  T &operator[](const AIFVertex::ptr v) const { return (*this)[v->GetIndex()]; }

  /// add operator[] for edge_descriptor
  T &operator[](const AIFEdge::ptr v) const { return (*this)[v->GetIndex()]; }

  /// add operator[] for face_descriptor
  T &operator[](const AIFFace::ptr v) const { return (*this)[v->GetIndex()]; }

  /// number of elements in the property map
  std::size_t size(void) const
  {
    return std::distance(this->storage_begin(), this->storage_end());
    // no access to underlying std::vector :(
  }

  /**
   * (Re)move one element in property map.
   * \param  idx  index of element to remove
   * \param  cLastIdx  index of the last element of the cell container
   */
  void remove(std::size_t idx, std::size_t cLastIdx)
  {
    // Remove element at #idx in the property map.
    // If #idx is not the last element of the cell container,
    // then replace #idx by the last element of the cell container.
    // element, then remove the last element.
    // It's like AIFCellContainer::remove(), except that we can NOT remove
    // the last element because it is not allowed by
    // boost::vector_property_map interface.
    // Hence the need to pass the cLastIdx parameter, because the size of
    // the property map and the size of the cell container can be different.

    if(this->size() == 0 || cLastIdx >= this->size())
      return;

    if(idx != cLastIdx)
      (*this)[idx] = (*this)[cLastIdx];

    // not allowed by boost::vector_property_map
    // this->pop_back();
  }
};


/*
 * An associative property map, needed for example to store texture
 * per HalfEdge.
 * Note1: boost::associative_property_map is an adaptor that converts an
 * associative container such as std::map into property map ; the adaptor only
 * retains a reference to the container, so the lifetime of the container must
 * encompass the use of the adaptor ; see the description of
 * boost::associative_property_map for more details. Note2: the present class
 * completes boost::associative_property_map by also storing the underlying
 * associative container.
 */
template< typename KeyType, typename ValueType >
class AssocPropertyMap
    : public boost::associative_property_map< std::map< KeyType, ValueType > >,
      public BasePropertyMap
{
public:
  AssocPropertyMap(void)
  {
    // at AssocPropertyMap creation, the boost::associative_property_map part
    // of the object does not point to a valid associative container because
    // none has been provided to  boost::associative_property_map constructor;
    // so the following line initializes the boost::associative_property_map
    // part of the object with the container stored in the object
    boost::associative_property_map< std::map< KeyType, ValueType > >::
    operator=(
        boost::associative_property_map< std::map< KeyType, ValueType > >(map));
  }

private:
  std::map< KeyType, ValueType > map; // the real associative container
};


/**
 * A property maps container.
 */
class PropertyMapContainer
{
public:
  ~PropertyMapContainer()
  {
    // delete property maps
    for(auto &m : m_PropertyMaps)
      if(m.second != NULL)
      {
        delete m.second;
        m.second = NULL; // in case of cross-reference
      }
  }

  /**
   * Test if a property map exists.
   *
   * \param  name  name of the property map
   * \return  true if the property map exists, else false
   */
  bool isPropertyMap(const std::string &name) const
  {
    if(m_PropertyMaps.count(name) > 0)
      return true;
    else
      return false;
  }
  /**
   * Test if there is at least one property map starting with prefix.
   *
   * \param  prefix  prefix of the property map
   * \return true if there is at least one property map starting 
   *         with given prefix, else false
   */  
  bool isAPropertyMapStartingWithPrefix(const std::string &prefix) const
  {
    auto iter = m_PropertyMaps.lower_bound(prefix);
	  if( iter != m_PropertyMaps.end() )
	  {
      return (iter->first.compare(0, prefix.size(), prefix) == 0);
	  }
	  return false;
  }

  /**
  * Get property map names starting with prefix.
  *
  * \param  prefix  prefix of the property map
  * \return std::vector of property map names starting with prefix
  */
  std::vector< std::string > GetPropertyMapNamesStartingWithPrefix(const std::string &prefix) const
  {
    std::vector< std::string > res;
    auto iter = m_PropertyMaps.lower_bound(prefix);
    auto iter_end = m_PropertyMaps.end();
    for( ; iter!=iter_end; ++iter)
    {
      if (iter->first.compare(0, prefix.size(), prefix) == 0)
        res.push_back(iter->first);
    }
    return res;
  }

  /**
   * Get the property map if it exists, else throw an exception.
   *
   * \param  name  name of the property map
   * \return  a pointer to the property map
   */
  template< typename T >
  PropertyMap< T > *getPropertyMap(const std::string &name) const
  {
    if(!isPropertyMap(name))
      throw std::runtime_error("Property map '" + name + "' doesn't exist.");

    PropertyMap< T > *pm =
        dynamic_cast< PropertyMap< T > * >(m_PropertyMaps.at(name));
    if(!pm)
      throw std::runtime_error("Property map '" + name +
                               "' is not of the given type.");

    return pm;
  }

  /**
   * Add a property map or return an existing one.
   *
   * \param  name  name of the property map
   * \return  a pointer to the property map
   */
  template< typename T >
  PropertyMap< T > *addPropertyMap(const std::string &name)
  {
    if(!isPropertyMap(name))
    {
      // create property map
      m_PropertyMaps[name] = new PropertyMap< T >;
    }

    return getPropertyMap< T >(name);
  }

  /**
   * Remove a property map.
   *
   * \param  name  name of the property map
   */
  void removePropertyMap(const std::string &name)
  {
    if(!isPropertyMap(name))
      return;

    delete m_PropertyMaps.at(name);
    m_PropertyMaps.erase(name);
  }

  /**
  * Remove all property maps.
  */
  void clear() {
    // delete property maps
    for (auto &m : m_PropertyMaps)
      if (m.second != NULL)
      {
        delete m.second;
        m.second = NULL; // in case of cross-reference
      }
    m_PropertyMaps.clear();
  }

  /**
   * Set value at index in property map.
   *
   * \param  name  name of the property map
   * \param  idx  index where to set the value at
   * \param  value  value at the given index
   */
  template< typename T >
  void setProperty(const std::string &name, std::size_t idx, const T &value)
  {
    PropertyMap< T > *prop_map = getPropertyMap< T >(name);
    (*prop_map)[idx] = value;
  }

  /**
   * Get the value of a property map at the given index.
   *
   * \param  name  name of the property map
   * \param  idx  index where to set the value from
   * \param  value  value at the given index
   */
  template< typename T >
  T &getProperty(const std::string &name, std::size_t idx)
  {
    PropertyMap< T > *prop_map = getPropertyMap< T >(name);
    return (*prop_map)[idx];
  }

  /**
   * Get the value of a property map at the given index.
   *
   * \param  name  name of the property map
   * \param  idx  index where to set the value from
   * \param  value  value at the given index
   */
  template< typename T >
  T &getProperty(const std::string &name, std::size_t idx) const
  {
    const PropertyMap< T > *prop_map = getPropertyMap< T >(name);
    return (*prop_map)[idx];
  }

  /**
   * Remove property at index, in all property maps.
   *
   * \param  idx  the index where to remove properties
   * \param  cLastIdx  index of the last element of the cell container
   */
  void removeProperties(std::size_t idx, std::size_t cLastIdx)
  {
    for(auto it : m_PropertyMaps)
    {
      auto map = it.second;
      map->remove(idx, cLastIdx);
    }
  }

  //------------------Associative prop maps -------------------------------

  /**
   * Get the property map if it exists, else throw an exception.
   *
   * \param  name  name of the property map
   * \return  a pointer to the property map
   */
  template< typename KeyType, typename ValueType >
  AssocPropertyMap< KeyType, ValueType > *
  getAssocPropertyMap(const std::string &name)
  {
    if(!isPropertyMap(name))
      throw std::runtime_error("Property map '" + name + "' doesn't exist.");

    AssocPropertyMap< KeyType, ValueType > *pm =
        dynamic_cast< AssocPropertyMap< KeyType, ValueType > * >(
            m_PropertyMaps.at(name));
    if(!pm)
      throw std::runtime_error("Property map '" + name +
                               "' is not of the given type.");

    return pm;
  }

  /**
   * Add a property map or return an existing one.
   *
   * \param  name  name of the property map
   * \return  a pointer to the property map
   */
  template< typename KeyType, typename ValueType >
  AssocPropertyMap< KeyType, ValueType > *
  addAssocPropertyMap(const std::string &name)
  {
    if(!isPropertyMap(name))
    {
      // create property map
      m_PropertyMaps[name] = new AssocPropertyMap< KeyType, ValueType >();
    }

    return getAssocPropertyMap< KeyType, ValueType >(name);
  }

private:
  std::map< const std::string, BasePropertyMap * > m_PropertyMaps;
};


} // namespace AIF
} // namespace DataStructures
} // namespace FEVV

