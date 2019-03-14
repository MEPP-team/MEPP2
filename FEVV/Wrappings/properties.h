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

#include <CGAL/boost/graph/properties.h>
#include <boost/shared_ptr.hpp>
#include <boost/any.hpp>

#include "FEVV/Types/Material.h"
#include "FEVV/Types/Gui_properties.h"
#include "FEVV/Wrappings/Geometry_traits.h"

/**
 * \file
 * Management of generic property maps.
 * Refer to \ref GenericPropertyMapConceptPage.
 */


namespace FEVV {

//---------------------------------------------------------
// list of the defaulty supported properties

/// the *normal* property of a vertex
/// (refer to \ref GenericPropertyMapConceptPage)
enum vertex_normal_t { vertex_normal };

/// the *tangent* property of a vertex
/// (refer to \ref GenericPropertyMapConceptPage)
enum vertex_tangent_t { vertex_tangent };

/// the *texcoord* property of a vertex
/// (refer to \ref GenericPropertyMapConceptPage)
enum vertex_texcoord_t { vertex_texcoord };

/// the *color* property of a vertex
/// (refer to \ref GenericPropertyMapConceptPage)
enum vertex_color_t { vertex_color };
// warning: boost::vertex_color_t is already
// defined in /usr/include/boost/graph/properties.hpp

/// the *normal* property of a halfedge (aka normal by corner)
/// (refer to \ref GenericPropertyMapConceptPage)
enum halfedge_normal_t { halfedge_normal };

/// the *texcoord* property of a halfedge (aka texcoord by corner)
/// (refer to \ref GenericPropertyMapConceptPage)
enum halfedge_texcoord_t { halfedge_texcoord };

/// the *color* property of an edge
/// (refer to \ref GenericPropertyMapConceptPage)
enum edge_color_t { edge_color };

/// the *normal* property of a face
/// (refer to \ref GenericPropertyMapConceptPage)
enum face_normal_t { face_normal };

/// the *color* property of a face
/// (refer to \ref GenericPropertyMapConceptPage)
enum face_color_t { face_color };

/// the *material* property of a face
/// (refer to \ref GenericPropertyMapConceptPage)
enum face_material_t { face_material };

/// the *materials* property of a mesh
/// (refer to \ref GenericPropertyMapConceptPage)
enum mesh_materials_t { mesh_materials };

/// the *guiproperties* property of a mesh
/// (refer to \ref GenericPropertyMapConceptPage)
enum mesh_guiproperties_t { mesh_guiproperties };


//---------------------------------------------------------
// property map container and related functions

typedef std::map< std::string, boost::any > PMapsContainer;

/// (refer to \ref GenericPropertyMapConceptPage)
inline bool
has_map(const PMapsContainer &pmaps, const std::string &map_name)
{
  return (pmaps.count(map_name) > 0);
}

/// (refer to \ref GenericPropertyMapConceptPage)
template< typename PropertyT >
bool
has_map(const PMapsContainer &pmaps, PropertyT p)
{
  return has_map(pmaps, get_property_map_name(p));
}


//---------------------------------------------------------
// associate a string to each property type in order to store
// the property maps in the container which use a string key

/// (refer to \ref GenericPropertyMapConceptPage)
inline std::string get_property_map_name(FEVV::vertex_normal_t)
{
  return std::string("v:normal");
}

/// (refer to \ref GenericPropertyMapConceptPage)
inline std::string get_property_map_name(FEVV::vertex_tangent_t)
{
  return std::string("v:tangent");
}

/// (refer to \ref GenericPropertyMapConceptPage)
inline std::string get_property_map_name(FEVV::vertex_texcoord_t)
{
  return std::string("v:texcoord");
}

/// (refer to \ref GenericPropertyMapConceptPage)
inline std::string get_property_map_name(FEVV::vertex_color_t)
{
  return std::string("v:color");
}

/// (refer to \ref GenericPropertyMapConceptPage)
inline std::string get_property_map_name(FEVV::halfedge_normal_t)
{
  return std::string("h:normal");
}

/// (refer to \ref GenericPropertyMapConceptPage)
inline std::string get_property_map_name(FEVV::halfedge_texcoord_t)
{
  return std::string("h:texcoord");
}

/// (refer to \ref GenericPropertyMapConceptPage)
inline std::string get_property_map_name(FEVV::edge_color_t)
{
  return std::string("e:color");
}

/// (refer to \ref GenericPropertyMapConceptPage)
inline std::string get_property_map_name(FEVV::face_normal_t)
{
  return std::string("f:normal");
}

/// (refer to \ref GenericPropertyMapConceptPage)
inline std::string get_property_map_name(FEVV::face_color_t)
{
  return std::string("f:color");
}

/// (refer to \ref GenericPropertyMapConceptPage)
inline std::string get_property_map_name(FEVV::face_material_t)
{
  return std::string("f:material");
}

/// (refer to \ref GenericPropertyMapConceptPage)
inline std::string get_property_map_name(FEVV::mesh_materials_t)
{
  return std::string("m:materials");
}

/// (refer to \ref GenericPropertyMapConceptPage)
inline std::string get_property_map_name(FEVV::mesh_guiproperties_t)
{
  return std::string("m:guiproperties");
}


//---------------------------------------------------------
// Associative property map type

/*
 * An associative property map, needed for example for CGAL::Polyhedron
 * which does not manage vertex/halfedge/face indices natively.
 * Notes:
 *    1) boost::associative_property_map is an adaptor that converts an
 * associative container such as std::map into property map ; the adaptor only
 * retains a reference to the container, so the lifetime of the container must
 * encompass the use of the adaptor ; see the description of
 * boost::associative_property_map for more details. 2) the present class
 * completes boost::associative_property_map by also storing the underlying
 * associative container. 3) the present class mimics the behavior of
 * boost::vector_property_map (see
 * http://www.boost.org/doc/libs/1_37_0/boost/vector_property_map.hpp), aka the
 * copy constructor and the assignment operator do a shallow copy of the
 * underlying std::map container ; as a consequence, like with get(vertex_point,
 * g), get(property_t, g) will return a 'shared copy' of the associative
 * property map.
 */
template< typename KeyType, typename ValueType >
class Assoc_property_map
    : public boost::associative_property_map< std::map< KeyType, ValueType > >
{
public:
  // default constructor
  Assoc_property_map(void) : map_ptr(new std::map< KeyType, ValueType >)
  {
    // at Assoc_property_map creation, the boost::associative_property_map part
    // of the object does not point to a valid associative container because
    // none has been provided to  boost::associative_property_map constructor;
    // so the following line initializes the boost::associative_property_map
    // part of the object with the container stored in the Assoc_property_map
    // object
    boost::associative_property_map< std::map< KeyType, ValueType > >::
    operator=(boost::associative_property_map< std::map< KeyType, ValueType > >(
        *map_ptr));
    // DBG std::cout << __FUNCTION__ << " " << this << std::endl;
  }

  // copy constructor (shallow copy !)
  Assoc_property_map(const Assoc_property_map &other)
  {
    // DBG std::cout << __FUNCTION__ << " " << this << " <- copy(" << &other <<
    // ")" << std::endl;
    map_ptr = other.map_ptr;
    boost::associative_property_map< std::map< KeyType, ValueType > >::
    operator=(boost::associative_property_map< std::map< KeyType, ValueType > >(
        *map_ptr));
  }

  // destructor
  ~Assoc_property_map()
  {
    // DBG std::cout << __FUNCTION__ << " " << this << std::endl;
  }

  // assignment operator (shallow copy !)
  Assoc_property_map &operator=(const Assoc_property_map &other)
  {
    // DBG std::cout << __FUNCTION__ << " " << this << " <- operator=(" <<
    // &other << ")" << std::endl;
    map_ptr = other.map_ptr;
    boost::associative_property_map< std::map< KeyType, ValueType > >::
    operator=(boost::associative_property_map< std::map< KeyType, ValueType > >(
        *map_ptr));

    return *this;
  }

  // debugging informations
  void infos(void)
  {
    std::cout << "this=" << this << "  map.size()=" << (*map_ptr).size()
              << std::endl;
  }

private:
  boost::shared_ptr< std::map< KeyType, ValueType > >
      map_ptr; // the real associative container
};


//---------------------------------------------------------
// property map traits and its specializations
//
// The property maps traits make it possible to deduce the
// kind of property map (associative, vector...) to use with
// the mesh type.
// It also make it possible to deduce the full property map type
// from the property type (vertex_normal, vertex_texcoord,
// face_normal...) and the mesh type (CGAL::Polyhedron,
// OpenMesh...).


// define default property map type (associative property map,
// vector property map...) for each cell type (vertex, face...) ;

// default vertex property map
template< typename MeshT, typename ValueT >
struct Vertex_pmap_traits
{
  typedef
      typename boost::graph_traits< MeshT >::vertex_descriptor vertex_key_type;
  typedef typename FEVV::Assoc_property_map< vertex_key_type, ValueT >
      pmap_type;

  static pmap_type create(const MeshT &m)
  {
    pmap_type pmap;
    return pmap;
  }
};

// default face property map
template< typename MeshT, typename ValueT >
struct Face_pmap_traits
{
  typedef typename boost::graph_traits< MeshT >::face_descriptor face_key_type;
  typedef typename FEVV::Assoc_property_map< face_key_type, ValueT > pmap_type;

  static pmap_type create(const MeshT &m)
  {
    pmap_type pmap;
    return pmap;
  }
};

// default edge property map
template< typename MeshT, typename ValueT >
struct Edge_pmap_traits
{
  typedef typename boost::graph_traits< MeshT >::edge_descriptor edge_key_type;
  typedef typename FEVV::Assoc_property_map< edge_key_type, ValueT > pmap_type;

  static pmap_type create(const MeshT &m)
  {
    pmap_type pmap;
    return pmap;
  }
};

// default halfedge property map
template< typename MeshT, typename ValueT >
struct Halfedge_pmap_traits
{
  typedef typename boost::graph_traits< MeshT >::halfedge_descriptor
      halfedge_key_type;
  typedef typename FEVV::Assoc_property_map< halfedge_key_type, ValueT >
      pmap_type;

  static pmap_type create(const MeshT &m)
  {
    pmap_type pmap;
    return pmap;
  }
};


// define standard property map types for (property,cell) pair,
// for example vertex-normal
template< typename MeshT, typename PropertyT >
struct _PMap_traits;

// specialize the property maps traits for vertex-normal
template< typename MeshT >
struct _PMap_traits< MeshT, FEVV::vertex_normal_t >
{
  typedef typename FEVV::Geometry_traits< MeshT >::Vector value_type;
  typedef typename Vertex_pmap_traits< MeshT, value_type >::pmap_type pmap_type;

  static pmap_type create(const MeshT &m)
  {
    pmap_type pmap;
    return pmap;
  }
};

// specialize the property maps traits for vertex-tangent
template< typename MeshT >
struct _PMap_traits< MeshT, FEVV::vertex_tangent_t >
{
  typedef typename FEVV::Geometry_traits< MeshT >::Vector value_type;
  typedef typename Vertex_pmap_traits< MeshT, value_type >::pmap_type pmap_type;

  static pmap_type create(const MeshT &m)
  {
    pmap_type pmap;
    return pmap;
  }
};

// specialize the property maps traits for vertex-texcoord
template< typename MeshT >
struct _PMap_traits< MeshT, FEVV::vertex_texcoord_t >
{
  typedef typename FEVV::Geometry_traits< MeshT >::Vector value_type;
  typedef typename Vertex_pmap_traits< MeshT, value_type >::pmap_type pmap_type;

  static pmap_type create(const MeshT &m)
  {
    pmap_type pmap;
    return pmap;
  }
};

// specialize the standard property map for vertex-color
template< typename MeshT >
struct _PMap_traits< MeshT, FEVV::vertex_color_t >
{
  typedef typename FEVV::Geometry_traits< MeshT >::Vector value_type;
  typedef typename Vertex_pmap_traits< MeshT, value_type >::pmap_type pmap_type;

  static pmap_type create(const MeshT &m)
  {
    pmap_type pmap;
    return pmap;
  }
};

// specialize the property maps traits for halfedge-normal
template< typename MeshT >
struct _PMap_traits< MeshT, FEVV::halfedge_normal_t >
{
  typedef typename FEVV::Geometry_traits< MeshT >::Vector value_type;
  typedef
      typename Halfedge_pmap_traits< MeshT, value_type >::pmap_type pmap_type;

  static pmap_type create(const MeshT &m)
  {
    pmap_type pmap;
    return pmap;
  }
};

// specialize the property maps traits for halfedge-texcoord
template< typename MeshT >
struct _PMap_traits< MeshT, FEVV::halfedge_texcoord_t >
{
  typedef typename FEVV::Geometry_traits< MeshT >::Vector value_type;
  typedef
      typename Halfedge_pmap_traits< MeshT, value_type >::pmap_type pmap_type;

  static pmap_type create(const MeshT &m)
  {
    pmap_type pmap;
    return pmap;
  }
};

// specialize the property maps traits for edge-color
template< typename MeshT >
struct _PMap_traits< MeshT, FEVV::edge_color_t >
{
  typedef typename FEVV::Geometry_traits< MeshT >::Vector value_type;
  typedef typename Edge_pmap_traits< MeshT, value_type >::pmap_type pmap_type;

  static pmap_type create(const MeshT &m)
  {
    pmap_type pmap;
    return pmap;
  }
};

// specialize the property maps traits for face-normal
template< typename MeshT >
struct _PMap_traits< MeshT, FEVV::face_normal_t >
{
  typedef typename FEVV::Geometry_traits< MeshT >::Vector value_type;
  typedef typename Face_pmap_traits< MeshT, value_type >::pmap_type pmap_type;

  static pmap_type create(const MeshT &m)
  {
    pmap_type pmap;
    return pmap;
  }
};

// specialize the property maps traits for face-color
template< typename MeshT >
struct _PMap_traits< MeshT, FEVV::face_color_t >
{
  typedef typename FEVV::Geometry_traits< MeshT >::Vector value_type;
  typedef typename Face_pmap_traits< MeshT, value_type >::pmap_type pmap_type;

  static pmap_type create(const MeshT &m)
  {
    pmap_type pmap;
    return pmap;
  }
};

// specialize the property maps traits for face-material
template< typename MeshT >
struct _PMap_traits< MeshT, FEVV::face_material_t >
{
  typedef size_t value_type;
  typedef typename Face_pmap_traits< MeshT, value_type >::pmap_type pmap_type;

  static pmap_type create(const MeshT &m)
  {
    pmap_type pmap;
    return pmap;
  }
};

// specialize the property maps traits for mesh materials
// beware: the mesh materials case is very specific,
//         don't use it as a model
template< typename MeshT >
struct _PMap_traits< MeshT, FEVV::mesh_materials_t >
{
  typedef FEVV::Types::Material value_type;
  typedef typename boost::identity_property_map index_map_type;
  typedef typename boost::vector_property_map< value_type, index_map_type >
      pmap_type;

  static pmap_type create(const MeshT &m)
  {
    pmap_type pmap;
    return pmap;
  }
};

// specialize the property maps traits for mesh gui properties
// beware: this case is very specific, don't use it as a model
template< typename MeshT >
struct _PMap_traits< MeshT, FEVV::mesh_guiproperties_t >
{
  typedef FEVV::Types::GuiProperties value_type;
  typedef typename boost::identity_property_map index_map_type;
  typedef typename boost::vector_property_map< value_type, index_map_type >
      pmap_type;

  static pmap_type create(const MeshT &m)
  {
    pmap_type pmap;
    return pmap;
  }
};


// below is the interface to use property maps

// some shortcuts to make it easier to use property maps

// TODO-elo-fixit->  template<typename MeshT, typename PropertyT>
template< typename PropertyT, typename MeshT >
using PMap_traits = _PMap_traits< MeshT, PropertyT >;

template< typename MeshT, typename ValueT >
using Vertex_pmap = typename Vertex_pmap_traits< MeshT, ValueT >::pmap_type;

template< typename MeshT, typename ValueT >
using Face_pmap = typename Face_pmap_traits< MeshT, ValueT >::pmap_type;

template< typename MeshT, typename ValueT >
using Edge_pmap = typename Edge_pmap_traits< MeshT, ValueT >::pmap_type;

template< typename MeshT, typename ValueT >
using Halfedge_pmap = typename Halfedge_pmap_traits< MeshT, ValueT >::pmap_type;


//---------------------------------------------------------
//---------------------------------------------------------
// Using the property maps traits, the functions that
// manage the property maps can be written in a generic way.

// --- property maps bag management ---

/// Make (aka build) a new property map for a given property
/// (refer to \ref GenericPropertyMapConceptPage)
template< typename PropertyT, typename MeshT >
typename PMap_traits< PropertyT, MeshT >::pmap_type
make_property_map(PropertyT p, const MeshT &m)
{
#if 0
	// DBG
	std::string map_name = get_property_map_name(p);
	std::cout << "make_property_map<PropertyT, MeshT>" << ", map name = " << map_name << std::endl;
#endif

  return PMap_traits< PropertyT, MeshT >::create(m);
}


/// Retrieve a property map from the property maps container
/// (refer to \ref GenericPropertyMapConceptPage)
template< typename PropertyT, typename MeshT >
typename PMap_traits< PropertyT, MeshT >::pmap_type
get_property_map(PropertyT p, const MeshT &, const PMapsContainer &pmaps)
{
#if 0
	// DBG
	std::cout << "get_property_map<PropertyT, MeshT>" << ", map name = " << get_property_map_name(p) << std::endl;
#endif

  std::string map_name = get_property_map_name(p);

  return boost::any_cast< typename PMap_traits< PropertyT, MeshT >::pmap_type >(
      pmaps.at(map_name));
}


/// Store a property map in the property maps container
/// (refer to \ref GenericPropertyMapConceptPage)
template< typename PropertyT, typename MeshT >
void
put_property_map(
    PropertyT p,
    const MeshT &,
    PMapsContainer &pmaps,
    const typename PMap_traits< PropertyT, MeshT >::pmap_type &pmap)
{
#if 0
	// DBG
	std::cout << "put_property_map<PropertyT, MeshT>" << ", map name = " << get_property_map_name(p) << std::endl;
#endif

  std::string map_name = get_property_map_name(p);

  pmaps[map_name] = pmap;
}


/// Remove a property map from the property maps container by name,
/// if it exists
/// (refer to \ref GenericPropertyMapConceptPage)
inline
void
remove_property_map_by_name(const std::string &name, PMapsContainer &pmaps)
{
  pmaps.erase(name);
}


/// Remove a property map from the property maps container by type,
/// if it exists
/// (refer to \ref GenericPropertyMapConceptPage)
template< typename PropertyT >
void
remove_property_map(PropertyT p, PMapsContainer &pmaps)
{
  std::string map_name = get_property_map_name(p);
  remove_property_map_by_name(map_name, pmaps);
}


/// List the property maps from the property maps container
/// (refer to \ref GenericPropertyMapConceptPage)
inline
std::vector< std::string >
list_property_maps(const PMapsContainer &pmaps)
{
  std::vector< std::string > names;
  for(auto it = pmaps.begin(); it != pmaps.end(); ++it)
    names.push_back(it->first);

  return names;
}


/// Print the list of the property maps from the property maps container
/// (refer to \ref GenericPropertyMapConceptPage)
inline
void
print_property_maps(const PMapsContainer &pmaps)
{
  std::vector< std::string > names = list_property_maps(pmaps);
  for(auto it = names.begin(); it != names.end(); ++it)
    std::cout << *it << std::endl;
}


// --- easy property maps creation ---

/// Make a new vertex property map to store data of type ValueT
/// (refer to \ref GenericPropertyMapConceptPage)
template< typename MeshT, typename ValueT >
typename Vertex_pmap_traits< MeshT, ValueT >::pmap_type
make_vertex_property_map(const MeshT &m)
{
  return Vertex_pmap_traits< MeshT, ValueT >::create(m);
}


/// Make a new vertex property map to store data of type ValueT
/// (refer to \ref GenericPropertyMapConceptPage)
template< typename MeshT, typename ValueT >
typename Face_pmap_traits< MeshT, ValueT >::pmap_type
make_face_property_map(const MeshT &m)
{
  return Face_pmap_traits< MeshT, ValueT >::create(m);
}


/// Make a new vertex property map to store data of type ValueT
/// (refer to \ref GenericPropertyMapConceptPage)
template< typename MeshT, typename ValueT >
typename Edge_pmap_traits< MeshT, ValueT >::pmap_type
make_edge_property_map(const MeshT &m)
{
  return Edge_pmap_traits< MeshT, ValueT >::create(m);
}


/// Make a new vertex property map to store data of type ValueT
/// (refer to \ref GenericPropertyMapConceptPage)
template< typename MeshT, typename ValueT >
typename Halfedge_pmap_traits< MeshT, ValueT >::pmap_type
make_halfedge_property_map(const MeshT &m)
{
  return Halfedge_pmap_traits< MeshT, ValueT >::create(m);
}


} // namespace FEVV

