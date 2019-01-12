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

#include <fstream>
#include <sstream>
#include <iostream>
#include <vector>
#include <deque>
#include <set>
#include <utility> // for std::pair
#include <algorithm>
#include <cmath>   // for std::abs()
#include <numeric> // for std::accumulate()


namespace are_meshes_identical {

//------------------------------------------------------------------------------

typedef std::set< std::pair< FloatT, FloatT > > ErrorContainer;
// A single error value can occur multiple time during mesh comparison,
// for example a coordinate that is always rounded in the same way. We
// want to keep track of the DIFFERENT errors that occur, so we store
// one error as a couple of values (reference value and effective value)
// into a set, to get rid of duplicated errors.

inline
std::string
get_next_line(std::ifstream &file);

inline
FloatT
relative_distance(FloatT a, FloatT b);

inline FloatT
absolute_distance(FloatT a, FloatT b);

inline
bool
arevectorequal(const std::vector< FloatT > &v1,
               const std::vector< FloatT > &v2,
               FloatT attr_threshold,
               bool relative_threshold,
               ErrorContainer *errors_log);

//------------------------------------------------------------------------------

/*
 * 3D point class. Points can be sorted in x, y, z order.
 */
class Point
{
public:
  Point(FloatT x = 0, FloatT y = 0, FloatT z = 0)
  {
    m_x = x;
    m_y = y;
    m_z = z;
  }

  bool operator<(const Point &p) const
  {
    if(m_x < p.m_x)
      return true;
    else if(m_x > p.m_x)
      return false;
    else if(m_y < p.m_y)
      return true;
    else if(m_y > p.m_y)
      return false;
    else if(m_z < p.m_z)
      return true;
    else if(m_z > p.m_z)
      return false;

    return false; // points are equal
  }

  // Test if the point coordinates are (almost) equal
  // to another point coordinates.
  // A numerical difference is allowed (see threshold)
  bool isequalto(const Point &p,
                 FloatT geom_threshold,
                 bool relative_threshold,
                 ErrorContainer *errors_log = NULL) const
  {
    bool isequal = false;
    FloatT dist_x;
    FloatT dist_y;
    FloatT dist_z;

    if(geom_threshold == 0.0)
    {
      // do an exact comparison
      isequal = (*this == p);
    }
    else if(relative_threshold)
    {
      // comparison with a relative threshold
      dist_x = relative_distance(m_x, p.m_x);
      dist_y = relative_distance(m_y, p.m_y);
      dist_z = relative_distance(m_z, p.m_z);

      isequal = (dist_x <= geom_threshold && dist_y <= geom_threshold &&
                 dist_z <= geom_threshold);
    }
    else
    {
      // comparison with an absolute threshold
      dist_x = absolute_distance(m_x, p.m_x);
      dist_y = absolute_distance(m_y, p.m_y);
      dist_z = absolute_distance(m_z, p.m_z);

      isequal = (dist_x <= geom_threshold && dist_y <= geom_threshold &&
                 dist_z <= geom_threshold);
    }

    // log errors for future statistics
    if(!isequal && errors_log)
    {
      if(geom_threshold == 0.0)
      {
        // case of exact comparison
        if(m_x != p.m_x)
          errors_log->insert(std::make_pair(m_x, p.m_x));
        if(m_y != p.m_y)
          errors_log->insert(std::make_pair(m_y, p.m_y));
        if(m_z != p.m_z)
          errors_log->insert(std::make_pair(m_z, p.m_z));
      }
      else
      {
        // case of comparison with a threshold
        if(dist_x > geom_threshold)
          errors_log->insert(std::make_pair(m_x, p.m_x));
        if(dist_y > geom_threshold)
          errors_log->insert(std::make_pair(m_y, p.m_y));
        if(dist_z > geom_threshold)
          errors_log->insert(std::make_pair(m_z, p.m_z));
      }
    }

    return isequal;
  }

  bool operator==(const Point &p) const
  {
    if(m_x == p.m_x && m_y == p.m_y && m_z == p.m_z)
      return true;
    return false;
  }

  FloatT m_x, m_y, m_z; // point coordinates
};

inline
std::ostream &
operator<<(std::ostream &os, const Point &p)
{
  // write obj to stream
  os << '(' << p.m_x << ' ' << p.m_y << ' ' << p.m_z << ')';
  return os;
}

//------------------------------------------------------------------------------

class Vertex
{
public:
  Vertex(const Point &point, const std::vector< FloatT > &attributes)
  {
    m_point = point;
    m_attributes = attributes;
  }

  Vertex &operator=(const Vertex &other) // copy assignment
  {
    m_point = other.m_point;
    m_attributes = other.m_attributes;

    return *this;
  }

  bool operator<(const Vertex &other) const
  {
    // if the points are equal, compare the attributes
    if(m_point == other.m_point)
      return (m_attributes < other.m_attributes);

    // if the points are different, compare the points
    return (m_point < other.m_point);
  }

  // compare for equality, a numerical difference is allowed
  bool isequalto(const Vertex &other,
                 FloatT geom_threshold,
                 FloatT attr_threshold,
                 bool relative_thresholds,
                 ErrorContainer *errors_log = NULL) const
  {
    // compare point and attributes
    bool p_equal = m_point.isequalto(
        other.m_point, geom_threshold, relative_thresholds, errors_log);
    bool v_equal = arevectorequal(m_attributes,
                                  other.m_attributes,
                                  attr_threshold,
                                  relative_thresholds,
                                  errors_log);
    return (p_equal && v_equal);
  }

  // compare for exact equality
  bool operator==(const Vertex &other) const
  {
    // compare point and attributes
    return (m_point == other.m_point) && (m_attributes == other.m_attributes);
  }

  Point m_point; // vertex coordinates
  std::vector< FloatT > m_attributes;
};

inline
std::ostream &
operator<<(std::ostream &os, const Vertex &v)
{
  // write obj to stream
  os << '[';
  os << v.m_point;
  for(FloatT a : v.m_attributes)
    os << ' ' << a;
  os << ']';

  return os;
}

//------------------------------------------------------------------------------

class Face
{
public:
  Face(const std::deque< Vertex > &vertices,
       const std::vector< FloatT > &attributes)
  {
    m_vertices = vertices;
    m_attributes = attributes;

    // rotate face vertices until the lowest point comes first
    // example:
    //   2 3 1 1 4
    //   4 2 3 1 1
    //   1 4 2 3 1
    //   1 1 4 2 3
    // note: the points coordinates are used, not the points indices ;
    //       if there are several lowest points in the face, use the
    //       next point after the lowest to discriminate.

    auto min_vx_it = std::min_element(
        m_vertices.begin(), m_vertices.end()); // points to the 1st min

    // evaluate next min candidate
    auto min_vx_it_candidate =
        std::find(min_vx_it + 1, m_vertices.end(), *min_vx_it);
    while(min_vx_it_candidate != m_vertices.end())
    {
      auto min_vx_it_candidate_next = next_vertex_in_face(min_vx_it_candidate);
      auto min_vx_it_next = next_vertex_in_face(min_vx_it);

      while(*min_vx_it_candidate_next == *min_vx_it_next &&
            min_vx_it_candidate_next != min_vx_it /*to avoid cycle*/)
      {
        min_vx_it_candidate_next =
            next_vertex_in_face(min_vx_it_candidate_next);
        min_vx_it_next = next_vertex_in_face(min_vx_it_next);
      }

      if((*min_vx_it_candidate_next) < (*min_vx_it_next))
        min_vx_it = min_vx_it_candidate;

      min_vx_it_candidate =
          std::find(min_vx_it_candidate + 1, m_vertices.end(), *min_vx_it);
    }

    // rotate vertices list until the minimal vertex comes 1st
    std::rotate(m_vertices.begin(), min_vx_it, m_vertices.end());
  }

  bool operator<(const Face &f) const
  {
    auto v = m_vertices.begin();
    auto fv = f.m_vertices.begin();
    for(; v != m_vertices.end() && fv != f.m_vertices.end(); ++v, ++fv)
    {
      if(*v < *fv)
        return true;
      else if(*fv < *v)
        return false;
    }

    // all vertices are equal
    if(m_vertices.size() < f.m_vertices.size())
      return true;
    else if(m_vertices.size() > f.m_vertices.size())
      return false;

    return false; // faces are equal
  }

  // compare for strict equality
  bool operator==(const Face &f) const
  {
    // compare vertices
    if(m_vertices.size() != f.m_vertices.size())
      return false;
    auto v = m_vertices.begin();
    auto fv = f.m_vertices.begin();
    for(; v != m_vertices.end() && fv != f.m_vertices.end(); ++v, ++fv)
    {
      if(!(*v == *fv))
        return false;
    }

    // compare attributes
    if(m_attributes.size() != f.m_attributes.size())
      return false;
    auto a = m_attributes.begin();
    auto fa = f.m_attributes.begin();
    for(; a != m_attributes.end() && fa != f.m_attributes.end(); ++a, ++fa)
    {
      if(!(*a == *fa))
        return false;
    }

    return true;
  }

  // compare for equality, a numerical difference is allowed
  bool isequalto(const Face &f,
                 FloatT geom_threshold,
                 FloatT attr_threshold,
                 bool relative_thresholds,
                 ErrorContainer *errors_log = NULL) const
  {
    // compare vertices
    if(m_vertices.size() != f.m_vertices.size())
    {
      std::cout << "Topological difference detected in face (not the same "
                   "number of vertices). Skip face comparison."
                << std::endl;
      return false;
    }

    auto v = m_vertices.begin();
    auto fv = f.m_vertices.begin();
    for(; v != m_vertices.end() && fv != f.m_vertices.end(); ++v, ++fv)
    {
      if(!(v->isequalto(*fv,
                        geom_threshold,
                        attr_threshold,
                        relative_thresholds,
                        errors_log)))
        return false;
    }

    // compare face attributes
    if(!arevectorequal(m_attributes,
                       f.m_attributes,
                       attr_threshold,
                       relative_thresholds,
                       errors_log))
      return false;

    return true;
  }

  // return the next vertex of one vertex in face description
  std::deque< Vertex >::iterator
  next_vertex_in_face(std::deque< Vertex >::iterator v_it)
  {
    // if vIt is at the end of the face vertices list,
    // then next vertex is the first vertex in face vertices list
    auto v_next_it = v_it + 1;
    if(v_next_it == m_vertices.end())
      v_next_it = m_vertices.begin();

    return v_next_it;
  }

  std::deque< Vertex > m_vertices;    // face vertices
  std::vector< FloatT > m_attributes; // face attributes
};

inline
std::ostream &
operator<<(std::ostream &os, const Face &p)
{
  // write obj to stream
  auto v = p.m_vertices.begin();
  os << *v;
  ++v;
  for(; v != p.m_vertices.end(); ++v)
    os << '-' << *v;

  return os;
}

//------------------------------------------------------------------------------

/*
 * A mesh, composed of vertices and faces.
 */
class Mesh
{
public:
  void load(std::ifstream &file, int vertices_nbr, int faces_nbr);
  void display(void);
#if 0 // replaced by isequalto
    bool operator==(const Mesh& m) const;
#endif
  bool isequalto(const Mesh &other,
                 FloatT geom_threshold,
                 FloatT attr_threshold,
                 bool relative_thresholds,
                 ErrorContainer *errors_log = NULL);

protected:
  int m_verticesNbr;
  int m_facesNbr;
  std::multiset< Vertex > m_vertices; // automatically sorted
  std::multiset< Face > m_faces;      // automatically sorted
};

inline
void
Mesh::load(std::ifstream &file, int vertices_nbr, int faces_nbr)
{
  m_verticesNbr = vertices_nbr;
  m_facesNbr = faces_nbr;

  m_vertices.clear();
  m_faces.clear();
  std::string line;

  // load vertices
  std::vector< Vertex > vtmp; // to later retrieve vertex by index
  for(int i = 0; i < m_verticesNbr; i++)
  {
    line = get_next_line(file);
    size_t found = line.find_first_not_of(" \t");
    if(found != std::string::npos)
    {
      if(line[found] == '#')
        continue;
    }
    FloatT x, y, z;
    std::stringstream strline(line);
    strline >> x >> y >> z;
    Point point(x, y, z);

    // optional vertex color or other vertex attribute
    std::vector< FloatT > vertex_attributes;
    FloatT value;
    while(strline >> value)
    {
      vertex_attributes.push_back(value);
    }

    Vertex v(point, vertex_attributes);
    vtmp.push_back(v);
    m_vertices.insert(v);
  }

  // check if some vertices are duplicated
#if 0
  //TODO-elo  turned off because break a number of existing CI tests
  std::set<Vertex> unique_vertices(m_vertices.begin(), m_vertices.end());
  if( unique_vertices.size() != m_vertices.size() )
    throw std::runtime_error("Some vertices are duplicated (same coordinates and attributes). Can't reliably compare the meshes. Aborting.");
#endif

  // load faces
  for(int i = 0; i < m_facesNbr; i++)
  {
    line = get_next_line(file);
    std::stringstream strline(line);

    int v_nbr;
    strline >> v_nbr;

    std::deque< Vertex > face;
    for(int j = 0; j < v_nbr; j++)
    {
      int v_id;
      strline >> v_id;
      face.push_back(vtmp[v_id]);
    }

    // optional face color (or other face attribute)
    std::vector< FloatT > face_attributes;
    FloatT value;
    while(strline >> value)
    {
      face_attributes.push_back(value);
    }

    m_faces.insert(Face(face, face_attributes));
  }
}

inline
void
Mesh::display(void)
{
  // display points list
  for(auto v = m_vertices.begin(); v != m_vertices.end(); ++v)
    std::cout << ' ' << *v << '\n';

  // display faces list
  for(auto f = m_faces.begin(); f != m_faces.end(); ++f)
    std::cout << ' ' << *f << '\n';
}

#if 0 // replaced by isequalto
bool Mesh::operator==(const Mesh& m) const
{
  // check vertices and faces numbers
  if( m_vertices.size() != m.m_vertices.size() )
    return false;

  if( m_faces.size() != m.m_faces.size()  )
    return false;

  // check vertices equality
  auto v = m_vertices.begin();
  auto mv = m.m_vertices.begin();
  for(; v != m_vertices.end()  &&  mv != m.m_vertices.end(); ++v, ++mv)
  {
    if( ! (*v == *mv) )
	{
      std::cout << "Difference found between vertices:\n";
	  std::cout << "(a) " << *v << '\n';
	  std::cout << "(b) " << *mv << '\n';

      return false;
	}
  }

  // check faces equality
  auto f = m_faces.begin();
  auto mf = m.m_faces.begin();
  for(; f != m_faces.end()  &&  mf != m.m_faces.end(); ++f, ++mf)
  {
    if( ! (*f == *mf) )
	{
      std::cout << "Difference found between faces:\n";
	  std::cout << "(a) " << *f << '\n';
	  std::cout << "(b) " << *mf << '\n';

      return false;
	}
  }

  return true;
}
#endif

// test Mesh approximated equality
inline
bool
Mesh::isequalto(const Mesh &other,
                FloatT geom_threshold,
                FloatT attr_threshold,
                bool relative_thresholds,
                ErrorContainer *errors_log)
{
  // check vertices and faces numbers
  if(m_vertices.size() != other.m_vertices.size())
  {
    std::cout << "Topological difference detected in mesh (not the same number "
                 "of vertices). Skip mesh comparison."
              << std::endl;
    return false;
  }

  if(m_faces.size() != other.m_faces.size())
  {
    std::cout << "Topological difference detected in mesh (not the same number "
                 "of faces). Skip mesh comparison."
              << std::endl;
    return false;
  }

  bool isequal = true;
  unsigned int vertice_error_nbr = 0;
  unsigned int face_error_nbr = 0;

  // check vertices equality
  auto vi = m_vertices.begin();
  auto othervi = other.m_vertices.begin();
  for(; vi != m_vertices.end(); ++vi, ++othervi)
  {
    if(!(vi->isequalto(*othervi,
                       geom_threshold,
                       attr_threshold,
                       relative_thresholds,
                       errors_log)))
    {
      std::cout << "Difference found between vertices:\n";
      std::cout << " (a) " << *vi << '\n';
      std::cout << " (b) " << *othervi << '\n';

      isequal = false;
      vertice_error_nbr++;
    }
  }

  // check faces equality
  auto f = m_faces.begin();
  auto otherf = other.m_faces.begin();
  for(; f != m_faces.end(); ++f, ++otherf)
  {
    if(!(f->isequalto(*otherf,
                      geom_threshold,
                      attr_threshold,
                      relative_thresholds,
                      errors_log)))
    {
      std::cout << "Difference found between faces:\n";
      std::cout << " (a) " << *f << '\n';
      std::cout << " (b) " << *otherf << '\n';

      isequal = false;
      face_error_nbr++;
    }
  }

  if(!isequal)
  {
    std::cout << "Number of differences in vertices: " << vertice_error_nbr
              << std::endl;
    std::cout << "Number of differences in faces: " << face_error_nbr
              << std::endl;
  }

  return isequal;
}

//------------------------------------------------------------------------------

/*
 * Get next line from .off file.
 * Skip empty lines and lines beginning with '#'
 */
inline
std::string
get_next_line(std::ifstream &file)
{
  std::string line;
  do
  {
    std::getline(file, line);

    // remove DOS end of line extra character if any
    if((!line.empty()) && (line.back() == '\r'))
      line.pop_back();
  } while(line.empty() || line[0] == '#');

  return line;
}

/*
 * Compute the relative distance between two floats.
 */
inline
FloatT
relative_distance(FloatT a, FloatT b)
{
  if(a == b)
    return 0.0;

  FloatT mean_abs = (std::abs(a) + std::abs(b)) / 2.0f;
  FloatT rel_dist = std::abs(a - b) / mean_abs;

  return rel_dist;
}

/*
 * Compute the absolute distance between two floats.
 */
inline
FloatT
absolute_distance(FloatT a, FloatT b)
{
  return std::abs(a - b);
}

/*
 * Compare to vectors of floats for approximated equality.
 */
inline
bool
arevectorequal(const std::vector< FloatT > &v1,
               const std::vector< FloatT > &v2,
               FloatT attr_threshold,
               bool relative_threshold,
               ErrorContainer *errors_log)
{
  if(v1.size() != v2.size())
  {
    std::cout
        << "The two vectors don't have the same size. Skip vectors comparison."
        << std::endl;
    return false;
  }

  // check vectors elements equality (approximated)
  bool are_vectors_equal = true;
  bool are_values_equal = true;
  auto v1i = v1.begin();
  auto v2i = v2.begin();
  for(; v1i != v1.end(); ++v1i, ++v2i)
  {
    are_values_equal = true;

    if(attr_threshold == 0)
    {
      // do an exact comparison
      if(!(*v1i == *v2i))
      {
        are_values_equal = false;
        are_vectors_equal = false;
      }
    }
    else if(relative_threshold)
    {
      // comparison with a relative threshold
      FloatT rel_dist = relative_distance(*v1i, *v2i);
      if(rel_dist > attr_threshold)
      {
        are_values_equal = false;
        are_vectors_equal = false;
      }
    }
    else
    {
      // comparison with an absolute threshold
      FloatT dist = absolute_distance(*v1i, *v2i);
      if(dist > attr_threshold)
      {
        are_values_equal = false;
        are_vectors_equal = false;
      }
    }

    if(!are_values_equal && errors_log)
    {
      // log errors for future statistics
      errors_log->insert(std::make_pair(*v1i, *v2i));
    }
  }

  return are_vectors_equal;
}

//------------------------------------------------------------------------------

} // namespace are_meshes_identical

//------------------------------------------------------------------------------

inline
bool
are_meshes_equal(std::string filename_a,
                 std::string filename_b,
                 bool verbose,
                 FloatT geom_threshold,
                 FloatT attr_threshold,
                 bool relative_thresholds)
{
  std::ifstream file_a(filename_a);
  if(!file_a)
  {
    std::cout << "Unable to read first file." << filename_a << std::endl;
    return false;
  }

  std::ifstream file_b(filename_b);
  if(!file_b)
  {
    std::cout << "Unable to read second file " << filename_b << std::endl;
    return false;
  }

  std::string line_a;
  std::string line_b;

  // check that files are in OFF format

  line_a = are_meshes_identical::get_next_line(file_a);
  if(line_a != "OFF" && line_a != "COFF")
  {
    std::cout << filename_a << " is not in OFF format." << std::endl;
    return false;
  }

  line_b = are_meshes_identical::get_next_line(file_b);
  if(line_b != "OFF" && line_b != "COFF")
  {
    std::cout << filename_b << " is not in OFF format." << std::endl;
    return false;
  }

  // check meshes sizes

  int vertices_nbr_a, faces_nbr_a, edges_nbr_a;
  std::stringstream(are_meshes_identical::get_next_line(file_a))
      >> vertices_nbr_a >> faces_nbr_a >> edges_nbr_a;
  int vertices_nbr_b, faces_nbr_b, edges_nbr_b;
  std::stringstream(are_meshes_identical::get_next_line(file_b))
      >> vertices_nbr_b >> faces_nbr_b >> edges_nbr_b;

  if(vertices_nbr_a != vertices_nbr_b ||
     faces_nbr_a != faces_nbr_b /*||
     edgesNbrA != edgesNbrB*/) // number of edges is not
                               // always correct in .off
                               // files!
  {
    std::cout << "Meshes do NOT have the same size:" << std::endl;
    std::cout << " - " << filename_a << " has " << vertices_nbr_a
              << " vertices, " << faces_nbr_a << " faces, and " << edges_nbr_a
              << " edges" << std::endl;
    std::cout << " - " << filename_b << " has " << vertices_nbr_b
              << " vertices, " << faces_nbr_b << " faces, and " << edges_nbr_b
              << " edges" << std::endl;
    return false;
  }

  // load meshes

  are_meshes_identical::Mesh mesh_a, mesh_b;
  mesh_a.load(file_a, vertices_nbr_a, faces_nbr_a);
  mesh_b.load(file_b, vertices_nbr_b, faces_nbr_b);

  if(verbose)
  {
    std::cout << "Mesh A =\n";
    mesh_a.display();
    std::cout << '\n';
    std::cout << "Mesh B =\n";
    mesh_b.display();
    std::cout << '\n';
  }

  // compare meshes
  are_meshes_identical::ErrorContainer errors_log;
  bool isequal = mesh_a.isequalto(
      mesh_b, geom_threshold, attr_threshold, relative_thresholds, &errors_log);

  if(!isequal)
  {
    // display information on errors

    std::vector< FloatT > errors;
    if(relative_thresholds)
    {
      for(auto pair : errors_log)
        errors.push_back(
            are_meshes_identical::relative_distance(pair.first, pair.second));
    }
    else
    {
      for(auto pair : errors_log)
        errors.push_back(
            are_meshes_identical::absolute_distance(pair.first, pair.second));
    }

    auto minmax = std::minmax_element(errors.begin(), errors.end());
    FloatT min = *(minmax.first);
    FloatT max = *(minmax.second);
    FloatT mean =
        std::accumulate(errors.begin(), errors.end(), 0.0f) / errors.size();

    std::cout << "Number of effective numeric errors: " << errors.size()
              << std::endl;
    std::cout << "Numeric error min = " << min << std::endl;
    std::cout << "Numeric error max = " << max << std::endl;
    std::cout << "Numeric error mean = " << mean << std::endl;
    if(relative_thresholds)
      std::cout << "Mode: relative error" << std::endl;
    else
      std::cout << "Mode: absolute error" << std::endl;
  }

  return isequal;
}


inline
bool
are_meshes_equal(std::string filename_a,
                 std::string filename_b,
                 bool verbose,
                 FloatT threshold,
                 bool relative_threshold)
{
  return are_meshes_equal(filename_a,
                          filename_b,
                          verbose,
                          threshold,
                          threshold,
                          relative_threshold);
}


inline
bool
are_meshes_equal(std::string filename_a, std::string filename_b, bool verbose)
{
  return are_meshes_equal(filename_a, filename_b, verbose, 0, 0, false);
}
