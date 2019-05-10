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
#include "FEVV/DataStructures/AIF/AIFTopologyHelpers.h" // for copy constructor
#include "FEVV/DataStructures/AIF/AIFPropertiesHelpers.h" // for copy constructor

namespace FEVV {
namespace DataStructures {
namespace AIF {

inline
AIFMesh::AIFMesh(void)
{
  AddPropertyMap< AIFVertex::ptr, Point >(
      "v:point"); // create vertex coordinates property map (mandatory, else
                  // cannot set vertex positions)
}

inline
AIFMesh::AIFMesh(const Self &other)
{
	this->operator=(other);
}

inline 
AIFMesh::~AIFMesh()
{
  this->clear();
}

inline
AIFMesh::Self &
AIFMesh::operator=(const Self &other)
{
  if(this==&other)
    return *this;
  ///////////////////////////////////////////////////////////////////////////
  this->clear();
  ///////////////////////////////////////////////////////////////////////////
  typedef AIFTopologyHelpers helpers;
  typedef AIFPropertiesHelpers PropHelpers;
  typedef helpers::vertex_container::const_iterator vertex_iterator;
  typedef boost::iterator_range< vertex_iterator > vertex_range;
  typedef helpers::edge_container::const_iterator edge_iterator;
  typedef boost::iterator_range< edge_iterator > edge_range;
  typedef helpers::face_container::const_iterator face_iterator;
  typedef boost::iterator_range< face_iterator > face_range;
  typedef helpers::vertex_container_in_face::const_iterator
      vertex_face_iterator;
  typedef boost::iterator_range< vertex_face_iterator > vertex_face_range;

  typedef AIFMesh::CoordinateType coord_type;        // position coordinate type
  typedef AIFMesh::NormalCoordinateType coordN_type; // normal coordinate type
  typedef AIFMesh::CoordinateType coordC_type;       // color coordinate type
  typedef AIFMesh::CoordinateType coordT_type;       // texture coordinate type
  typedef unsigned long index_type;

  std::vector< std::vector< coord_type > > points_coords;
  std::vector< std::vector< coordN_type > > normals_coords;
  std::vector< std::vector< coordT_type > > texture_coords;
  std::vector< std::vector< coordC_type > > vertex_color_coords;
  std::vector< std::vector< coordC_type > > face_color_coords;
  std::vector< std::vector< index_type > > lines_indices, faces_indices;
  std::vector< std::vector< index_type > > texture_face_indices;
  std::vector< std::vector< index_type > > normal_face_indices;
  std::vector< std::vector< coordC_type > > line_color_coords;
  std::vector< std::vector< std::vector< double > > > field_attributes ;
  std::vector< std::string > field_names ;
  // std::vector< Material > materials;
  std::string texture_file_name;

  bool useVertexNormal = false;
  bool useVertexColor = false;
  bool useVertexTextureCoord = false;
  bool useVertexDatafield = false;
  bool useFaceDatafield = false;
  
  AddPropertyMap< AIFVertex::ptr, Point >(
      "v:point"); // create vertex coordinates property map (mandatory, else
                  // cannot set vertex positions)

  if(other.isPropertyMap< AIFVertex::ptr >("v:normal"))
    useVertexNormal = true;
  if(other.isPropertyMap< AIFVertex::ptr >("v:color"))
    useVertexColor = true;
  if(other.isPropertyMap< AIFVertex::ptr >("v:texcoord"))
    useVertexTextureCoord = true;

  bool useLineColor = false;
  if(other.isPropertyMap< AIFEdge::ptr >("e:color"))
    useLineColor = true;

  bool useFaceColor = false;
  if(other.isPropertyMap< AIFFace::ptr >("f:color"))
    useFaceColor = true;
	
  if (other.isAPropertyMapStartingWithPrefix< AIFVertex::ptr >("v:datafield:"))
    useVertexDatafield = true;
  if (other.isAPropertyMapStartingWithPrefix< AIFFace::ptr >("f:datafield:"))
    useFaceDatafield = true;
  
  std::vector<std::string> dfv_names = other.GetPropertyMapNamesStartingWithPrefix< AIFVertex::ptr >("v:datafield:");
  std::vector<std::string> dff_names = other.GetPropertyMapNamesStartingWithPrefix< AIFFace::ptr >("f:datafield:");
  
  unsigned int pointDim = 3;
  long vertexIndex = 0; // to be sure to start to 0
  std::map< helpers::vertex_descriptor, long > indexMap;

  // VERTICES IN CONTAINER
  vertex_range vRange = helpers::vertices(const_cast< Self * >(&other));

  for(vertex_iterator itV = vRange.begin(); itV != vRange.end(); ++itV)
  {
    PropHelpers::Point p =
        PropHelpers::get_point(const_cast< Self * >(&other), *itV);
    std::vector< coord_type > point;

    for(unsigned int i = 0; i < p.size(); ++i)
      point.push_back(p[i]);

    points_coords.push_back(point);

    indexMap[*itV] = vertexIndex;
    ++vertexIndex;

    // deal with vertex attributes here

    if(useVertexNormal)
    {
      // store vertex normal in standard container
      const AIFMesh::Normal &vn =
          other.GetProperty< AIFVertex::ptr, AIFMesh::Normal >(
              "v:normal", (*itV)->GetIndex());
      std::vector< coordN_type > normal(vn.cbegin(), vn.cend());
      normals_coords.push_back(normal);
    }

    if(useVertexColor)
    {
      // store vertex color in standard container
      const AIFMesh::Vector &vc =
          other.GetProperty< AIFVertex::ptr, AIFMesh::Vector >(
              "v:color", (*itV)->GetIndex());
      std::vector< coordC_type > color(vc.cbegin(), vc.cend());
      vertex_color_coords.push_back(color);
    }

    if(useVertexTextureCoord)
    {
      // store vertex texture-coordinates in standard container
      const AIFMesh::PointUV &uv =
          other.GetProperty< AIFVertex::ptr, AIFMesh::PointUV >(
              "v:texcoord", (*itV)->GetIndex());
      std::vector< coordT_type > vt(uv.cbegin(), uv.cend());
      texture_coords.push_back(vt);
    }
	
    if (useVertexDatafield)
    {
      auto it_s = dfv_names.begin(), it_se = dfv_names.end();
      size_t i = 0;
      for (; it_s != it_se; ++it_s, ++i)
      {
        const std::vector< double >&vdata =
          other.GetProperty< AIFVertex::ptr, std::vector< double > >(
            *it_s, (*itV)->GetIndex());

        if (itV == vRange.begin())
        {
          field_names.push_back(*it_s);
          field_attributes.resize(field_attributes.size() + 1);
        }

        field_attributes[i].push_back(vdata);
      }	
	}
  }

  // EDGES IN CONTAINER (needed to take into accound dangling or isolated edges)
  edge_range eRange = helpers::edges(const_cast< Self * >(&other));

  for(edge_iterator itE = eRange.begin(); itE != eRange.end(); ++itE)
  {
    if(helpers::is_dangling_edge(*itE) || helpers::is_isolated_edge(*itE))
    {
      std::vector< index_type > line;

      line.push_back(indexMap[(*itE)->get_first_vertex()]);
      line.push_back(indexMap[(*itE)->get_second_vertex()]);

      lines_indices.push_back(line);

      if(useLineColor)
      {
        // store vertex color in standard container
        const AIFMesh::Vector &vc =
            other.GetProperty< AIFEdge::ptr, AIFMesh::Vector >(
                "e:color", (*itE)->GetIndex());
        std::vector< coordC_type > color(vc.cbegin(), vc.cend());
        line_color_coords.push_back(color);
      }
    }
  }

  // FACES IN CONTAINER
  face_range fRange = helpers::faces(const_cast< Self * >(&other));

  for(face_iterator itF = fRange.begin(); itF != fRange.end(); ++itF)
  {
    std::vector< index_type > face;
    std::vector< index_type > normal_indices;
    std::vector< index_type > texture_indices;

    vertex_face_range vfRange = helpers::incident_vertices(*itF);

    for(vertex_face_iterator itVF = vfRange.begin(); itVF != vfRange.end();
        ++itVF)
    {
      face.push_back(indexMap[*itVF]);

      // face's vertices attributes (normal & texture)
      if(useVertexNormal)
        normal_indices.push_back(indexMap[*itVF]);
      if(useVertexTextureCoord)
        texture_indices.push_back(indexMap[*itVF]);
    }

    faces_indices.push_back(face);
    if(useVertexNormal)
      normal_face_indices.push_back(normal_indices);
    if(useVertexTextureCoord)
      texture_face_indices.push_back(texture_indices);

    if(useFaceColor)
    {
      // store vertex color in standard container
      const AIFMesh::Vector &vc =
          other.GetProperty< AIFFace::ptr, AIFMesh::Vector >(
              "f:color", (*itF)->GetIndex());
      std::vector< coordC_type > color(vc.cbegin(), vc.cend());
      face_color_coords.push_back(color);
    }

    if (useFaceDatafield)
    {
      auto it_s = dff_names.begin(), it_se = dff_names.end();
      size_t i = 0;
      for (; it_s != it_se; ++it_s, ++i)
      {
        const std::vector< double >&vdata =
          other.GetProperty< AIFFace::ptr, std::vector< double > >(
            *it_s, (*itF)->GetIndex());

        if (itF == fRange.begin())
        {
          field_names.push_back(*it_s);
          field_attributes.resize(field_attributes.size() + 1);
        }
        if (useVertexDatafield && (it_s == dff_names.begin()))
          i = other.GetPropertyMapNamesStartingWithPrefix< AIFVertex::ptr >("v:datafield:").size();

        field_attributes[i].push_back(vdata);
      }
    }
  }

  ///////////////////////////////////////////////////////////////////////
  if(vertex_color_coords.size() != points_coords.size())
    vertex_color_coords.clear();
  for(auto &color : vertex_color_coords)
  {
    if(color.size() < 3 || color.size() > 4)
    {
      std::cout << "AIFMesh copy constructor: found vertex color with size != "
                   "3 and !=4. Disabling vertex color for all vertices."
                << std::endl;
      vertex_color_coords.clear();
      break;
    }
  }

  if(face_color_coords.size() != faces_indices.size())
    face_color_coords.clear();
  for(auto &color : face_color_coords)
  {
    if(color.size() < 3 || color.size() > 4)
    {
      std::cout << "AIFMesh copy constructor: found face color with size != 3 "
                   "and !=4. Disabling face color for all faces."
                << std::endl;
      face_color_coords.clear();
      break;
    }
  }

  if(!texture_file_name.empty())
  {
    // create a property map to store the texture name
    auto tfn_pm = this->AddAssocPropertyMap< helpers::ptr_mesh, std::string >(
        "m:texturefilename");
    if(!this->isAssocPropertyMap("m:texturefilename"))
      throw std::runtime_error(
          "Failed to create texture-filename property map.");

    (*tfn_pm)[this] = texture_file_name;

    std::cout << "AIF property map for texture-filename created." << std::endl;
  }

  // import data into AIF datastructure
  bool useCornerTextureCoord = false;
  bool useFaceNormal = false; // per-face vertices normals
  if(!texture_coords.empty())
  {
    if(texture_coords.size() == points_coords.size())
    {
      // texture coordinates are given by vertex
      useVertexTextureCoord = true;

      // create property map to store vertex texture-coord
      this->AddPropertyMap< AIFVertex::ptr, AIFMesh::PointUV >("v:texcoord");
      if(!this->isPropertyMap< AIFVertex::ptr >("v:texcoord"))
        throw std::runtime_error(
            "Failed to create vertex-texture-coordinate property map.");

      std::cout << "AIF property map for vertex-texture-coordinate created."
                << std::endl;
    }
    else if(texture_face_indices.size() == faces_indices.size())
    {
      // texture coordinates are given by corner
      useCornerTextureCoord = true;

      // create property map to store corner texture-coord
      this->AddAssocPropertyMap< helpers::halfedge_descriptor,
                                 AIFMesh::PointUV >("h:texcoord");
      if(!this->isAssocPropertyMap("h:texcoord"))
        throw std::runtime_error(
            "Failed to create halfedge-texture-coordinate property map.");

      std::cout << "AIF property map for halfedge-texture-coordinate created."
                << std::endl;
    }
  }

  if(!normals_coords.empty())
  {
    // per-vertex normal or per-face-vertex normal ?
    // if there are as many normals as vertices, then it's per-vertex normal

    if(normals_coords.size() == points_coords.size())
    {
      useVertexNormal = true;

      // create property map to store vertex normal
      this->AddPropertyMap< AIFVertex::ptr, AIFMesh::Normal >("v:normal");
      if(!this->isPropertyMap< AIFVertex::ptr >("v:normal"))
        throw std::runtime_error(
            "Failed to create vertex-normal property map.");

      std::cout << "AIF property map for vertex-normal created." << std::endl;
    }
    else if(!normal_face_indices[0].empty())
    {
      // per-face-vertex normal
      useFaceNormal = true;

      // create property map to store face-vertex normal
      this->AddPropertyMap< AIFFace::ptr, AIFMesh::Normal >("f:normal");
      if(!this->isPropertyMap< AIFFace::ptr >("f:normal"))
        throw std::runtime_error("Failed to create face-normal property map.");

      std::cout << "AIF property map for face-normal created." << std::endl;
    }
  }

  if(!vertex_color_coords.empty())
  {
    // create property map to store vertex color
    this->AddPropertyMap< AIFVertex::ptr, AIFMesh::Vector >("v:color");
    if(!this->isPropertyMap< AIFVertex::ptr >("v:color"))
      throw std::runtime_error("Failed to create vertex-color property map.");

    std::cout << "AIF property map for vertex-color created." << std::endl;
  }
  if(!face_color_coords.empty())
  {
    // create property map to store vertex color
    this->AddPropertyMap< AIFFace::ptr, AIFMesh::Vector >("f:color");
    if(!this->isPropertyMap< AIFFace::ptr >("f:color"))
      throw std::runtime_error("Failed to create face-color property map.");

    std::cout << "AIF property map for face-color created." << std::endl;
  }

  if (!field_names.empty())
  {
    auto it = field_names.begin();
    auto ite = field_names.end();
    for (; it != ite; ++it)
    {
      if (it->find("v:datafield:") != std::string::npos)
      { // property map associated to vertices

        // create property map to store vertex current data field
        this->AddPropertyMap< AIFVertex::ptr, std::vector<double> >(
         *it);
        if (!this->isPropertyMap< AIFVertex::ptr >(*it))
          throw std::runtime_error(
            "Failed to create a vertex data field property map.");
      }
      else if (it->find("f:datafield:") != std::string::npos)
      { // property map associated to faces

        // create property map to store face current data field
        this->AddPropertyMap< AIFFace::ptr, std::vector<double> >(
          *it);
        if (!this->isPropertyMap< AIFFace::ptr >(*it))
          throw std::runtime_error(
            "Failed to create a face data field property map.");
      }
    }
  }
  
  
  /////////////////////////////// VERTICES ADDITION
  //////////////////////////////////////
  std::vector< helpers::vertex_type::ptr > vertices;
  vertices.reserve(points_coords.size());

  typedef std::vector< std::vector< coord_type > >::const_iterator point_it;
  int index = 0;
  for(point_it itP = points_coords.begin(); itP != points_coords.end(); ++itP)
  {
    helpers::vertex_descriptor currentVertex = helpers::add_vertex(this);
    index++;

    PropHelpers::set_point(
        this, currentVertex, (*itP)[0], (*itP)[1], (*itP)[2]);

    vertices.push_back(
        currentVertex); // Warning : must keep the vertex in a vector to handle
                        // vertex indices for faces

    // Vertex normal/color/texture-coords addition when present
    if(useVertexNormal || useVertexColor || useVertexTextureCoord)
    {
      if(useVertexNormal)
      {
        // DBG std::cout << "store normal by vertex #" << vertexNbr <<
        // std::endl;
        std::vector< coordN_type > &coords = normals_coords[index - 1];

        if(coords.size() != 3)
          throw std::runtime_error("Support of normal coordinates of size != 3 "
                                   "not yet implemented!");
        // store normal in property map
        AIFMesh::Normal vn(coords[0], coords[1], coords[2]);
        this->SetProperty< AIFVertex::ptr, AIFMesh::Normal >(
            "v:normal", currentVertex->GetIndex(), vn);
      }
      if(useVertexColor)
      {
        // DBG std::cout << "store color for vertex #" << index << std::endl;
        std::vector< coordC_type > &coords = vertex_color_coords[index - 1];
        // at this point we are sure that coords.size() == 3 because of
        // previous sanity checks

        // store color in property map
        AIFMesh::Vector vc(coords[0], coords[1], coords[2]);
        this->SetProperty< AIFVertex::ptr, AIFMesh::Vector >(
            "v:color", currentVertex->GetIndex(), vc);
      }
      if(useVertexTextureCoord)
      {
        // DBG std::cout << "store tex coord for vertex #" << vertexNbr <<
        // std::endl;
        std::vector< coordT_type > &coords = texture_coords[index - 1];

        if(coords.size() != 2)
          throw std::runtime_error("Support of texture coordinates of size != "
                                   "2 not yet implemented!");

        // store texture coordinates in property map
        AIFMesh::PointUV texCoords = {coords[0], coords[1]};
        this->SetProperty< AIFVertex::ptr, AIFMesh::PointUV >(
            "v:texcoord", currentVertex->GetIndex(), texCoords);
      }
    }
    if (!field_names.empty())
    {
      auto it = field_names.begin();
      auto ite = field_names.end();
      auto it_d = field_attributes.begin();
      for (; it != ite; ++it, ++it_d)
      {
        if (it->find("v:datafield:") != std::string::npos)
        { // property map associated to vertices

          this->SetProperty< AIFVertex::ptr, std::vector<double> >(
            *it, currentVertex->GetIndex(), (*it_d)[index - 1]);
        }
      }
    }
  }


  /////////////////////////////// LINES ADDITION
  //////////////////////////////////////
  typedef std::vector< std::vector< index_type > >::const_iterator line_it;
  index = 0;
  for(line_it itL = lines_indices.begin(); itL != lines_indices.end(); ++itL)
  {
    // Dangling edge
    // create the edge in the mesh
    helpers::edge_type::ptr danglingEdge = helpers::add_edge(this);

    // find the corresponding vertices in the vertices std::vector (care about
    // the potential decrementation)
    helpers::vertex_descriptor firstVertex = vertices[(*itL)[0]];
    helpers::vertex_descriptor secondVertex = vertices[(*itL)[1]];
    // link the edge and the vertices
    helpers::link_vertex_and_edge(
        firstVertex, danglingEdge, helpers::vertex_pos::FIRST);
    helpers::link_vertex_and_edge(
        secondVertex, danglingEdge, helpers::vertex_pos::SECOND);

    if(useLineColor)
    {
      std::vector< coordC_type > &coords = line_color_coords[index++];

      // store color in property map
      AIFMesh::Vector vc(coords[0], coords[1], coords[2]);
      this->SetProperty< AIFEdge::ptr, AIFMesh::Vector >(
          "e:color", danglingEdge->GetIndex(), vc);
    }
  }


  /////////////////////////////// FACES ADDITION
  //////////////////////////////////////
  typedef std::vector< std::vector< index_type > >::const_iterator face_it;
  typedef std::vector< index_type >::const_iterator index_it;
  size_t faceId = 0;
  for(face_it itF = faces_indices.begin(); itF != faces_indices.end();
      ++itF, ++faceId)
  {
    if(itF->size() < 2)
    {
      throw std::runtime_error("AIFMesh copy constructor -> find a face with "
                               "less than two vertices : not complying case.");
    }
    else if(itF->size() == 2)
    {
      // Dangling edge
      // create the edge in the mesh
      helpers::edge_type::ptr danglingEdge = helpers::add_edge(this);
      // find the corresponding vertices in the vertices std::vector (care about
      // the potential decrementation)
      helpers::vertex_descriptor firstVertex = vertices[(*itF)[0]];
      helpers::vertex_descriptor secondVertex = vertices[(*itF)[1]];
      // link the edge and the vertices
      helpers::link_vertex_and_edge(
          firstVertex, danglingEdge, helpers::vertex_pos::FIRST);
      helpers::link_vertex_and_edge(
          secondVertex, danglingEdge, helpers::vertex_pos::SECOND);
    }
    else
    {
      // create the face in the mesh
      helpers::face_descriptor currentFace = helpers::add_face(this);

      // get the first vertex of the face to successively constructs the face's
      // edges (we can use back() here as vertices is explicitly an std::vector)
      helpers::vertex_descriptor prevVertex = vertices[itF->back()];
      helpers::vertex_descriptor currentVertex;
      helpers::edge_descriptor currentEdge;

      // loop over face's vertices
      index_type faceVertexId = 0;
      for(index_it it = itF->begin(); it != itF->end(); ++it, ++faceVertexId)
      {
        currentVertex = vertices[(*it)];

        // check if the 2 vertices ain't already link by an edge
        helpers::edge_descriptor currentEdge;
        if((currentEdge = helpers::common_edge(prevVertex, currentVertex)) ==
           helpers::null_edge())
        {
          // create the edge in the mesh
          currentEdge = helpers::add_edge(this);
          // link the edge and the vertices
          helpers::link_vertex_and_edge(
              prevVertex, currentEdge, helpers::vertex_pos::FIRST);
          helpers::link_vertex_and_edge(
              currentVertex, currentEdge, helpers::vertex_pos::SECOND);
        }
        // link the edge and the face
        helpers::link_edge_and_face(currentEdge, currentFace);

        // corner (aka halfedge) texture
        if(useCornerTextureCoord)
        {
          // DBG std::cout << "store tex coord for face #" << faceId << " and
          // vertex #" << faceVertexId << std::endl;

          // retrieve texture coords
          index_type textureId = texture_face_indices[faceId][faceVertexId];
          std::vector< coordT_type > &coords = texture_coords[textureId];
          if(coords.size() != 2)
            throw std::runtime_error("Support of texture coordinates of size "
                                     "!= 2 not yet implemented!");
          AIFMesh::PointUV texCoords = {coords[0], coords[1]};

          // retrieve halfedge
          helpers::halfedge_descriptor he =
              helpers::AIFHalfEdge(prevVertex, currentVertex, currentFace);

          // fill in property map
          auto pm = this->GetAssocPropertyMap< helpers::halfedge_descriptor,
                                               AIFMesh::PointUV >("h:texcoord");
          (*pm)[he] = texCoords;

          // DBG std::cout << __FILE__ << ":" << __LINE__ << " he=" << he << "
          // (*pm)[he][0]=" << (*pm)[he][0] << "   (*pm)[he][1]=" << (*pm)[he][1]
          // << "   coords[0]=" << coords[0] << "   coords[1]=" << coords[1] << "
          // texCoords[0]=" << texCoords[0] << "   texCoords[1]=" << texCoords[1]
          // << std::endl;
        }

        // face-vertex-normal
        if(useFaceNormal)
        {
          // DBG std::cout << "compute normal for face #" << faceId << " (vertex
          // #" << faceVertexId << ")" << std::endl;

          // retrieve current vertex normal
          index_type normalId = normal_face_indices[faceId][faceVertexId];
          std::vector< coordN_type > &normal = normals_coords[normalId];
          if(normal.size() != 3)
            throw std::runtime_error(
                "Support of normals of size != 3 not yet implemented!");
          AIFMesh::Normal vnormal = {normal[0], normal[1], normal[2]};

          // retrieve face normal
          auto pm =
              this->GetPropertyMap< AIFFace::ptr, AIFMesh::Normal >("f:normal");
          AIFMesh::Normal fnormal = (*pm)[currentFace];

          // computation of face normal (mean of the vertices normals)
          fnormal = fnormal + vnormal;
          if(it == itF->end() - 1)
            fnormal =
                fnormal / (faceVertexId +
                           1); // why to compute face normal not based on its
                               // underlying plane and just use the vertex
                               // normals to ensure correct orientation?

          // update face-normal map
          (*pm)[currentFace] = fnormal;
        }

        prevVertex = currentVertex;
      }

      // face-color
      if(useFaceColor)
      {
        std::vector< coordC_type > &coords = face_color_coords[faceId];

        // store color in property map
        AIFMesh::Vector vc(coords[0], coords[1], coords[2]);
        this->SetProperty< AIFFace::ptr, AIFMesh::Vector >(
            "f:color", currentFace->GetIndex(), vc);
      }
      if (!field_names.empty())
      {
        auto it = field_names.begin();
        auto ite = field_names.end();
        auto it_d = field_attributes.begin();
        for (; it != ite; ++it, ++it_d)
        {
          if (it->find("f:datafield:") != std::string::npos)
          { // property map associated to faces

            this->SetProperty< AIFFace::ptr, std::vector<double> >(
              *it, currentFace->GetIndex(), (*it_d)[faceId]);
          }
        }
      }
    }
  }
  ///////////////////////////////////////////////////////////////////////////
  return *this;
}

inline
void 
AIFMesh::clear() {
  while (m_Faces.begin() != m_Faces.end())
    EraseIsolatedFace(*(m_Faces.begin()));
  while (m_Edges.begin() != m_Edges.end())
    EraseIsolatedEdge(*(m_Edges.begin()));
  while (m_Vertices.begin() != m_Vertices.end())
    EraseIsolatedVertex(*(m_Vertices.begin()));

  m_VertexPropertyMaps.clear(); m_EdgePropertyMaps.clear();
  m_FacePropertyMaps.clear(); m_AssocPropertyMaps.clear();
}

inline
AIFMesh::ptr_mesh
AIFMesh::New()
{
  ptr_mesh ptr(new Self());
  return ptr;
}

inline
void
AIFMesh::Print() const
{
  std::cout << " Nb of vertices = " << GetNumberOfVertices() << std::endl;
  std::cout << " Nb of edges = " << GetNumberOfEdges() << std::endl;
  std::cout << " Nb of faces = " << GetNumberOfFaces() << std::endl;

  if(GetNumberOfFaces() <= 30)
  {
    for(VertexContainerType::const_iterator it = GetVertices().begin();
        it != GetVertices().end();
        ++it)
      (*it)->Print();

    std::cout << "-----------------------------------------------------------"
              << std::endl;

    for(EdgeContainerType::const_iterator it = GetEdges().begin();
        it != GetEdges().end();
        ++it)
      (*it)->Print();

    std::cout << "-----------------------------------------------------------"
              << std::endl;

    for(FaceContainerType::const_iterator it = GetFaces().begin();
        it != GetFaces().end();
        ++it)
      (*it)->Print();
  }
}


// specializations of GetPropertyMapContainer()

template<>
inline
PropertyMapContainer *
AIFMesh::GetPropertyMapContainer< AIFVertex::ptr >(void)
{
  return &m_VertexPropertyMaps;
}

template<>
inline
PropertyMapContainer *
AIFMesh::GetPropertyMapContainer< AIFEdge::ptr >(void)
{
  return &m_EdgePropertyMaps;
}

template<>
inline
PropertyMapContainer *
AIFMesh::GetPropertyMapContainer< AIFFace::ptr >(void)
{
  return &m_FacePropertyMaps;
}

template<>
inline
const PropertyMapContainer *
AIFMesh::GetPropertyMapContainer< AIFVertex::ptr >(void) const
{
  return &m_VertexPropertyMaps;
}

template<>
inline
const PropertyMapContainer *
AIFMesh::GetPropertyMapContainer< AIFEdge::ptr >(void) const
{
  return &m_EdgePropertyMaps;
}

template<>
inline
const PropertyMapContainer *
AIFMesh::GetPropertyMapContainer< AIFFace::ptr >(void) const
{
  return &m_FacePropertyMaps;
}

} // namespace AIF
} // namespace DataStructures
} // namespace FEVV
