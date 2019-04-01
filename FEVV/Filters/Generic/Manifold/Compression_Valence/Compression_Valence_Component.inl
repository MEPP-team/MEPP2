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

#include "Compression_Valence_Component.h"
#include "Compression_Valence_Common.h"
#include <FEVV/Tools/Math/MatrixOperations.hpp>
#include <FEVV/Wrappings/Graph_traits_extension.h>
#include <FEVV/Operators/Geometry/triangles.hpp>

#include <CGAL/boost/graph/Euler_operations.h>
#include <CGAL/boost/graph/helpers.h>         // for CGAL::clear(mesh)
#include <CGAL/boost/graph/copy_face_graph.h> // for CGAL::copy_face_graph()
#include <CGAL/circulator.h>                  // for CGAL_For_all()
#include <Eigen/Dense>
#include <map>
#include <set>
#include <bitset>
#include <sstream>
#include <chrono>
#include <string>
#include <algorithm>                          // for std::min() and std::max()

#include "FEVV/Filters/Generic/generic_writer.hpp" // for FEVV::Filters::write_mesh()

#ifdef _MSC_VER
// disable some warnings on Windows
#pragma warning(push)
#pragma warning(disable : 4244)
#pragma warning(disable : 4267)
// 4244 & 4267: converting type A to type B, possible data loss
#endif


#define COLOR_NUMBER 10000
#define USE_COLOR_METRIC
#define AC_BUFFER 1024 * 10000

const int MINIMUM_PREDICTION_NUMBER = 3;
const int LIMIT_NUMBER = 50;


//#define DBG_Main_Function

template< typename HalfedgeGraph, typename PointMap, typename VertexColorMap >
std::string
Compression_Valence_Component< HalfedgeGraph, PointMap, VertexColorMap >::
    Main_Function(HalfedgeGraph &_pMesh,
                  PointMap *_pm,
                  VertexColorMap *_v_cm,
                  const std::string &_Input_File_Name,
                  const std::string &_File_Name, // output
                  const int &_Qbit,
                  const int &_NVertices,
                  const bool _Normal_flipping,
                  const bool _Use_metric,
                  const float &_Metric_thread,
                  const bool _Use_forget_metric,
                  const int &_Forget_value,
                  const bool _Compression_selected,
                  const bool _Adaptive_quantization,
                  const bool _Is_bijection_selected)
{
#ifdef DBG_Main_Function
  // std::cout << std::hexfloat;
  DBG_print_mesh_geometry(_pMesh, _pm, "beginning of Main_Function()");
  DBG_print_mesh_vertexcolor(_pMesh, _v_cm, "beginning of Main_Function()");
#endif

  // ensure the mesh has triangular faces only
  if(!CGAL::is_triangle_mesh(_pMesh))
    return "Some faces of the mesh are not triangles. No compression was done.";

  // initialize property maps
  vertex_color_int =
      FEVV::make_vertex_property_map< HalfedgeGraph, Color_Unit >(_pMesh);
  vertex_Seed_Edge =
      FEVV::make_vertex_property_map< HalfedgeGraph, int >(_pMesh);
  vertex_Component_Number =
      FEVV::make_vertex_property_map< HalfedgeGraph, int >(_pMesh);
  Vertex_Flag = FEVV::make_vertex_property_map< HalfedgeGraph, int >(_pMesh);
  Vertex_Number = FEVV::make_vertex_property_map< HalfedgeGraph, int >(_pMesh);
  Vertex_Sign = FEVV::make_vertex_property_map< HalfedgeGraph, int >(_pMesh);
  vertex_normal =
      FEVV::make_vertex_property_map< HalfedgeGraph, Vector >(_pMesh);
  vertex_Q_Index = FEVV::make_vertex_property_map< HalfedgeGraph, int >(_pMesh);
  vertex_Region_Number =
      FEVV::make_vertex_property_map< HalfedgeGraph, int >(_pMesh);
  vertex_Removal_Order =
      FEVV::make_vertex_property_map< HalfedgeGraph, int >(_pMesh);

  facet_tag = FEVV::make_face_property_map< HalfedgeGraph, int >(_pMesh);
  facet_Component_Number =
      FEVV::make_face_property_map< HalfedgeGraph, int >(_pMesh);
  Facet_Flag = FEVV::make_face_property_map< HalfedgeGraph, int >(_pMesh);
  facet_normal = FEVV::make_face_property_map< HalfedgeGraph, Vector >(_pMesh);


  // ELO  in Mepp1, compute_normals() is called when the mesh is loaded ;
  // ELO  added here to ensure normals are computed before the processing
  // ELO  begins ; maybe this line should go in Global_Initialization() ?
  compute_normals(_pMesh, _pm); // ELO+

  // start time measurement
  auto time_start = std::chrono::steady_clock::now();

  // read size of input mesh (file size)
  if(FILE *file = fopen(_Input_File_Name.c_str(), "r"))
  {
    fseek(file, 0, SEEK_END);
    this->Initial_file_size = ftell(file);
    fclose(file);
  }
  else
    this->Initial_file_size = 0;

  // Use of bijection or not for the geometry encoding
  this->Is_Bijection_Enabled = _Is_bijection_selected;

  unsigned Init_number_vertices = (unsigned)FEVV::size_of_vertices(_pMesh);

  // Initialization - Quantization, Color, Multiple components
  this->Global_Initialization(_pMesh, _Qbit, _File_Name.c_str(), _pm, _v_cm);

#ifdef DBG_Main_Function
  DBG_print_mesh_vertexcolor(
      _pMesh, _v_cm, "in Main_Function() after Global_Initialization()");
#endif

  // When use of adaptive quantization is selected
  if(_Adaptive_quantization)
  {
    this->Adaptive_Quantization(_pMesh,
                                _pm,
                                _v_cm,
                                _NVertices,
                                _Normal_flipping,
                                _Use_metric,
                                _Metric_thread,
                                _Use_forget_metric,
                                _Forget_value,
                                _Qbit);
  }
  else
  {
    this->Simplification(_pMesh,
                         _pm,
                         _NVertices,
                         _Normal_flipping,
                         _Use_metric,
                         _Metric_thread,
                         _Use_forget_metric,
                         _Forget_value);
  }

  unsigned Connectivity_size = 0, Color_size = 0, Total_size = 0;

  // Compression
  if(_Compression_selected)
  {
    this->Compression(_pMesh,
                      _File_Name.c_str(),
                      _Qbit,
                      Connectivity_size,
                      Color_size,
                      Total_size,
                      _pm); //, this->Initial_file_size);
  }

  unsigned Number_layers = this->GlobalCountOperation;
  unsigned Final_number_vertices = (unsigned)FEVV::size_of_vertices(_pMesh);

  // stop time measurement
  auto time_end = std::chrono::steady_clock::now();
  std::chrono::duration< double > time_diff = time_end - time_start;

  // To show result
  double Connectivity_rate = (double)Connectivity_size / Init_number_vertices;
  double Color_rate = (double)Color_size / Init_number_vertices;
  double Total_rate = (double)Total_size * 8 / Init_number_vertices;
  double Geometry_rate = Total_rate - Connectivity_rate - Color_rate;

  std::ostringstream Res_tmp;
  Res_tmp << "Base mesh : " << std::setw(3) << Final_number_vertices
          << " vertices \n";
  Res_tmp << "Connectivity : " << float(Connectivity_rate) << " b/v \n";
  Res_tmp << "Geometry : " << float(Geometry_rate) << " b/v\n";
  Res_tmp << "Color : " << float(Color_rate) << " b/v\n";
  Res_tmp << "Total size : " << float(Total_rate) << " b/v\n";
  if(this->Initial_file_size != 0)
    Res_tmp << "Ratio : " << (float)Total_size / this->Initial_file_size * 100
            << " % \n\n";
  else
    Res_tmp << "Ratio : "
            << "unavailable\n";
  Res_tmp << "Number of layers : " << Number_layers << "\n";
  Res_tmp << "Calculation time : " << (float)time_diff.count() << " seconds \n";
  std::string Res = Res_tmp.str();

  // ELO+beg
  std::cout << Res << std::endl;
#ifdef DBG_Main_Function
  DBG_print_mesh_geometry(_pMesh, _pm, "end of Main_Function()");
  DBG_print_mesh_vertexcolor(_pMesh, _v_cm, "end of Main_Function()");
#endif
  // ELO+end

  return Res;
}


// Description : To select the input gate.
template< typename HalfedgeGraph, typename PointMap, typename VertexColorMap >
void
Compression_Valence_Component< HalfedgeGraph, PointMap, VertexColorMap >::
    Global_Initialization(HalfedgeGraph &_pMesh,
                          const int &_Qbit,
                          const char *_File_Name,
                          const PointMap *pm,
                          VertexColorMap *v_cm)
{
  // (1) Determination if the mesh is colored.
  // (2) Conversion of color space.
  // (3) Quantization of converted colors into "vertex->color_int()"
  // (4) Re-paint mesh with re_calculated rgb colors from quantized converted
  // color. (5) Establish color palette. (6) Color quantization - Descrease of
  // possible color numbers
  this->Color_Initialization(_pMesh, v_cm);

  // Initialization of multiple components (quantization is performed separately
  // for each component)
  this->Multiple_Components_Initialization(_pMesh, pm, _Qbit);

  // Quantization of each component
  this->Quantization(_pMesh, pm);
}


//#define DBG_Multiple_Components_Initialization

template< typename HalfedgeGraph, typename PointMap, typename VertexColorMap >
void
Compression_Valence_Component< HalfedgeGraph, PointMap, VertexColorMap >::
    Multiple_Components_Initialization(HalfedgeGraph &_pMesh,
                                       const PointMap *pm,
                                       int const &_Qbit)
{
  typedef boost::graph_traits< HalfedgeGraph > GraphTraits;
  typedef typename GraphTraits::vertex_iterator vertex_iterator;
  typedef typename GraphTraits::face_iterator face_iterator;
  typedef typename FEVV::Geometry_traits< HalfedgeGraph >::Point Point3d;
  typedef typename FEVV::Geometry_traits< HalfedgeGraph >::Vector Vector;
  FEVV::Geometry_traits< HalfedgeGraph > gt(_pMesh);

  // Initialize vertex flags
  auto vertex_iterator_pair = vertices(_pMesh);
  vertex_iterator pVertex = vertex_iterator_pair.first;
  vertex_iterator pVertex_end = vertex_iterator_pair.second;
  for(; pVertex != pVertex_end; ++pVertex)
  {
    put(this->vertex_Seed_Edge, *pVertex, OTHER_COORDINATE);
    put(this->vertex_Component_Number, *pVertex, -1);
  }

  // Several components
  // (1) To get the number of components;
  // (2) To get the volume, area and the number of vertices of each component;
  // (3) To know if each component is closed or not.
  // (4) Tag vertices and facets to correspoding component number

  // Initialize face tag
  auto face_iterator_pair = faces(_pMesh);
  face_iterator pFacet_beg = face_iterator_pair.first;
  face_iterator pFacet_end = face_iterator_pair.second;
  for(face_iterator pFacet = pFacet_beg; pFacet != pFacet_end; ++pFacet)
  {
    put(this->facet_tag, *pFacet, -1);
  }

  int Component_index = 0;

  for(face_iterator pFacet = pFacet_beg; pFacet != pFacet_end; ++pFacet)
  {

    if(get(this->facet_tag, *pFacet) == -1)
    {
      bool Is_closed = true;
      bool Is_seed_edge_found = false;

      float xmin = 50000., ymin = 50000., zmin = 50000.;
      float xmax = -50000., ymax = -50000., zmax = -50000.;
      double area = 0;
      int Number_vertices = 0;

      put(this->facet_tag, *pFacet, Component_index);

      std::list< face_descriptor > facets;
      facets.push_front(*pFacet);

      while(!facets.empty())
      {
        face_descriptor F = facets.front();
        facets.pop_front();

        put(this->facet_tag, F, Component_index);

        // tag component number to facet
        put(this->facet_Component_Number, F, Component_index);

        area += Area_Facet_Triangle(halfedge(F, _pMesh), _pMesh, pm);


        CGAL::Halfedge_around_face_circulator< HalfedgeGraph > pHalfedge(
            halfedge(F, _pMesh), _pMesh);
        CGAL::Halfedge_around_face_circulator< HalfedgeGraph > end(pHalfedge);
        do
        {
          // tag the vertex to its corresponding component number
          if(get(this->vertex_Component_Number, target(*pHalfedge, _pMesh)) ==
             -1)
          {
            put(this->vertex_Component_Number,
                target(*pHalfedge, _pMesh),
                Component_index);
            Number_vertices++;
          }

          // To calculate the bounding box of each component
          Point3d p = get(*pm, target(*pHalfedge, _pMesh));
          if(p[0] > xmax)
            xmax = p[0];

          if(p[1] > ymax)
            ymax = p[1];

          if(p[2] > zmax)
            zmax = p[2];

          if(p[0] < xmin)
            xmin = p[0];

          if(p[1] < ymin)
            ymin = p[1];

          if(p[2] < zmin)
            zmin = p[2];

          // To know if the component is closed or not
          if(CGAL::is_border_edge(*pHalfedge, _pMesh))
            Is_closed = false;

          if(!Is_seed_edge_found)
          {
            if(((!CGAL::is_border_edge(*pHalfedge, _pMesh)) &&
                (!Is_Border_Vertex(*pHalfedge, _pMesh)) &&
                (!Is_Border_Vertex(opposite(*pHalfedge, _pMesh), _pMesh))) ||
               (out_degree(target(next(*pHalfedge, _pMesh), _pMesh), _pMesh) !=
                6))
            {
              Is_seed_edge_found = true;

              // Seed edge of each component
              put(this->vertex_Seed_Edge,
                  target(*pHalfedge, _pMesh),
                  2 * Component_index);
              put(this->vertex_Seed_Edge,
                  target(opposite(*pHalfedge, _pMesh), _pMesh),
                  2 * Component_index + 1);
            }
          }

          face_descriptor pNFacet = face(opposite(*pHalfedge, _pMesh), _pMesh);
          if(pNFacet != boost::graph_traits< HalfedgeGraph >::null_face() &&
             get(this->facet_tag, pNFacet) == -1)
          {
            facets.push_front(pNFacet);
            put(this->facet_tag, pNFacet, Component_index);
            put(this->facet_Component_Number, pNFacet, Component_index);
          }
        } while(++pHalfedge != end);
      }

      this->xmin.push_back(xmin);
      this->ymin.push_back(ymin);
      this->zmin.push_back(zmin);

      this->xmax.push_back(xmax);
      this->ymax.push_back(ymax);
      this->zmax.push_back(zmax);

      double HighestBB = -5000;
      if(xmax - xmin > HighestBB)
        HighestBB = xmax - xmin;
      if(ymax - ymin > HighestBB)
        HighestBB = ymax - ymin;
      if(zmax - zmin > HighestBB)
        HighestBB = zmax - zmin;

      this->HighestLengthBB.push_back(HighestBB);

      // volume, area, number of vertices of each component
      Vector e1(xmax - xmin, 0., 0.);
      Vector e2(0., ymax - ymin, 0.);

      Vector normal = FEVV::Math::Vector::cross_product<
          FEVV::Geometry_traits< HalfedgeGraph > >(e1, e2, gt);
      double base_area = FEVV::Math::Vector::l2_distance<
          FEVV::Geometry_traits< HalfedgeGraph > >(normal, gt);
      double volume = base_area * (zmax - zmin);

      // For quasi-optimal determination of quantization precision.
      // Originally, the highest dimension of bounding box was considered as 10.
      // We shoud take this point into account.
      volume *= std::pow(((double)10.0 / (double)HighestBB), 3.0);
      area *= std::pow(((double)10.0 / (double)HighestBB), 2.0);

      // Stock information for each component
      this->ComponentVolume.push_back(volume);
      this->ComponentArea.push_back(area);
      this->ComponentNumberVertices.push_back(Number_vertices);

      // To know if each component is open or closed.
      this->IsClosed.push_back(Is_closed);

      Component_index++;
    }
  }

  this->NumberComponents = Component_index;

  std::list< int > li;
  std::list< Point_Int > pi;
  std::list< Color_Unit > cu;

  // Initilization of containers for each component
  for(int i = 0; i < this->NumberComponents; i++)
  {
    // operation (decimation or quantization)
    this->ListOperation.push_back(li);

    // connectivity
    this->Connectivity.push_back(li);

    // geometry
    this->Geometry.push_back(pi);

    // vertex color
    this->VertexColor.push_back(cu);

    // Number of connectivity and geometry symbols
    this->NumberSymbol.push_back(li);
    this->NumberVertices.push_back(li);

    // Range of alpha and gamma
    this->AlphaRange.push_back(li);
    this->AlphaOffset.push_back(li);
    this->GammaRange.push_back(li);
    this->GammaOffset.push_back(li);

    // Number of operations for each component
    this->ComponentOperations.push_back(0);

    // Number of decimations for each component
    this->NumberDecimation.push_back(0);

    // Number of under quantization
    this->NumberChangeQuantization.push_back(0);
    this->NumberColorQuantization.push_back(0);

    // Qbit for each component
    this->Qbit.push_back(_Qbit);

    // Displacement vector for under quantization for each component
    this->QuantizationCorrectVector.push_back(li);

    // The number of vertices for each under_quantization
    this->NumberQuantizationLayer.push_back(li);

    // For color quantization change
    this->NumberProcessedVertices.push_back(li);
    this->ColorChildcellIndex.push_back(li);
    this->ColorEncoderIndex.push_back(li);
  }

#ifdef DBG_Multiple_Components_Initialization
  {
    std::cout << "DBG "
              << "--- function " << __func__ << std::endl;
    std::cout << "DBG "
              << "this->NumberComponents=" << this->NumberComponents
              << std::endl;

    {
      std::cout << "DBG "
                << "vertex->Seed_Edge=" << std::endl;

      auto vertex_iterator_pair = vertices(_pMesh);
      vertex_iterator pVertex = vertex_iterator_pair.first;
      vertex_iterator pVertex_end = vertex_iterator_pair.second;
      for(; pVertex != pVertex_end; ++pVertex)
      {
        auto s = get(this->vertex_Seed_Edge, *pVertex);
        std::cout << "DBG   " << s << std::endl;
      }
    }

    {
      std::cout << "DBG "
                << "vertex->Component_Number=" << std::endl;

      auto vertex_iterator_pair = vertices(_pMesh);
      vertex_iterator pVertex = vertex_iterator_pair.first;
      vertex_iterator pVertex_end = vertex_iterator_pair.second;
      for(; pVertex != pVertex_end; ++pVertex)
      {
        auto n = get(this->vertex_Component_Number, *pVertex);
        std::cout << "DBG   " << n << std::endl;
      }
    }

    {
      std::cout << "DBG "
                << "face->tag=" << std::endl;

      auto face_iterator_pair = faces(_pMesh);
      face_iterator pFacet = face_iterator_pair.first;
      face_iterator pFacet_end = face_iterator_pair.second;
      for(; pFacet != pFacet_end; ++pFacet)
      {
        auto t = get(this->facet_tag, *pFacet);
        std::cout << "DBG   " << t << std::endl;
      }
    }

    {
      std::cout << "DBG "
                << "face->Component_Number=" << std::endl;

      auto face_iterator_pair = faces(_pMesh);
      face_iterator pFacet = face_iterator_pair.first;
      face_iterator pFacet_end = face_iterator_pair.second;
      for(; pFacet != pFacet_end; ++pFacet)
      {
        auto n = get(this->facet_Component_Number, *pFacet);
        std::cout << "DBG   " << n << std::endl;
      }
    }
  }
#endif
}

//#define DBG_Quantization

/*
        Description : Quantize all vertices so that the new positions
        are reguliraly spaced in the 3D space. */
template< typename HalfedgeGraph, typename PointMap, typename VertexColorMap >
void
Compression_Valence_Component< HalfedgeGraph, PointMap, VertexColorMap >::
    Quantization(HalfedgeGraph &_pMesh, const PointMap *pm)
{
  typedef boost::graph_traits< HalfedgeGraph > GraphTraits;
  typedef typename GraphTraits::vertex_iterator vertex_iterator;
  typedef typename FEVV::Geometry_traits< HalfedgeGraph >::Point Point3d;

  // Quantization step for each component
  for(int i = 0; i < this->NumberComponents; i++)
  {
    float max = this->xmax[i] - this->xmin[i];

    if(this->ymax[i] - this->ymin[i] > max)
      max = this->ymax[i] - this->ymin[i];

    if(this->zmax[i] - this->zmin[i] > max)
      max = this->zmax[i] - this->zmin[i];

    int NbInteraval = std::pow(2., (int)this->Qbit[i]);

    float Q_Step = max / (float)NbInteraval;
    this->Quantization_Step.push_back(Q_Step);
  }

  // Vertex quantization
  auto vertex_iterator_pair = vertices(_pMesh);
  vertex_iterator pVert_beg = vertex_iterator_pair.first;
  vertex_iterator pVert_end = vertex_iterator_pair.second;
  for(vertex_iterator pVert = pVert_beg; pVert != pVert_end; ++pVert)
  {
    Point3d point = get(*pm, *pVert);
#if 0
		//TODO-elo-restore ?
		double x = point[0];
		double y = point[1];
		double z = point[2];
#else
    // ELO+
    // ELO  simulate reading points coordinates in floats
    // ELO  then converting to double as it is done by Mepp1
    float xtmp = point[0];
    float ytmp = point[1];
    float ztmp = point[2];
    double x = xtmp;
    double y = ytmp;
    double z = ztmp;
#endif

    int Component_ID = get(this->vertex_Component_Number, *pVert);

    int Qx = (int)(ceil((x - (double)this->xmin[Component_ID]) /
                        (double)this->Quantization_Step[Component_ID])) -
             1;
    if(Qx == -1)
      Qx = 0;
    int Qy = (int)(ceil((y - (double)this->ymin[Component_ID]) /
                        (double)this->Quantization_Step[Component_ID])) -
             1;
    if(Qy == -1)
      Qy = 0;
    int Qz = (int)(ceil((z - (double)this->zmin[Component_ID]) /
                        (double)this->Quantization_Step[Component_ID])) -
             1;
    if(Qz == -1)
      Qz = 0;

    put(*pm,
        *pVert,
        Point3d(this->xmin[Component_ID] +
                    (Qx + 0.5) * this->Quantization_Step[Component_ID],
                this->ymin[Component_ID] +
                    (Qy + 0.5) * this->Quantization_Step[Component_ID],
                this->zmin[Component_ID] +
                    (Qz + 0.5) * this->Quantization_Step[Component_ID]));
  }

#ifdef DBG_Quantization
  {
    std::cout << "DBG "
              << "--- function " << __func__ << std::endl;

    {
      std::cout << "DBG " << __func__ << " "
                << "vertex->point()=" << std::endl;

      auto vertex_iterator_pair = vertices(_pMesh);
      vertex_iterator pVertex = vertex_iterator_pair.first;
      vertex_iterator pVertex_end = vertex_iterator_pair.second;
      for(; pVertex != pVertex_end; ++pVertex)
      {
        auto p3d = get(*pm, *pVertex);
        std::cout << "DBG " << __func__ << " " << p3d[0] << " " << p3d[1] << " "
                  << p3d[2] << std::endl;
      }
    }
  }
#endif
}

//#define DBG_Color_Initialization

// this->ColorArray -> contains all initial colors present in the input mesh.
template< typename HalfedgeGraph, typename PointMap, typename VertexColorMap >
void
Compression_Valence_Component< HalfedgeGraph, PointMap, VertexColorMap >::
    Color_Initialization(HalfedgeGraph &_pMesh, VertexColorMap *v_cm)
{
  typedef boost::graph_traits< HalfedgeGraph > GraphTraits;
  typedef typename GraphTraits::vertex_iterator vertex_iterator;
  typedef typename boost::property_traits< VertexColorMap >::value_type Color;

  // (1) To determine if the mesh is colored.
  //     We consider the mesh is colored if a colormap is provided.
  if(v_cm != nullptr)
  {
    this->IsColored = true;

    auto vertex_iterator_pair = vertices(_pMesh);
    vertex_iterator vi = vertex_iterator_pair.first;
    vertex_iterator vi_end = vertex_iterator_pair.second;

    // keep the first color in case there is a unique color
    Color c = get(*v_cm, *vi);
    this->OnlyColor[0] = c[0];
    this->OnlyColor[1] = c[1];
    this->OnlyColor[2] = c[2];

    // check if color is unique for the whole mesh
    this->IsOneColor = true;
    for(; vi != vi_end; ++vi) // loop over vertices
    {
      Color c = get(*v_cm, *vi);
      if(c[0] != this->OnlyColor[0] || c[1] != this->OnlyColor[1] ||
         c[2] != this->OnlyColor[2])
      {
        this->IsOneColor = false;
        break;
      }
    }
  }

  // (2) If the mesh is colored -> Color initialization.
  //     Conversion to Lab colors and determine max and min for each components.
  if((this->IsColored) && (!this->IsOneColor))
  {
    float Temp_color[3];

    float C0_min = 5000, C1_min = 5000, C2_min = 5000;
    float C0_max = -5000, C1_max = -5000, C2_max = -5000;

    // calculate new color and find min values.
    auto vertex_iterator_pair = vertices(_pMesh);
    vertex_iterator pVertex = vertex_iterator_pair.first;
    vertex_iterator pVertex_end = vertex_iterator_pair.second;
    for(; pVertex != pVertex_end; ++pVertex)
    {
      Color c = get(*v_cm, *pVertex);
      Temp_color[0] = c[0];
      Temp_color[1] = c[1];
      Temp_color[2] = c[2];

      // Color space conversion.
      float New_color[3];
      RGB_To_LAB(Temp_color[0], Temp_color[1], Temp_color[2], New_color);

      // assignment of color in new color space in order to quantize later.
      put(*v_cm, *pVertex, Color(New_color[0], New_color[1], New_color[2]));

      // calculate min and max to identify quantization step.
      if(New_color[0] > C0_max)
        C0_max = New_color[0];
      if(New_color[0] < C0_min)
        C0_min = New_color[0];

      if(New_color[1] > C1_max)
        C1_max = New_color[1];
      if(New_color[1] < C1_min)
        C1_min = New_color[1];

      if(New_color[2] > C2_max)
        C2_max = New_color[2];
      if(New_color[2] < C2_min)
        C2_min = New_color[2];
    }

    const int Nb_interaval = (int)std::pow(2.0, C0_QUANTIZATION) - 1;

    float Color_max = C0_max - C0_min;
    if(C1_max - C1_min > Color_max)
      Color_max = C1_max - C1_min;
    if(C2_max - C2_min > Color_max)
      Color_max = C2_max - C2_min;

    // Information needed to quantize converted color.
    this->C0_Min = C0_min;
    this->C1_Min = C1_min;
    this->C2_Min = C2_min;

    // Step size of color quantization
    this->Color_Quantization_Step = (float)((Color_max) / Nb_interaval);

    // Enter quantized color valued into "vertex->color_int" to use lated.
    // Also, recalculate vertex color using these quantized color values.
    Color_Unit Resulting_color;
    float New_vertex_color[3];
    float Reconstructed_color[3];

    // int Color_index = 0;
    pVertex = vertex_iterator_pair.first;
    for(; pVertex != pVertex_end; ++pVertex)
    {
      Color c = get(*v_cm, *pVertex);
      Temp_color[0] = c[0];
      Temp_color[1] = c[1];
      Temp_color[2] = c[2];

      int Qc0 = (int)(floor(
          (Temp_color[0] - C0_min) / this->Color_Quantization_Step + 0.5));
      int Qc1 = (int)(floor(
          (Temp_color[1] - C1_min) / this->Color_Quantization_Step + 0.5));
      int Qc2 = (int)(floor(
          (Temp_color[2] - C2_min) / this->Color_Quantization_Step + 0.5));

      Resulting_color.c0 = Qc0;
      Resulting_color.c1 = Qc1;
      Resulting_color.c2 = Qc2;
      put(this->vertex_color_int,
          *pVertex,
          Color_Unit(
              Resulting_color.c0, Resulting_color.c1, Resulting_color.c2));

      New_vertex_color[0] =
          this->C0_Min + Resulting_color.c0 * this->Color_Quantization_Step;
      New_vertex_color[1] =
          this->C1_Min + Resulting_color.c1 * this->Color_Quantization_Step;
      New_vertex_color[2] =
          this->C2_Min + Resulting_color.c2 * this->Color_Quantization_Step;

      LAB_To_RGB(New_vertex_color[0],
                 New_vertex_color[1],
                 New_vertex_color[2],
                 Reconstructed_color);
      for(int i = 0; i < 3; i++)
      {
        if(Reconstructed_color[i] < 0.)
          Reconstructed_color[i] = 0.;
        if(Reconstructed_color[i] > 1.)
          Reconstructed_color[i] = 1.;
      }

      // re-paint the input mesh with reconstructed colors from Lab to RGB
      // transformation.
      put(*v_cm,
          *pVertex,
          Color(Reconstructed_color[0],
                Reconstructed_color[1],
                Reconstructed_color[2]));
    }
  }

#ifdef DBG_Color_Initialization
  {
    std::cout << "DBG "
              << "--- function " << __func__ << std::endl;
    std::cout << "DBG "
              << "this->IsColored=" << this->IsColored << std::endl;
    std::cout << "DBG "
              << "this->IsOneColor=" << this->IsOneColor << std::endl;

    if(this->IsColored)
    {
      std::cout << "DBG "
                << "vertex->color()=" << std::endl;

      auto vertex_iterator_pair = vertices(_pMesh);
      vertex_iterator pVertex = vertex_iterator_pair.first;
      vertex_iterator pVertex_end = vertex_iterator_pair.second;
      for(; pVertex != pVertex_end; ++pVertex)
      {
        Color c = get(*v_cm, *pVertex);
        std::cout << "DBG   " << c[0] << " " << c[1] << " " << c[2]
                  << std::endl;
      }
    }

    if((this->IsColored) && (!this->IsOneColor))
    {
      std::cout << "DBG "
                << "vertex->color_int()=" << std::endl;

      auto vertex_iterator_pair = vertices(_pMesh);
      vertex_iterator pVertex = vertex_iterator_pair.first;
      vertex_iterator pVertex_end = vertex_iterator_pair.second;
      for(; pVertex != pVertex_end; ++pVertex)
      {
        Color_Unit c = get(this->vertex_color_int, *pVertex);
        std::cout << "DBG   " << c.c0 << " " << c.c1 << " " << c.c2
                  << std::endl;
      }
    }
  }
#endif
}


template< typename HalfedgeGraph, typename PointMap, typename VertexColorMap >
void
Compression_Valence_Component< HalfedgeGraph, PointMap, VertexColorMap >::
    Adaptive_Quantization(HalfedgeGraph &_pMesh,
                          PointMap *_pm,
                          VertexColorMap *_v_cm,
                          const int &_NVertices,
                          const bool _Normal_flipping,
                          const bool _Use_metric,
                          const float &_Metric_thread,
                          const bool _Use_forget_metric,
                          const int &_Forget_value,
                          const int &_Qbit)
{


  if(!this->IsColored)
  {
    // bool Is_any_vertex_removed = true;
    unsigned Last_Number = 0;
    unsigned Current_Number = FEVV::size_of_vertices(_pMesh);
    int Operation_choice = -1;

    // int Old_Q;
    // double mrms, mrmswrtBB, hausdorff, hausdorffwrtBB;
    bool Continue;
    do
    {
      Last_Number = Current_Number;
      Continue = false;

      for(int Component_ID = 0; Component_ID < this->NumberComponents;
          Component_ID++)
      {
        // if the ith component did not remove any vertex in last loop, it is
        // not necessary to simplifiy it.
        if(this->ComponentOperations[Component_ID] ==
           this->GlobalCountOperation)
        {
          int temp = 0;

          // Calculate area of each component.
          this->Recalculate_Component_Area(_pMesh, _pm, Component_ID, temp);

          // Estimated quantization precision of geometry
          int QG = Estimate_Geometry_Quantization(
              this->ComponentVolume[Component_ID],
              this->ComponentArea[Component_ID],
              this->ComponentNumberVertices[Component_ID]);

          if(QG < 4)
            QG = 4;

          // if the current precision is > QG then we apply decrease of
          // quantization precision
          if(QG < (int)this->Qbit[Component_ID]) // MT
          {
            Continue = true;
            Operation_choice = 1;

            // Reducing of quantization precision of 1 bit.
            this->Diminush_Geometry_Quantization_Precision(
                _pMesh, Component_ID, _pm);

            this->Qbit[Component_ID]--;
            this->NumberChangeQuantization[Component_ID]++;

            this->ComponentOperations[Component_ID]++;
            this->ListOperation[Component_ID].push_front(Operation_choice);
          }
          // else -> Decimation
          else
          {
            Operation_choice = 0;

            unsigned Initial_number_vertices = FEVV::size_of_vertices(_pMesh);

            // Decimation and regulation conquests.
            this->Decimation_Conquest(_pMesh,
                                      _pm,
                                      _Normal_flipping,
                                      _Use_metric,
                                      _Metric_thread,
                                      _Use_forget_metric,
                                      _Forget_value,
                                      Component_ID);
            this->Regulation(_pMesh,
                             _Normal_flipping,
                             _Use_metric,
                             _Metric_thread,
                             _Use_forget_metric,
                             _Forget_value,
                             Component_ID,
                             _pm);

            this->NumberDecimation[Component_ID] += 1;
            this->ComponentOperations[Component_ID] += 1;
            this->ListOperation[Component_ID].push_front(Operation_choice);

            unsigned Diff_number_vertices =
                FEVV::size_of_vertices(_pMesh) - Initial_number_vertices;
            this->ComponentNumberVertices[Component_ID] += Diff_number_vertices;

            if(Diff_number_vertices == 0)
              this->Remove_Last_Phase_Elements(Component_ID);

            else
              Continue = true;
          }
        }
      }

      Current_Number = FEVV::size_of_vertices(_pMesh);

      if(Continue)
        this->GlobalCountOperation++;


      if(Current_Number < (unsigned)_NVertices) // MT
        break;

    } while((Current_Number != Last_Number) || (Continue));

    compute_normals(_pMesh, _pm);
  }

  // if mesh is colored
  else
  {
    FILE *Operation_order = fopen("Operation_order.txt", "w");
    fclose(Operation_order);

    // bool Is_any_vertex_removed = true;
    unsigned Last_Number = 0;
    unsigned Current_Number = FEVV::size_of_vertices(_pMesh);
    int Operation_choice = -1;

    std::vector< int > QC_Initials;
    std::vector< int > QC_Finals;

    // To get QC_Initial and QC_Final
    for(int Component_ID = 0; Component_ID < this->NumberComponents;
        Component_ID++)
    {
      double Max_color, Mean_color;
      int Number_vertices = 0;

      // Calculate max and mean color error.
      this->Calculate_Edge_Color_Difference(
          _pMesh, Component_ID, Max_color, Mean_color, Number_vertices);
      double Mean_Max = Mean_color / Max_color;

      // int QC_Initial = floor(-56.334*Mean_Max*Mean_Max +4.6838*Mean_Max
      // +7.8632 + 0.5); int QC_Final = floor(-0.893 * log(Mean_Max*Area)
      // + 8.0957 + 0.5);

      int QC_init = floor(-57.966 * Mean_Max * Mean_Max + 5.311 * Mean_Max +
                          7.8062 + 0.5);
      int QC_fin = floor(84.548 * Mean_Max * Mean_Max - 33.723 * Mean_Max +
                         7.8222 + 0.5);

      QC_Initials.push_back(QC_init);
      QC_Finals.push_back(QC_fin);
    }


    bool Continue;
    do
    {
      Last_Number = Current_Number;
      Continue = false;

      for(int Component_ID = 0; Component_ID < this->NumberComponents;
          Component_ID++)
      {
        // if the i-th component did not remove any vertex in last loop, it is
        // not necessary to simplifiy it.
        if(this->ComponentOperations[Component_ID] ==
           this->GlobalCountOperation)
        {
          int temp = 0;
          this->Recalculate_Component_Area(_pMesh, _pm, Component_ID, temp);
          int Number_of_vertices = 0;
          double Max_color = 0., Mean_color = 0.;
          this->Calculate_Edge_Color_Difference(
              _pMesh, Component_ID, Max_color, Mean_color, Number_of_vertices);
          double Mean_Max = Mean_color / Max_color;

          // Estimated color quantization(QC) and geometry quantization (QG)
          int QC =
              floor(-1.181 * log(Mean_Max * this->ComponentArea[Component_ID] /
                                 (double)Number_of_vertices) +
                    0.3281 + 0.5);
          int QG = Estimate_Geometry_Quantization(
              this->ComponentVolume[Component_ID],
              this->ComponentArea[Component_ID],
              this->ComponentNumberVertices[Component_ID]);

          // If the current color quantization precision > QC_INIT -> Decrease
          // one bit from color quantization resolution.
          if((8 - this->NumberColorQuantization[Component_ID] >
              QC_Initials[Component_ID]) &&
             (!this->IsOneColor)) // && (QC_Final < 8 -
                                  // this->NumberColorQuantization))
          {
            Operation_choice = 2;
            Continue = true;

            // Reducing color quantization precision of 1 bit.
            this->Diminush_Color_Quantization_Precision(
                _pMesh, Component_ID, _v_cm);

            this->NumberColorQuantization[Component_ID]++;
            this->ComponentOperations[Component_ID]++;
            this->ListOperation[Component_ID].push_front(Operation_choice);
          }
          // If the current color quantization precision > Estimated precision
          // -> Decrease
          else if(((8 - this->NumberColorQuantization[Component_ID] > QC) &&
                   (QC >= QC_Finals[Component_ID])) &&
                  (!this->IsOneColor))
          {
            Operation_choice = 2;
            Continue = true;

            // Reducing color quantization precision of 1 bit.
            this->Diminush_Color_Quantization_Precision(
                _pMesh, Component_ID, _v_cm);

            this->NumberColorQuantization[Component_ID]++;
            this->ComponentOperations[Component_ID]++;
            this->ListOperation[Component_ID].push_front(Operation_choice);
          }

          // Geometry quantization precision decrease
          else if(((int)this->Qbit[Component_ID] > QG) && (QG >= 5))
          {
            Continue = true;
            Operation_choice = 1;

            // Reducuing geometry quantization precision of 1 bit.
            this->Diminush_Geometry_Quantization_Precision(
                _pMesh, Component_ID, _pm);

            this->Qbit[Component_ID]--;
            this->NumberChangeQuantization[Component_ID]++;

            this->ComponentOperations[Component_ID]++;
            this->ListOperation[Component_ID].push_front(Operation_choice);
          }

          // Else we perform decimation
          else
          {
            Operation_choice = 0;

            unsigned Initial_number_vertices = FEVV::size_of_vertices(_pMesh);

            // Decimation and regulation conquests.
            this->Decimation_Conquest(_pMesh,
                                      _pm,
                                      _Normal_flipping,
                                      _Use_metric,
                                      _Metric_thread,
                                      _Use_forget_metric,
                                      _Forget_value,
                                      Component_ID);
            this->Regulation(_pMesh,
                             _Normal_flipping,
                             _Use_metric,
                             _Metric_thread,
                             _Use_forget_metric,
                             _Forget_value,
                             Component_ID,
                             _pm);

            this->NumberDecimation[Component_ID] += 1;
            unsigned Diff_number_vertices =
                FEVV::size_of_vertices(_pMesh) - Initial_number_vertices;

            this->ComponentOperations[Component_ID] += 1;
            this->ListOperation[Component_ID].push_front(Operation_choice);

            this->ComponentNumberVertices[Component_ID] += Diff_number_vertices;

            if(Diff_number_vertices == 0)
              this->Remove_Last_Phase_Elements(Component_ID);

            else
              Continue = true;
          }
          FILE *Operation_order = fopen("Operation_order.txt", "a");
          fprintf(Operation_order,
                  "Component_ID = %d    Operation = %d \n",
                  Component_ID,
                  Operation_choice);
          fclose(Operation_order);
        }
      }

      Current_Number = FEVV::size_of_vertices(_pMesh);

      if(Continue)
        this->GlobalCountOperation++;


      if(Current_Number < (unsigned)_NVertices) // MT
        break;

    } while((Current_Number != Last_Number) || (Continue));

    compute_normals(_pMesh, _pm);
  }
}


template< typename HalfedgeGraph, typename PointMap, typename VertexColorMap >
void
Compression_Valence_Component< HalfedgeGraph, PointMap, VertexColorMap >::Init(
    HalfedgeGraph &mesh)
{
  typedef boost::graph_traits< HalfedgeGraph > GraphTraits;
  typedef typename GraphTraits::vertex_iterator vertex_iterator;
  typedef typename GraphTraits::face_iterator face_iterator;

  // vertices flags initialization
  int i = 0;
  // Initialize vertex flags
  auto vertex_iterator_pair = vertices(mesh);
  vertex_iterator pVertex_beg = vertex_iterator_pair.first;
  vertex_iterator pVertex_end = vertex_iterator_pair.second;
  for(vertex_iterator pVertex = pVertex_beg; pVertex != pVertex_end;
      ++pVertex, ++i)
  {
    put(this->Vertex_Flag, *pVertex, FREE);
    put(this->Vertex_Number, *pVertex, i);
    put(this->Vertex_Sign, *pVertex, NOSIGN);
  }

  // facets flag initialization.
  auto face_iterator_pair = faces(mesh);
  face_iterator pFacet_beg = face_iterator_pair.first;
  face_iterator pFacet_end = face_iterator_pair.second;
  for(face_iterator pFacet = pFacet_beg; pFacet != pFacet_end; ++pFacet)
  {
    put(this->Facet_Flag, *pFacet, FREE);
  }
}


#if 1 // DBG-ELO+beg
template< typename HalfedgeGraph, typename PointMap, typename VertexColorMap >
void
Compression_Valence_Component< HalfedgeGraph, PointMap, VertexColorMap >::
    print_halfedge(const std::string &title,
                   const halfedge_descriptor &h,
                   const HalfedgeGraph &_pMesh,
                   const PointMap *_pm)
{
  vertex_descriptor vtarget = target(h, _pMesh);
  vertex_descriptor vsource = source(h, _pMesh);

  Point3d psource = get(*_pm, vsource);
  Point3d ptarget = get(*_pm, vtarget);

  std::cout << title << " :" << std::endl;
  std::cout << "  psource = " << psource[0] << ", " << psource[1] << ", "
            << psource[2] << std::endl;
  std::cout << "  ptarget = " << ptarget[0] << ", " << ptarget[1] << ", "
            << ptarget[2] << std::endl;
  bool is_border_edge =
      (face(h, _pMesh) == boost::graph_traits< HalfedgeGraph >::null_face());
  std::cout << "  is_border = " << std::boolalpha << is_border_edge
            << std::endl;
}

template< typename HalfedgeGraph, typename PointMap, typename VertexColorMap >
void
Compression_Valence_Component< HalfedgeGraph, PointMap, VertexColorMap >::
    print_vertex(const std::string &title,
                 const vertex_descriptor &v,
                 const PointMap *_pm)
{
  Point3d p = get(*_pm, v);
  std::cout << title << " = " << p[0] << ", " << p[1] << ", " << p[2]
            << std::endl;
}

template< typename HalfedgeGraph, typename PointMap, typename VertexColorMap >
std::string
Compression_Valence_Component< HalfedgeGraph, PointMap, VertexColorMap >::
    vertex_to_string(const vertex_descriptor &v, const PointMap *_pm)
{
  Point3d p = get(*_pm, v);
  std::ostringstream s;
  s << "(" << p[0] << ", " << p[1] << ", " << p[2] << ")";
  return s.str();
}
#endif // DBG-ELO+end


#if 1 // DBG-ELO+beg
template< typename HalfedgeGraph, typename PointMap, typename VertexColorMap >
std::string
Compression_Valence_Component< HalfedgeGraph, PointMap, VertexColorMap >::
    halfedge_to_string(const halfedge_descriptor &h,
                       const HalfedgeGraph &_pMesh,
                       const PointMap *_pm)
{
  vertex_descriptor vtarget = target(h, _pMesh);
  vertex_descriptor vsource = source(h, _pMesh);

  std::ostringstream s;
  s << vertex_to_string(vsource, _pm) << "->" << vertex_to_string(vtarget, _pm);
  return s.str();
}
#endif // DBG-ELO+end


#if 1 // DBG-ELO+beg
template< typename HalfedgeGraph, typename PointMap, typename VertexColorMap >
bool
Compression_Valence_Component< HalfedgeGraph, PointMap, VertexColorMap >::
    v_inf_to_v(const vertex_descriptor &v1,
               const vertex_descriptor &v2,
               const PointMap *_pm)
{
  Point3d p1 = get(*_pm, v1);
  Point3d p2 = get(*_pm, v2);

  if(p1[0] < p2[0])
    return true;
  else if((p1[0] == p2[0]) && (p1[1] < p2[1]))
    return true;
  else if((p1[0] == p2[0]) && (p1[1] == p2[1]) && (p1[2] < p2[2]))
    return true;
  else
    return false;
}

template< typename HalfedgeGraph, typename PointMap, typename VertexColorMap >
std::string
Compression_Valence_Component< HalfedgeGraph, PointMap, VertexColorMap >::
    edge_to_string(const halfedge_descriptor &h,
                   const HalfedgeGraph &_pMesh,
                   const PointMap *_pm)
{
  vertex_descriptor vtarget = target(h, _pMesh);
  vertex_descriptor vsource = source(h, _pMesh);

  std::ostringstream s;
  if(v_inf_to_v(vsource, vtarget, _pm))
    s << vertex_to_string(vsource, _pm) << "->"
      << vertex_to_string(vtarget, _pm);
  else
    s << vertex_to_string(vtarget, _pm) << "->"
      << vertex_to_string(vsource, _pm);

  return s.str();
}
#endif // DBG-ELO+end


//#define DBG_Decimation_Conquest

// Description : This function select a set of independent vertices to be
// removed
template< typename HalfedgeGraph, typename PointMap, typename VertexColorMap >
int
Compression_Valence_Component< HalfedgeGraph, PointMap, VertexColorMap >::
    Decimation_Conquest(HalfedgeGraph &_pMesh,
                        const PointMap *_pm,
                        const bool Normal_flipping,
                        const bool Use_metric,
                        const float &Metric_thread,
                        const bool Use_forget_metric,
                        const int &Forget_value,
                        const int &Component_ID)
{
  typedef boost::graph_traits< HalfedgeGraph > GraphTraits;
  typedef typename GraphTraits::halfedge_iterator halfedge_iterator;

#ifdef DBG_Decimation_Conquest
  static unsigned int Decimation_Conquest_call_cnt = 0;
  Decimation_Conquest_call_cnt++;
  std::cout << "Decimation_Conquest_call_cnt = " << Decimation_Conquest_call_cnt
            << std::endl;
#endif

  // Calculate mean color and meah area for color metric
  double Max_color, Mean_color;
  int Temp_NV = 0;
  int Number_facets;
  this->Calculate_Edge_Color_Difference(
      _pMesh, Component_ID, Max_color, Mean_color, Temp_NV);
  this->Recalculate_Component_Area(_pMesh, _pm, Component_ID, Number_facets);
  double Mean_area =
      (double)this->ComponentArea[Component_ID] / (double)Number_facets;


  // Initialize vertex and face flags.
  this->Init(_pMesh);

  // to count number of independent vertices and number of symbol of
  // connectivity.
  int Number_vertices = 0;
  int Number_symbol = 0;

  // To find first edge.
  auto halfedge_iterator_pair = halfedges(_pMesh);
  halfedge_iterator hi = halfedge_iterator_pair.first;

  while(
      (get(this->vertex_Seed_Edge, target(*hi, _pMesh)) != 2 * Component_ID) ||
      (get(this->vertex_Seed_Edge, target(opposite(*hi, _pMesh), _pMesh)) !=
       2 * Component_ID + 1))
    hi++;

  halfedge_descriptor First_halfedge = *hi;
#ifdef DBG_Decimation_Conquest
  std::cout << "marker #1" << std::endl;
  print_halfedge("First_halfedge", First_halfedge, _pMesh, _pm);
  int loop_counter = 0;
#endif

  // Two vertices of seed edge are flaged CONQUERED
  put(this->Vertex_Flag, target(First_halfedge, _pMesh), CONQUERED);
  put(this->Vertex_Flag,
      target(opposite(First_halfedge, _pMesh), _pMesh),
      CONQUERED);

  // These vertices are also flaged with sign flags for retriangulation
  put(this->Vertex_Sign, target(First_halfedge, _pMesh), PLUS);
  put(this->Vertex_Sign,
      target(opposite(First_halfedge, _pMesh), _pMesh),
      MINUS);

  std::queue< halfedge_descriptor > Halfedges; // halfedge queue.
  Halfedges.push(First_halfedge); // push the first halfedge in the queue.

  halfedge_descriptor h; // The current gate

  /// Main loop
  while(!Halfedges.empty())
  {
#ifdef DBG_Decimation_Conquest
    std::cout << "marker #2" << std::endl;
    std::cout << "loop " << ++loop_counter << std::endl;
#endif
    h = Halfedges.front();
#ifdef DBG_Decimation_Conquest
    print_halfedge("h #2", h, _pMesh, _pm);
#endif
    Halfedges.pop();

    unsigned type = 0; // define type of retriangulation

    int valence = (int)out_degree(target(next(h, _pMesh), _pMesh),
                                  _pMesh); // valence

    if(CGAL::is_border(h, _pMesh) == true)
      continue;
#ifdef DBG_Decimation_Conquest
    std::cout << "marker #3" << std::endl;
    std::cout << "(h->next()->vertex()->Vertex_Flag == FREE)"
              << (get(this->Vertex_Flag, target(next(h, _pMesh), _pMesh)) ==
                  FREE)
              << std::endl;
    std::cout << "(valence >= 3)" << (valence >= 3) << std::endl;
    std::cout << "(valence <= 4)" << (valence <= 4) << std::endl;
    std::cout << "(Is_Border_Vertex(h->next()) == true)"
              << (Is_Border_Vertex(next(h, _pMesh), _pMesh) == true)
              << std::endl;
#endif

    // if its front face is not tagged CONQUERED nor TO_BE_REMOVED, do nothing!!
    if((get(this->Facet_Flag, face(h, _pMesh)) == CONQUERED) ||
       (get(this->Facet_Flag, face(h, _pMesh)) == TO_BE_REMOVED))
    {
#ifdef DBG_Decimation_Conquest
      std::cout << "marker #4" << std::endl;
#endif
      continue;
    }
    // if its front vertex is free and has a valence <= 6 and it is not a border
    // vertex.
    else if((get(this->Vertex_Flag, target(next(h, _pMesh), _pMesh)) == FREE) &&
            (valence >= 3) && (valence <= 6) &&
            (!Is_Border_Vertex(next(h, _pMesh), _pMesh)))
    {
#ifdef DBG_Decimation_Conquest
      std::cout << "marker #5" << std::endl;
#endif

      type = Find_Type(_pMesh, h, valence);

      // Check if the manifold property is violated.
      bool Manifold_condition =
          Is_Manifold_Property_Violated(_pMesh, h, type, valence);

      bool Normal_flipping_condition = false;

      if(Normal_flipping)
        Normal_flipping_condition =
            Is_Normal_Flipping_Occured(_pMesh, _pm, h, valence);

      // calculate error caused by the removal. This metric decides if the
      // vertex can be removed or not.
      bool Geometric_metric_condition = false;

      if(Use_metric)
      {
        if(Use_forget_metric)
        {
          if(FEVV::size_of_vertices(_pMesh) > (unsigned)Forget_value)
            Geometric_metric_condition = false;
          else
            Geometric_metric_condition = Is_Geometric_Metric_Violated(
                _pMesh, _pm, h, type, valence, Metric_thread);
        }
        else
          Geometric_metric_condition = Is_Geometric_Metric_Violated(
              _pMesh, _pm, h, type, valence, Metric_thread);
      }
      bool Is_Color_Too_Important = false;

#ifdef USE_COLOR_METRIC
      if(this->IsColored)
        Is_Color_Too_Important = this->Error_Projected_Surface(
            _pMesh, _pm, h, Component_ID, Mean_color, Mean_area);
#endif
#ifdef DBG_Decimation_Conquest
      std::cout << "Is_Color_Too_Important = " << Is_Color_Too_Important
                << std::endl;
#endif

      // remove the front vertex if its removal does not viloate the manifold
      // property and some metrics
      bool Check_all_condition = false;

      if((!Manifold_condition) && (!Geometric_metric_condition) &&
         (!Normal_flipping_condition) && (!Is_Color_Too_Important))
        Check_all_condition = true;

      if(Check_all_condition) // All conditions are good. -> Remove the center
                              // vertex.
      {
#ifdef DBG_Decimation_Conquest
        std::cout << "marker #6" << std::endl;
#endif

        // increase number of vertices and symbols
        Number_vertices++;
        Number_symbol++;

        halfedge_descriptor g = h;

        Point3d Coordinates_removed_vertex =
            get(*_pm, target(next(g, _pMesh), _pMesh));
        Point_Int CRV =
            Change_Real_Int(Coordinates_removed_vertex, Component_ID);

        // encoding of vertex color
        if((this->IsColored) && (!this->IsOneColor))
        {
          Color_Unit Removed_vertex_color;
          Removed_vertex_color = Get_Vertex_Color(next(g, _pMesh), _pMesh);

          Color_Unit Average_color;

          Average_color = Get_Average_Vertex_Color_Lee(g, valence, _pMesh);


          // Color difference from average color of neighbors
          Color_Unit Color_diff = Removed_vertex_color - Average_color;
          this->InterVertexColor.push_front(Color_diff);
        }


        // Enter symbol 'VALENCE CODE' into the list of symbols
        this->InterConnectivity.push_front(valence - 3);

        // Calculate the position of barycenter.
        Point3d Barycenter = Barycenter_Patch_Before_Removal(g, _pMesh, _pm);

        Point_Int BC = Change_Real_Int(Barycenter, Component_ID);

        // remove the front vertex
#ifdef DBG_Decimation_Conquest
        print_halfedge("remove vertex", next(g, _pMesh), _pMesh, _pm);
#endif
        CGAL::Euler::remove_center_vertex(next(g, _pMesh), _pMesh);

        g = h;
        put(this->Facet_Flag, face(g, _pMesh), TEMP_FLAG);

        // Flags and fill queue
        for(int j = 0; j < (valence - 1); j++)
        {
          g = next(g, _pMesh);
          put(this->Vertex_Flag, target(g, _pMesh), CONQUERED);

          put(this->Vertex_Flag,
              target(opposite(g, _pMesh), _pMesh),
              CONQUERED);

          if(!(CGAL::is_border_edge(g, _pMesh)))
            Halfedges.push(opposite(g, _pMesh));
        }
        g = h;

        Retriangulation(_pMesh, g, valence, 0, Component_ID, _pm);

        Vector normal = Normal_Patch(g, valence, _pMesh, _pm);
        Vector T2 = Vector(0.0, 0.0, 0.0);
        Vector T1 = Calculate_T1_T2(h, normal, T2, _pMesh, _pm);

        if(normal == Vector(0.0, 0.0, 0.0))
        {
          T1 = Vector(1., 0., 0.);
          T2 = Vector(0., 1., 0.);
          normal = Vector(0., 0., 1.);
        }
        else if(T1 == Vector(0.0, 0.0, 0.0))
        {
          T1 = Vector(1., 0., 0.);
          T2 = Vector(0., 1., 0.);
          normal = Vector(0., 0., 1.);
        }
        else if(T2 == Vector(0.0, 0.0, 0.0))
        {
          T1 = Vector(1., 0., 0.);
          T2 = Vector(0., 1., 0.);
          normal = Vector(0., 0., 1.);
        }

        Point_Int Dist = CRV - BC;
        Point_Int Frenet_Coordinates;
        // Bijection
        if(this->Is_Bijection_Enabled)
          Frenet_Coordinates = Frenet_Rotation(Dist, T1, T2, normal);
        else
          Frenet_Coordinates = Dist;

        this->InterGeometry.push_front(Frenet_Coordinates);
      }

      // Conditions are not satisfied -> NULL CODE
      else
      {
#ifdef DBG_Decimation_Conquest
        std::cout << "marker #7" << std::endl;
#endif
        // Enter symbol 'NULL PATCH' into the list of symbols
        this->InterConnectivity.push_front(4);

        Number_symbol++;
        put(this->Facet_Flag, face(h, _pMesh), CONQUERED);
        put(this->Vertex_Flag, target(next(h, _pMesh), _pMesh), CONQUERED);

        if((get(this->Vertex_Sign, target(h, _pMesh)) == PLUS) &&
           (get(this->Vertex_Sign, target(opposite(h, _pMesh), _pMesh)) ==
            MINUS))
        {
          if(get(this->Vertex_Sign, target(next(h, _pMesh), _pMesh)) == NOSIGN)
            put(this->Vertex_Sign, target(next(h, _pMesh), _pMesh), PLUS);
        }
        else if((get(this->Vertex_Sign, target(h, _pMesh)) == MINUS) &&
                (get(this->Vertex_Sign, target(opposite(h, _pMesh), _pMesh)) ==
                 PLUS))
        {
          if(get(this->Vertex_Sign, target(next(h, _pMesh), _pMesh)) == NOSIGN)
            put(this->Vertex_Sign, target(next(h, _pMesh), _pMesh), PLUS);
        }
        else if((get(this->Vertex_Sign, target(h, _pMesh)) == PLUS) &&
                (get(this->Vertex_Sign, target(opposite(h, _pMesh), _pMesh)) ==
                 PLUS))
        {
          if(get(this->Vertex_Sign, target(next(h, _pMesh), _pMesh)) == NOSIGN)
            put(this->Vertex_Sign, target(next(h, _pMesh), _pMesh), MINUS);
        }
        else if((get(this->Vertex_Sign, target(h, _pMesh)) == MINUS) &&
                (get(this->Vertex_Sign, target(opposite(h, _pMesh), _pMesh)) ==
                 MINUS))
        {
          if(get(this->Vertex_Sign, target(next(h, _pMesh), _pMesh)) == NOSIGN)
            put(this->Vertex_Sign, target(next(h, _pMesh), _pMesh), PLUS);
        }

        if(!CGAL::is_border_edge(next(h, _pMesh), _pMesh))
          Halfedges.push(opposite(next(h, _pMesh), _pMesh));

        if(!CGAL::is_border_edge(prev(h, _pMesh), _pMesh))
          Halfedges.push(opposite(prev(h, _pMesh), _pMesh));
      }
    }


    // border edge.
    else if((get(this->Vertex_Flag, target(next(h, _pMesh), _pMesh)) == FREE) &&
            (valence >= 3) && (valence <= 4) &&
            (Is_Border_Vertex(next(h, _pMesh), _pMesh) == true))
    {
#ifdef DBG_Decimation_Conquest
      std::cout << "marker #8" << std::endl;
#endif
      /*****    conditions of vertex removal (based on area) will be added
       * *****/
      /*****    conditions of vertex removal (based on area) will be added
       * *****/

      type = Find_Type(_pMesh, h, valence);
#ifdef DBG_Decimation_Conquest
      print_halfedge("h #8", h, _pMesh, _pm);
#endif

      halfedge_descriptor g = h;

      bool Check_border_structure = true;

      halfedge_descriptor First_border_edge;
      halfedge_descriptor Standard_edge;


      // To find first border edge.
      CGAL::Halfedge_around_target_circulator< HalfedgeGraph > Hvc(
          target(next(g, _pMesh), _pMesh), _pMesh);
      CGAL::Halfedge_around_target_circulator< HalfedgeGraph > Hvc_end = Hvc;

      int Number_neighboring_border_vertices = 0;
      CGAL_For_all(Hvc, Hvc_end)
      {
        // ELO  in this loop the halfedges are not traversed in the same order
        // ELO  as in Mepp1 ; the direction of rotation is the same

        if(Is_Border_Vertex(opposite(*Hvc, _pMesh), _pMesh))
        {
          Number_neighboring_border_vertices++;
        }

        if(CGAL::is_border(*Hvc, _pMesh))
        {
          First_border_edge = *Hvc;
        }
      }
#ifdef DBG_Decimation_Conquest
      print_halfedge("First_border_edge", First_border_edge, _pMesh, _pm);
      std::cout << "Number_neighboring_border_vertices = "
                << Number_neighboring_border_vertices << std::endl;
#endif

      if(Number_neighboring_border_vertices > 2)
        Check_border_structure = false;

      if(Is_Manifold_Property_Violated(_pMesh, g, type, valence))
        Check_border_structure = false;
      // Added by Adrien Maglo.
      // Test if we are not in the case of a mesh pimple.
      if(Check_border_structure)
      {
#ifdef DBG_Decimation_Conquest
        std::cout << "marker #9" << std::endl;
        // Get the two patch border vertices.
        print_halfedge("First_border_edge", First_border_edge, _pMesh, _pm);
#endif
        vertex_descriptor vh1 = target(next(First_border_edge, _pMesh), _pMesh);
        vertex_descriptor vh2 =
            target(opposite(First_border_edge, _pMesh), _pMesh);
        CGAL::Halfedge_around_target_circulator< HalfedgeGraph > vh_it =
            CGAL::Halfedge_around_target_circulator< HalfedgeGraph >(vh1,
                                                                     _pMesh);
        CGAL::Halfedge_around_target_circulator< HalfedgeGraph > vh_it_end =
            vh_it;

        // Test if the two patch border vertices are not connected
        // by an edge.
        CGAL_For_all(vh_it, vh_it_end)
        {
          if(target(opposite(*vh_it, _pMesh), _pMesh) == vh2)
          {
            Check_border_structure = false;
            break;
          }
        }
#ifdef DBG_Decimation_Conquest
        std::cout << "Check_border_structure = " << std::boolalpha
                  << Check_border_structure << std::endl;
#endif
      }


      // if the hole is a triangle(Count_border_edges == 3), we should not
      // remove the vertex.
      int Count_border_edges = 1;
      Standard_edge = First_border_edge;
      First_border_edge = next(First_border_edge, _pMesh);

      for(; First_border_edge != Standard_edge;
          First_border_edge = next(First_border_edge, _pMesh))
      {
        Count_border_edges++;
        if(Count_border_edges > 5)
          break;
      }

      if(Count_border_edges <= 3)
        Check_border_structure = false;

      // Exception which causes a manifold violation
      if(Count_border_edges == 4)
        if(out_degree(
               target(opposite(prev(Standard_edge, _pMesh), _pMesh), _pMesh),
               _pMesh) == 2)
          Check_border_structure = false;


      /*if (Check_border_structure)
      {
              if(valence == 3)
              {
                      if(Is_Border_Vertex(g))
                      {
                              if(g->opposite()->vertex_degree() == 3)
                                      Check_border_structure = false;
                      }
                      else
                      {
                              if(g->vertex_degree() == 3)
                                      Check_border_structure = false;
                      }
              }
      }*/


      if(Check_border_structure)
      {
#ifdef DBG_Decimation_Conquest
        std::cout << "marker #10" << std::endl;
#endif
        Number_vertices++;
        Number_symbol++;

        // Connectivity Code = 5 for border valence 3 vertex
        //         "         = 6 for border valence 4 vertex
        this->InterConnectivity.push_front(valence + 2);

        Point3d Real_vertex_position =
            get(*_pm, target(next(h, _pMesh), _pMesh));
        Point_Int Vertex_position =
            Change_Real_Int(Real_vertex_position, Component_ID);

        Color_Unit Removed_vertex_color;
        // int Vertex_color_index = -1;
        if((this->IsColored) && (!this->IsOneColor))
        {
          Removed_vertex_color = Get_Vertex_Color(next(h, _pMesh), _pMesh);
        }


        std::vector< halfedge_descriptor > Border_edges;
        int Number_jump = 0;

        while(CGAL::is_border_edge(next(g, _pMesh), _pMesh) == false)
        {
          g = next(opposite(next(g, _pMesh), _pMesh), _pMesh);
          Number_jump++;
        }

        Border_edges.push_back(opposite(g, _pMesh));
#ifdef DBG_Decimation_Conquest
        print_halfedge(
            "Border_edges.push_back 1", opposite(g, _pMesh), _pMesh, _pm);
#endif
        for(int i = 0; i < (valence - 2); i++)
        {
          g = prev(opposite(prev(g, _pMesh), _pMesh), _pMesh);
          Border_edges.push_back(opposite(g, _pMesh));
#ifdef DBG_Decimation_Conquest
          print_halfedge(
              "Border_edges.push_back 2", opposite(g, _pMesh), _pMesh, _pm);
#endif
        }

        for(int i = 0; i < (valence - 1); i++)
        {
          halfedge_descriptor Temp = Border_edges[i];
          Temp = opposite(Temp, _pMesh);

          CGAL_assertion(!CGAL::is_border(Temp, _pMesh));
#ifdef DBG_Decimation_Conquest
          std::cout << "remove face" << std::endl;
#endif
          CGAL::Euler::remove_face(Temp, _pMesh);
        }

        if(valence == 3)
        {
#ifdef DBG_Decimation_Conquest
          std::cout << "marker #11" << std::endl;
#endif
          halfedge_descriptor Retriangulation_edge =
              prev(opposite(Border_edges[valence - 2], _pMesh), _pMesh);


          // One triangle has to be created to smooth the mesh boundary
#ifdef DBG_Decimation_Conquest
          {
            bool is_border_Retriangulation_edge =
                (face(Retriangulation_edge, _pMesh) ==
                 boost::graph_traits< HalfedgeGraph >::null_face());

            print_halfedge("Border_edges[valence - 2]",
                           Border_edges[valence - 2],
                           _pMesh,
                           _pm);
            print_halfedge(
                "Retriangulation_edge", Retriangulation_edge, _pMesh, _pm);
          }
          std::cout << "add face" << std::endl;
#endif
          CGAL::Euler::add_face_to_border(
              Retriangulation_edge, opposite(Border_edges[0], _pMesh), _pMesh);

          halfedge_descriptor Input_gate =
              opposite(Border_edges[Number_jump], _pMesh);

          // its front face is tagged CONQUERED
          put(this->Facet_Flag, face(Input_gate, _pMesh), CONQUERED);
          put(this->Vertex_Flag,
              target(next(Input_gate, _pMesh), _pMesh),
              CONQUERED);

          if((type == 1) || (type == 2) || (type == 4))
          {
            if(get(this->Vertex_Sign,
                   target(next(Input_gate, _pMesh), _pMesh)) == NOSIGN)
              put(this->Vertex_Sign,
                  target(next(Input_gate, _pMesh), _pMesh),
                  PLUS);
          }
          else if(type == 3)
          {
            if(get(this->Vertex_Sign,
                   target(next(Input_gate, _pMesh), _pMesh)) == NOSIGN)
              put(this->Vertex_Sign,
                  target(next(Input_gate, _pMesh), _pMesh),
                  MINUS);
          }


          Point3d Barycenter =
              Barycenter_Patch_After_Removal(Input_gate, valence, _pMesh, _pm);
          Point_Int BC = Change_Real_Int(Barycenter, Component_ID);

          if((this->IsColored) && (!this->IsOneColor))
          {
            Color_Unit Average_color;

            Average_color = Get_Average_Vertex_Color_After_Removal(
                g, valence, _pMesh); // g is the most left placed edge
            Color_Unit Color_diff = Removed_vertex_color - Average_color;

            this->InterVertexColor.push_front(Color_diff);
          }

          halfedge_descriptor g = Input_gate;
          Vector normal = Normal_Patch(g, valence, _pMesh, _pm);
          Vector T2 = Vector(0.0, 0.0, 0.0);

          Vector T1 = Calculate_T1_T2(g, normal, T2, _pMesh, _pm);

          if(normal == Vector(0.0, 0.0, 0.0))
          {
            normal = Vector(0, 0, 1);
            T1 = Vector(1, 0, 0);
            T2 = Vector(0, 1, 0);
          }

          Point_Int Dist = Vertex_position - BC;

          // bijection

          Point_Int Frenet_Coordinates;
          if(this->Is_Bijection_Enabled)
            Frenet_Coordinates = Frenet_Rotation(Dist, T1, T2, normal);
          else
            Frenet_Coordinates = Dist;


          // Point_Int Frenet_Coordinates = Frenet_Rotation(Dist,T1,T2,normal);

          this->InterGeometry.push_front(Frenet_Coordinates);

          if(Number_jump == 0)
            Halfedges.push(Border_edges[1]);
          else
            Halfedges.push(Border_edges[0]);
        }

        else if(valence == 4)
        {
#ifdef DBG_Decimation_Conquest
          std::cout << "marker #12" << std::endl;
#endif
          halfedge_descriptor Retriangulation_edge =
              prev(opposite(Border_edges[valence - 2], _pMesh), _pMesh);

          if(((Number_jump == 0) && ((type == 5) || (type == 8))) ||
             ((Number_jump == 1) && ((type == 6) || (type == 7))) ||
             ((Number_jump == 2) && ((type == 5) || (type == 8))))

          {
#ifdef DBG_Decimation_Conquest
            std::cout << "add face" << std::endl;
#endif
            CGAL::Euler::add_face_to_border(Retriangulation_edge,
                                            opposite(Border_edges[1], _pMesh),
                                            _pMesh);
            put(this->Facet_Flag,
                face(opposite(Border_edges[1], _pMesh), _pMesh),
                CONQUERED);

#ifdef DBG_Decimation_Conquest
            std::cout << "add face" << std::endl;
#endif
            CGAL::Euler::add_face_to_border(Retriangulation_edge,
                                            opposite(Border_edges[0], _pMesh),
                                            _pMesh);
            put(this->Facet_Flag,
                face(opposite(Border_edges[0], _pMesh), _pMesh),
                CONQUERED);
          }

          else if(((Number_jump == 0) && ((type == 6) || (type == 7))) ||
                  ((Number_jump == 1) && ((type == 5) || (type == 8))) ||
                  ((Number_jump == 2) && ((type == 6) || (type == 7))))
          {
#ifdef DBG_Decimation_Conquest
            std::cout << "add face" << std::endl;
#endif
            CGAL::Euler::add_face_to_border(opposite(Border_edges[2], _pMesh),
                                            opposite(Border_edges[0], _pMesh),
                                            _pMesh);
            put(this->Facet_Flag,
                face(opposite(Border_edges[1], _pMesh), _pMesh),
                CONQUERED);
            halfedge_descriptor Temp_border =
                next(opposite(Border_edges[2], _pMesh), _pMesh);

#ifdef DBG_Decimation_Conquest
            std::cout << "add face" << std::endl;
#endif
            CGAL::Euler::add_face_to_border(
                Retriangulation_edge, Temp_border, _pMesh);
            put(this->Facet_Flag,
                face(opposite(Border_edges[2], _pMesh), _pMesh),
                CONQUERED);
          }

          halfedge_descriptor Input_gate =
              opposite(Border_edges[Number_jump], _pMesh);


          // Vertex Signs
          if((type == 5) || (type == 8))
          {
            halfedge_descriptor g = Input_gate;
            g = next(opposite(prev(g, _pMesh), _pMesh), _pMesh);
            if(get(this->Vertex_Sign, target(g, _pMesh)) == NOSIGN)
              put(this->Vertex_Sign, target(g, _pMesh), MINUS);
            if(get(this->Vertex_Sign, target(opposite(g, _pMesh), _pMesh)) ==
               NOSIGN)
              put(this->Vertex_Sign, target(opposite(g, _pMesh), _pMesh), PLUS);
          }
          else if((type == 6) || (type == 7))
          {
            halfedge_descriptor g = Input_gate;
            g = prev(opposite(next(g, _pMesh), _pMesh), _pMesh);
            if(get(this->Vertex_Sign, target(g, _pMesh)) == NOSIGN)
              put(this->Vertex_Sign, target(g, _pMesh), PLUS);
            if(get(this->Vertex_Sign, target(opposite(g, _pMesh), _pMesh)) ==
               NOSIGN)
              put(this->Vertex_Sign,
                  target(opposite(g, _pMesh), _pMesh),
                  MINUS);
          }

          // Vertex Flags + Fill Halfedges queue.
          for(int i = 0; i < valence - 1; i++)
          {
            put(this->Vertex_Flag, target(Border_edges[i], _pMesh), CONQUERED);
            put(this->Vertex_Flag,
                target(opposite(Border_edges[i], _pMesh), _pMesh),
                CONQUERED);

            if((unsigned)i != (unsigned)Number_jump)
              Halfedges.push(Border_edges[i]);
          }

          if((this->IsColored) && (!this->IsOneColor))
          {
            Color_Unit Average_color;

            Color_Unit Color_0 =
                Get_Vertex_Color(opposite(Border_edges[0], _pMesh), _pMesh);
            Color_Unit Color_1 = Get_Vertex_Color(Border_edges[0], _pMesh);
            Color_Unit Color_2 = Get_Vertex_Color(Border_edges[1], _pMesh);
            Color_Unit Color_3 = Get_Vertex_Color(Border_edges[2], _pMesh);
            Average_color.c0 =
                (int)(Color_0.c0 + Color_1.c0 + Color_2.c0 + Color_3.c0) / 4.0;
            Average_color.c1 =
                (int)(Color_0.c1 + Color_1.c1 + Color_2.c1 + Color_3.c1) / 4.0;
            Average_color.c2 =
                (int)(Color_0.c2 + Color_1.c2 + Color_2.c2 + Color_3.c2) / 4.0;

            Color_Unit Color_diff = Removed_vertex_color - Average_color;

            this->InterVertexColor.push_front(Color_diff);
          }


          // encoding of geometry
          Point3d P0, P1, P2, P3;

          P0 = get(*_pm, target(opposite(Border_edges[0], _pMesh), _pMesh));
          P1 = get(*_pm, target(Border_edges[0], _pMesh));
          P2 = get(*_pm, target(Border_edges[1], _pMesh));
          P3 = get(*_pm, target(Border_edges[2], _pMesh));

          Point3d Barycenter((P0[0] + P1[0] + P2[0] + P3[0]) / 4,
                             (P0[1] + P1[1] + P2[1] + P3[1]) / 4,
                             (P0[2] + P1[2] + P2[2] + P3[2]) / 4);

          Point_Int BC = Change_Real_Int(Barycenter, Component_ID);

          halfedge_descriptor g = Input_gate;
          Vector normal = Normal_Patch(g, valence, _pMesh, _pm);

          Vector T2 = Vector(0.0, 0.0, 0.0);
          Vector T1 = Calculate_T1_T2(g, normal, T2, _pMesh, _pm);

          if(normal == Vector(0.0, 0.0, 0.0))
          {
            T1 = Vector(1, 0, 0);
            T2 = Vector(0, 1, 0);
            normal = Vector(0, 0, 1);
          }

          Point_Int Dist = Vertex_position - BC;
          Point_Int Frenet_Coordinates;
          if(this->Is_Bijection_Enabled)
            Frenet_Coordinates = Frenet_Rotation(Dist, T1, T2, normal);
          else
            Frenet_Coordinates = Dist;

          this->InterGeometry.push_front(Frenet_Coordinates);
        }
      }

      else // Border vertex can not be removed
      {
#ifdef DBG_Decimation_Conquest
        std::cout << "marker #13" << std::endl;
#endif
        Number_symbol++;
        this->InterConnectivity.push_front(4);

        put(this->Facet_Flag, face(h, _pMesh), CONQUERED);
        put(this->Vertex_Flag, target(next(h, _pMesh), _pMesh), CONQUERED);

        if((get(this->Vertex_Sign, target(h, _pMesh)) == PLUS) &&
           (get(this->Vertex_Sign, target(opposite(h, _pMesh), _pMesh)) ==
            MINUS))
        {
          if(get(this->Vertex_Sign, target(next(h, _pMesh), _pMesh)) == NOSIGN)
            put(this->Vertex_Sign, target(next(h, _pMesh), _pMesh), PLUS);
        }
        else if((get(this->Vertex_Sign, target(h, _pMesh)) == MINUS) &&
                (get(this->Vertex_Sign, target(opposite(h, _pMesh), _pMesh)) ==
                 PLUS))
        {
          if(get(this->Vertex_Sign, target(next(h, _pMesh), _pMesh)) == NOSIGN)
            put(this->Vertex_Sign, target(next(h, _pMesh), _pMesh), PLUS);
        }
        else if((get(this->Vertex_Sign, target(h, _pMesh)) == PLUS) &&
                (get(this->Vertex_Sign, target(opposite(h, _pMesh), _pMesh)) ==
                 PLUS))
        {
          if(get(this->Vertex_Sign, target(next(h, _pMesh), _pMesh)) == NOSIGN)
            put(this->Vertex_Sign, target(next(h, _pMesh), _pMesh), MINUS);
        }
        else if((get(this->Vertex_Sign, target(h, _pMesh)) == MINUS) &&
                (get(this->Vertex_Sign, target(opposite(h, _pMesh), _pMesh)) ==
                 MINUS))
        {
          if(get(this->Vertex_Sign, target(next(h, _pMesh), _pMesh)) == NOSIGN)
            put(this->Vertex_Sign, target(next(h, _pMesh), _pMesh), PLUS);
        }

        if(!CGAL::is_border_edge(next(h, _pMesh), _pMesh))
          Halfedges.push(opposite(next(h, _pMesh), _pMesh));

        if(!CGAL::is_border_edge(prev(h, _pMesh), _pMesh))
          Halfedges.push(opposite(prev(h, _pMesh), _pMesh));
      }
    }


    // if its front face sholud be labelled as NULL PATCH
    else
    {
#ifdef DBG_Decimation_Conquest
      std::cout << "marker #14" << std::endl;
#endif
      // Enter symbol 'NULL PATCH' into the list of symbols
      this->InterConnectivity.push_front(4);
      Number_symbol++;

      // its front face is tagged CONQUERED
      put(this->Facet_Flag, face(h, _pMesh), CONQUERED);
      put(this->Vertex_Flag,
          target(next(h, _pMesh), _pMesh),
          CONQUERED); ///////////////////////////

      if((get(this->Vertex_Sign, target(h, _pMesh)) == PLUS) &&
         (get(this->Vertex_Sign, target(opposite(h, _pMesh), _pMesh)) == MINUS))
      {
        if(get(this->Vertex_Sign, target(next(h, _pMesh), _pMesh)) == NOSIGN)
          put(this->Vertex_Sign, target(next(h, _pMesh), _pMesh), PLUS);
      }
      else if((get(this->Vertex_Sign, target(h, _pMesh)) == MINUS) &&
              (get(this->Vertex_Sign, target(opposite(h, _pMesh), _pMesh)) ==
               PLUS))
      {
        if(get(this->Vertex_Sign, target(next(h, _pMesh), _pMesh)) == NOSIGN)
          put(this->Vertex_Sign, target(next(h, _pMesh), _pMesh), PLUS);
      }
      else if((get(this->Vertex_Sign, target(h, _pMesh)) == PLUS) &&
              (get(this->Vertex_Sign, target(opposite(h, _pMesh), _pMesh)) ==
               PLUS))
      {
        if(get(this->Vertex_Sign, target(next(h, _pMesh), _pMesh)) == NOSIGN)
          put(this->Vertex_Sign, target(next(h, _pMesh), _pMesh), MINUS);
      }
      else if((get(this->Vertex_Sign, target(h, _pMesh)) == MINUS) &&
              (get(this->Vertex_Sign, target(opposite(h, _pMesh), _pMesh)) ==
               MINUS))
      {
        if(get(this->Vertex_Sign, target(next(h, _pMesh), _pMesh)) == NOSIGN)
          put(this->Vertex_Sign, target(next(h, _pMesh), _pMesh), PLUS);
      }

      // two other output gates are pused into the fifo queue.
      if(CGAL::is_border_edge(next(h, _pMesh), _pMesh) == false)
        Halfedges.push(opposite(next(h, _pMesh), _pMesh));

      if(CGAL::is_border_edge(prev(h, _pMesh), _pMesh) == false)
        Halfedges.push(opposite(prev(h, _pMesh), _pMesh));
    }
  }
  if((this->IsColored) && (!this->IsOneColor))
  {
#ifdef DBG_Decimation_Conquest
    std::cout << "marker #15" << std::endl;
#endif
    while(!this->InterVertexColor.empty())
    {
      Color_Unit Col = this->InterVertexColor.front();
      this->InterVertexColor.pop_front();
      this->VertexColor[Component_ID].push_front(Col);
    }
  }
  // InterConnectivity is the intermediate connectivity symbol container.
  // We put all symbols in the main container which is this->Connectivity
  while(!InterConnectivity.empty())
  {
#ifdef DBG_Decimation_Conquest
    std::cout << "marker #16" << std::endl;
#endif
    int Symbol = InterConnectivity.front();
    InterConnectivity.pop_front();
    this->Connectivity[Component_ID].push_front(Symbol);
  }

  // same operation than connectivity.
  while(!InterGeometry.empty())
  {
#ifdef DBG_Decimation_Conquest
    std::cout << "marker #17" << std::endl;
#endif
    Point_Int Geo = InterGeometry.front();
    InterGeometry.pop_front();
    this->Geometry[Component_ID].push_front(Geo);
  }

  this->NumberVertices[Component_ID].push_front(Number_vertices);
  this->NumberSymbol[Component_ID].push_front(Number_symbol);

  // if this decimation didn't remove any vertex, we should remove all
  // connectivity symbols. for that we store number of symbols.
  this->DumpSymbolDecimation = Number_symbol;


  return Number_vertices;
}


// Description : Decoding of the regulation conquest
template< typename HalfedgeGraph, typename PointMap, typename VertexColorMap >
void
Compression_Valence_Component< HalfedgeGraph, PointMap, VertexColorMap >::
    Un_Regulation(HalfedgeGraph &_pMesh,
                  Arithmetic_Codec &Decoder,
                  const int &Component_ID,
                  PointMap *_pm,
                  VertexColorMap *_v_cm)
{
  Init(_pMesh);

  Adaptive_Data_Model Connectivity(2);

  halfedge_iterator hi = halfedges(_pMesh).first;
  std::queue< halfedge_descriptor > Halfedges;

  unsigned Qbit =
      this->Qbit[Component_ID] + this->NumberChangeQuantization[Component_ID];

  int Alpha_range = Decoder.get_bits(Qbit + 1);
  int Alpha_offset = Decoder.get_bits(Qbit + 1);
  int Gamma_range = Decoder.get_bits(Qbit + 1);
  int Gamma_offset = Decoder.get_bits(Qbit + 1);

  if(this->Smallest_Alpha < 0)
    Alpha_offset = Alpha_offset + this->Smallest_Alpha;

  if(this->Smallest_Gamma < 0)
    Gamma_offset = Gamma_offset + this->Smallest_Gamma;

  bool check_alpha = false;
  bool check_gamma = false;

  if((Alpha_range == 0) || (Alpha_range == 1))
  {
    check_alpha = true;
    Alpha_range = 2;
  }

  if((Gamma_range == 0) || (Gamma_range == 1))
  {
    check_gamma = true;
    Gamma_range = 2;
  }

  float Color_step = 0.0;
  if(this->NumberColorQuantization[Component_ID] == 0)
    Color_step = this->Color_Quantization_Step;
  else
    Color_step = this->Color_Quantization_Step *
                 std::pow(2.0, this->NumberColorQuantization[Component_ID]);

  Adaptive_Data_Model alpha(Alpha_range);
  Adaptive_Data_Model gamma(Gamma_range);

  while(
      (get(this->vertex_Seed_Edge, target(*hi, _pMesh)) != 2 * Component_ID) ||
      (get(this->vertex_Seed_Edge, target(opposite(*hi, _pMesh), _pMesh)) !=
       2 * Component_ID + 1))
    hi++;

  put(this->Vertex_Flag, target(*hi, _pMesh), CONQUERED);
  put(this->Vertex_Flag, target(opposite(*hi, _pMesh), _pMesh), CONQUERED);
  Halfedges.push(*hi);

  halfedge_descriptor h;

  while(!Halfedges.empty())
  {
    h = Halfedges.front();
    Halfedges.pop();

    if((get(this->Facet_Flag, face(h, _pMesh)) == CONQUERED) ||
       (get(this->Facet_Flag, face(h, _pMesh)) ==
        TO_BE_REMOVED)) // already visited.
      continue;

    // read connectivity information
    int valence = Decoder.decode(Connectivity) + 3;

    // if valence is 3, remove the front vertex.
    if(valence == 3)
    {
      halfedge_descriptor g = h;

      // Insertion of a vertex
      halfedge_descriptor pass = h;

      std::vector< Point3d > Vertices; // contains 1-ring and 2-rings vertices
      Point3d Barycenter = Barycenter_Patch_After_Removal(pass, 3, _pMesh, _pm);
      Point_Int BC = Change_Real_Int(Barycenter, Component_ID);

      Vector normal = Triangle_Normal(_pMesh, _pm, pass);
      Vector T2 = Vector(0.0, 0.0, 0.0);
      Vector T1 = Calculate_T1_T2(pass, normal, T2, _pMesh, _pm);

      if(T1 == Vector(0.0, 0.0, 0.0))
      {
        T1 = Vector(1, 0, 0);
        T2 = Vector(0, 1, 0);
        normal = Vector(0, 0, 1);
      }

      else if(normal == Vector(0.0, 0.0, 0.0))
      {
        T1 = Vector(1, 0, 0);
        T2 = Vector(0, 1, 0);
        normal = Vector(0, 0, 1);
      }
      else if(T2 == Vector(0.0, 0.0, 0.0))
      {
        T1 = Vector(1, 0, 0);
        T2 = Vector(0, 1, 0);
        normal = Vector(0, 0, 1);
      }

      Point_Int Frenet;

      if(check_alpha == false)
      {
        Frenet.x = Decoder.decode(alpha);
        Frenet.y = Decoder.decode(alpha);
      }
      else
      {
        Frenet.x = 0;
        Frenet.y = 0;
      }

      if(check_gamma == false)
        Frenet.z = Decoder.decode(gamma);

      else
        Frenet.z = 0;

      Frenet.x -= Alpha_offset;
      Frenet.y -= Alpha_offset;
      Frenet.z -= Gamma_offset;

      Point_Int Diff;
      if(this->Is_Bijection_Enabled)
        Diff = Inverse_Frenet_Rotation(Frenet, T1, T2, normal);
      else
        Diff = Frenet;
      Point_Int Center = BC + Diff;

      Point3d Center_vertex = this->Change_Int_Real(Center, Component_ID);


      // Assign the region number to inserted vertex
      halfedge_descriptor reg = h;

      int Selected_region = 500000;
      // int Number_vertices = 500000;
      std::vector< int > T_Bin;
      std::vector< int > T_Number;

      for(int i = 0; i < valence; i++)
      {
        int N1 = get(this->vertex_Region_Number, target(reg, _pMesh));
        bool Is_existed = false;
        for(unsigned int j = 0; j < T_Bin.size(); j++)
        {
          if(N1 == T_Bin[j])
          {
            T_Number[j]++;
            Is_existed = true;
          }
        }
        if(!Is_existed)
        {
          T_Bin.push_back(N1);
          T_Number.push_back(1);
        }
        reg = next(reg, _pMesh);
      }
      int Max = -5000;
      for(unsigned int i = 0; i < T_Number.size(); i++)
      {
        if(T_Number[i] > Max)
          Max = T_Number[i];
      }
      std::vector< int > T_possible_bin;
      for(unsigned int i = 0; i < T_Number.size(); i++)
      {
        if(T_Number[i] == Max)
          T_possible_bin.push_back(T_Bin[i]);
      }

      if(T_possible_bin.size() == 1)
      {
        Selected_region = T_possible_bin[0];
      }
      else
      {
        Selected_region = 5000;
        for(unsigned int i = 0; i < T_possible_bin.size(); i++)
        {
          if(T_possible_bin[i] < Selected_region)
            Selected_region = T_possible_bin[i];
        }
      }


      // Vertex insertion
      g = CGAL::Euler::add_center_vertex(g, _pMesh);
      init_vertex_attributes(target(g, _pMesh));

      put(*_pm, target(g, _pMesh), Center_vertex);
      put(this->vertex_Seed_Edge, target(g, _pMesh), -1);

      int RO = this->GlobalCountOperation - Decompress_count;

      put(this->vertex_Removal_Order, target(g, _pMesh), RO);

      put(this->vertex_Region_Number, target(g, _pMesh), Selected_region);

      if(Selected_region != -1)
        this->m_Number_Vertices_Per_Regions[Selected_region]++;

      // Vertex flags
      put(this->Vertex_Flag, target(g, _pMesh), CONQUERED);
      g = h;
      put(this->Facet_Flag, face(g, _pMesh), CONQUERED);
      g = next(g, _pMesh);
      g = prev(opposite(g, _pMesh), _pMesh);
      put(this->Vertex_Flag, target(opposite(g, _pMesh), _pMesh), CONQUERED);
      put(this->Facet_Flag, face(g, _pMesh), CONQUERED);
      if(!CGAL::is_border_edge(prev(g, _pMesh), _pMesh))
      {
        halfedge_descriptor h1 = opposite(prev(g, _pMesh), _pMesh);
        put(this->Facet_Flag, face(h1, _pMesh), CONQUERED);
        put(this->Vertex_Flag, target(next(h1, _pMesh), _pMesh), CONQUERED);
        if(!CGAL::is_border_edge(next(h1, _pMesh), _pMesh))
          Halfedges.push(opposite(next(h1, _pMesh), _pMesh));
        if(!CGAL::is_border_edge(prev(h1, _pMesh), _pMesh))
          Halfedges.push(opposite(prev(h1, _pMesh), _pMesh));
      }
      g = prev(opposite(g, _pMesh), _pMesh);
      put(this->Facet_Flag, face(g, _pMesh), CONQUERED);
      put(this->Vertex_Flag, target(opposite(g, _pMesh), _pMesh), CONQUERED);
      if(!CGAL::is_border_edge(prev(g, _pMesh), _pMesh))
      {
        halfedge_descriptor h2 = opposite(prev(g, _pMesh), _pMesh);
        put(this->Facet_Flag, face(h2, _pMesh), CONQUERED);
        put(this->Vertex_Flag, target(next(h2, _pMesh), _pMesh), CONQUERED);
        if(!CGAL::is_border_edge(next(h2, _pMesh), _pMesh))
          Halfedges.push(opposite(next(h2, _pMesh), _pMesh));
        if(!CGAL::is_border_edge(prev(h2, _pMesh), _pMesh))
          Halfedges.push(opposite(prev(h2, _pMesh), _pMesh));
      }


      if((this->IsColored) && (!this->IsOneColor))
      {
        g = h;

        Color_Unit Predicted_color;

        Predicted_color = Get_Average_Vertex_Color_Lee(g, valence, _pMesh);

        Color_Unit Color_difference;
        Color_difference.c0 =
            this->Decoder.decode(this->Color_0_Model) + this->Smallest_C0;
        Color_difference.c1 =
            this->Decoder.decode(this->Color_1_Model) + this->Smallest_C1;
        Color_difference.c2 =
            this->Decoder.decode(this->Color_2_Model) + this->Smallest_C2;

        Color_Unit CV = Predicted_color + Color_difference;

        put(this->vertex_color_int,
            target(next(g, _pMesh), _pMesh),
            Color_Unit(CV.c0, CV.c1, CV.c2));

        float LAB[3];
        LAB[0] = this->C0_Min + CV.c0 * Color_step;
        LAB[1] = this->C1_Min + CV.c1 * Color_step;
        LAB[2] = this->C2_Min + CV.c2 * Color_step;

        float RGB[3];
        LAB_To_RGB(LAB[0], LAB[1], LAB[2], RGB);

        typedef
            typename boost::property_traits< VertexColorMap >::value_type Color;
        put(*_v_cm,
            target(next(g, _pMesh), _pMesh),
            Color(RGB[0], RGB[1], RGB[2]));
      }

      if((this->IsColored) && (this->IsOneColor))
      {
        g = next(h, _pMesh);
        typedef
            typename boost::property_traits< VertexColorMap >::value_type Color;
        put(*_v_cm,
            target(g, _pMesh),
            Color(this->OnlyColor[0], this->OnlyColor[1], this->OnlyColor[2]));
      }
    }

    else
    {
      put(this->Facet_Flag, face(h, _pMesh), CONQUERED);
      put(this->Vertex_Flag, target(next(h, _pMesh), _pMesh), CONQUERED);
      if(!CGAL::is_border_edge(next(h, _pMesh), _pMesh))
        Halfedges.push(opposite(next(h, _pMesh), _pMesh));
      if(!CGAL::is_border_edge(prev(h, _pMesh), _pMesh))
        Halfedges.push(opposite(prev(h, _pMesh), _pMesh));
    }
  }
}


//#define DBG_Un_Decimation_Conquest

// Description : Decoding function of decimation conquest
template< typename HalfedgeGraph, typename PointMap, typename VertexColorMap >
void
Compression_Valence_Component< HalfedgeGraph, PointMap, VertexColorMap >::
    Un_Decimation_Conquest(HalfedgeGraph &_pMesh,
                           Arithmetic_Codec &Decoder,
                           const int &Component_ID,
                           PointMap *_pm,
                           VertexColorMap *_v_cm)
{
#ifdef DBG_Un_Decimation_Conquest
  static unsigned int Un_Decimation_Conquest_call_cnt = 0;
  std::cout << "__func__"
            << "  Un_Decimation_Conquest_call_cnt = "
            << ++Un_Decimation_Conquest_call_cnt << std::endl;
#endif

  Init(_pMesh);

  int Number_connectivity_symbols;
  if(this->IsClosed[Component_ID])
    Number_connectivity_symbols = 5;
  else
    Number_connectivity_symbols = 7;

  Adaptive_Data_Model Connectivity(Number_connectivity_symbols);

  halfedge_iterator hi = halfedges(_pMesh).first;
  std::queue< halfedge_descriptor > Halfedges;

  unsigned Qbit =
      this->Qbit[Component_ID] + this->NumberChangeQuantization[Component_ID];

  int Alpha_range = Decoder.get_bits(Qbit + 1);
  int Alpha_offset = Decoder.get_bits(Qbit + 1);
  int Gamma_range = Decoder.get_bits(Qbit + 1);
  int Gamma_offset = Decoder.get_bits(Qbit + 1);

  if(this->Smallest_Alpha < 0)
    Alpha_offset = Alpha_offset + this->Smallest_Alpha;

  if(this->Smallest_Gamma < 0)
    Gamma_offset = Gamma_offset + this->Smallest_Gamma;

  bool check_alpha = false;
  bool check_gamma = false;

  if((Alpha_range == 0) || (Alpha_range == 1))
  {
    check_alpha = true;
    Alpha_range = 2;
  }

  if((Gamma_range == 0) || (Gamma_range == 1))
  {
    check_gamma = true;
    Gamma_range = 2;
  }


  float Color_step = 0.0;
  if(this->NumberColorQuantization[Component_ID] == 0)
    Color_step = this->Color_Quantization_Step;
  else
    Color_step = this->Color_Quantization_Step *
                 std::pow(2.0, this->NumberColorQuantization[Component_ID]);


  Adaptive_Data_Model alpha(Alpha_range);
  Adaptive_Data_Model gamma(Gamma_range);

  while(
      (get(this->vertex_Seed_Edge, target(*hi, _pMesh)) != 2 * Component_ID) ||
      (get(this->vertex_Seed_Edge, target(opposite(*hi, _pMesh), _pMesh)) !=
       2 * Component_ID + 1))
    hi++;

  // Two vertices of seed edges are flaged CONQUERED
  put(this->Vertex_Flag, target(*hi, _pMesh), CONQUERED);
  put(this->Vertex_Flag, target(opposite(*hi, _pMesh), _pMesh), CONQUERED);

  // These vertices are also flages with sign flags for retriangulation
  put(this->Vertex_Sign, target(*hi, _pMesh), PLUS);
  put(this->Vertex_Sign, target(opposite(*hi, _pMesh), _pMesh), MINUS);

  // Two vertices of seed edges are flaged CONQUERED
  Halfedges.push(*hi);

  halfedge_descriptor h;

#ifdef DBG_Un_Decimation_Conquest
  int dbg_loop_counter = 0;
#endif

  while(!Halfedges.empty())
  {
#ifdef DBG_Un_Decimation_Conquest
    std::cout << __func__ << "  marker #2"
              << "  loop " << ++dbg_loop_counter << std::endl;
    std::cout << __func__ << "  Halfedges.size()=" << Halfedges.size()
              << std::endl;
#endif

    h = Halfedges.front();
    Halfedges.pop();

#ifdef DBG_Un_Decimation_Conquest
    if(Un_Decimation_Conquest_call_cnt == 2)
      print_halfedge("h", h, _pMesh, _pm);
#endif

    unsigned int valence = 0, type = 0; // define type of retriangulation

    if(CGAL::is_border(h, _pMesh) == true)
      continue;

    // if its front face is not tagged CONQUERED nor TO_BE_REMOVED
    if((get(this->Facet_Flag, face(h, _pMesh)) == CONQUERED) ||
       (get(this->Facet_Flag, face(h, _pMesh)) == TO_BE_REMOVED))
      continue;

    valence = Decoder.decode(Connectivity) + 3;

    // if its front vertex is free
    if((valence >= 3) && (valence <= 6))
    {
      type = Find_Type(_pMesh, h, valence);

      // remove the front vertex if its removal does not viloate the manifold
      // property
      halfedge_descriptor pass = h;
      halfedge_descriptor g = h;

      Vector normal = Normal_Patch(pass, valence, _pMesh, _pm);
      Vector T2 = Vector(0.0, 0.0, 0.0);
      Vector T1 = Calculate_T1_T2(h, normal, T2, _pMesh, _pm);

      if(T1 == Vector(0.0, 0.0, 0.0))
      {
        T1 = Vector(1, 0, 0);
        T2 = Vector(0, 1, 0);
        normal = Vector(0, 0, 1);
      }

      else if(normal == Vector(0.0, 0.0, 0.0))
      {
        T1 = Vector(1, 0, 0);
        T2 = Vector(0, 1, 0);
        normal = Vector(0, 0, 1);
      }
      else if(T2 == Vector(0.0, 0.0, 0.0))
      {
        T1 = Vector(1, 0, 0);
        T2 = Vector(0, 1, 0);
        normal = Vector(0, 0, 1);
      }

      // remove edges to re_find polygon and (check. and attribute sign flag.)
      bool Check_Validity = Remove_Edges(_pMesh, g, type);
      Check_Validity = false;

      // Check_Validity = false;// * * * * * * * * * * *////
      if(Check_Validity == false)
      {
#ifdef DBG_Un_Decimation_Conquest
        if(Un_Decimation_Conquest_call_cnt == 9 && dbg_loop_counter == 1044)
          bool dbg_dummy = true; // just to be able to set a breakpoint
#endif
        g = h;
        halfedge_descriptor pass = h;

        std::vector< Point3d > Vertices; // contains 1-ring and 2-ring vertices;
        Point3d Barycenter =
            Barycenter_Patch_After_Removal(pass, valence, _pMesh, _pm);
        Point_Int BC = Change_Real_Int(Barycenter, Component_ID);

        Point_Int Frenet;
        if(check_alpha == false)
        {
          Frenet.x = Decoder.decode(alpha);
          Frenet.y = Decoder.decode(alpha);
        }
        else if(check_alpha == true)
        {
          Frenet.x = 0;
          Frenet.y = 0;
        }
        if(check_gamma == false)
          Frenet.z = Decoder.decode(gamma);
        else if(check_gamma == true)
          Frenet.z = 0;

        Frenet.x -= Alpha_offset;
        Frenet.y -= Alpha_offset;
        Frenet.z -= Gamma_offset;

        Point_Int Diff;
        if(this->Is_Bijection_Enabled)
          Diff = Inverse_Frenet_Rotation(Frenet, T1, T2, normal);
        else
          Diff = Frenet;
        /*Point_Int Center = BC + Diff;
        Point_Int Diff = Inverse_Frenet_Rotation(Frenet, T1, T2, normal);*/

        Point_Int Center = BC + Diff;

        Point3d Center_vertex = this->Change_Int_Real(Center, Component_ID);

        // Assign the region number to inserted vertex
        halfedge_descriptor reg = h;

        int Selected_region = 500000;
        // int Number_vertices = 500000;
        std::vector< int > T_Bin;
        std::vector< int > T_Number;

        for(int i = 0; i < (int)valence; i++)
        {
          int N1 = get(this->vertex_Region_Number, target(reg, _pMesh));
          bool Is_existed = false;
          for(unsigned int j = 0; j < T_Bin.size(); j++)
          {
            if(N1 == T_Bin[j])
            {
              T_Number[j]++;
              Is_existed = true;
            }
          }
          if(!Is_existed)
          {
            T_Bin.push_back(N1);
            T_Number.push_back(1);
          }
          reg = next(reg, _pMesh);
        }
        int Max = -5000;
        for(unsigned int i = 0; i < T_Number.size(); i++)
        {
          if(T_Number[i] > Max)
            Max = T_Number[i];
        }
        std::vector< int > T_possible_bin;
        for(unsigned int i = 0; i < T_Number.size(); i++)
        {
          if(T_Number[i] == Max)
            T_possible_bin.push_back(T_Bin[i]);
        }

        if(T_possible_bin.size() == 1)
        {
          Selected_region = T_possible_bin[0];
        }
        else
        {
          Selected_region = 5000;
          for(unsigned int i = 0; i < T_possible_bin.size(); i++)
          {
            if(T_possible_bin[i] < Selected_region)
              Selected_region = T_possible_bin[i];
          }
        }

        g = CGAL::Euler::add_center_vertex(g, _pMesh);
        init_vertex_attributes(target(g, _pMesh));
        put(*_pm, target(g, _pMesh), Center_vertex);
#ifdef DBG_Un_Decimation_Conquest
        std::cout << __func__ << "  marker PM#1"
                  << "  loop " << dbg_loop_counter << std::endl;
        std::cout << __func__
                  << "  vertex=" << vertex_to_string(target(g, _pMesh), _pm)
                  << std::endl;
#endif

        put(this->vertex_Region_Number, target(g, _pMesh), Selected_region);

        if(Selected_region != -1)
          this->m_Number_Vertices_Per_Regions[Selected_region]++;

        int RO = this->GlobalCountOperation - this->Decompress_count;
        put(this->vertex_Removal_Order, target(g, _pMesh), RO);

        put(this->Vertex_Flag, target(g, _pMesh), CONQUERED);

        g = h;
        put(this->Facet_Flag, face(g, _pMesh), TO_BE_REMOVED);
        put(this->Vertex_Flag, target(g, _pMesh), CONQUERED);
        put(this->Vertex_Flag, target(opposite(g, _pMesh), _pMesh), CONQUERED);

        for(unsigned int k = 0; k < (valence - 1); k++)
        {
          g = next(opposite(next(g, _pMesh), _pMesh), _pMesh);
          put(this->Facet_Flag, face(g, _pMesh), TO_BE_REMOVED);
          put(this->Vertex_Flag, target(g, _pMesh), CONQUERED);
          put(this->Vertex_Flag,
              target(opposite(g, _pMesh), _pMesh),
              CONQUERED);
          if(CGAL::is_border_edge(g, _pMesh) == false)
            Halfedges.push(opposite(g, _pMesh));
        }

        put(this->vertex_Seed_Edge, target(next(g, _pMesh), _pMesh), -1);

        if((this->IsColored) && (!this->IsOneColor))
        {
          g = h;

          Color_Unit Predicted_color;
          Predicted_color = Get_Average_Vertex_Color_Lee(g, valence, _pMesh);

          Color_Unit Color_difference;
          Color_difference.c0 =
              this->Decoder.decode(this->Color_0_Model) + this->Smallest_C0;
          Color_difference.c1 =
              this->Decoder.decode(this->Color_1_Model) + this->Smallest_C1;
          Color_difference.c2 =
              this->Decoder.decode(this->Color_2_Model) + this->Smallest_C2;

          Color_Unit CV = Predicted_color + Color_difference;

          put(this->vertex_color_int,
              target(next(g, _pMesh), _pMesh),
              Color_Unit(CV.c0, CV.c1, CV.c2));

          float LAB[3];
          LAB[0] = this->C0_Min + CV.c0 * Color_step;
          LAB[1] = this->C1_Min + CV.c1 * Color_step;
          LAB[2] = this->C2_Min + CV.c2 * Color_step;

          float RGB[3];
          LAB_To_RGB(LAB[0], LAB[1], LAB[2], RGB);

          typedef typename boost::property_traits< VertexColorMap >::value_type
              Color;
          put(*_v_cm,
              target(next(g, _pMesh), _pMesh),
              Color(RGB[0], RGB[1], RGB[2]));

#ifdef DBG_Un_Decimation_Conquest
          std::cout << __func__ << "  marker #3"
                    << "  loop " << dbg_loop_counter << std::endl;
          std::cout << __func__ << "  vertex="
                    << vertex_to_string(target(next(g, _pMesh), _pMesh), _pm)
                    << std::endl;
          std::cout << __func__ << "  Color=" << RGB[0] << " " << RGB[1] << " "
                    << RGB[2] << std::endl;
#endif

          //#endif
        }
        if((this->IsColored) && (this->IsOneColor))
        {
          g = next(h, _pMesh);
          typedef typename boost::property_traits< VertexColorMap >::value_type
              Color;
          put(*_v_cm,
              target(g, _pMesh),
              Color(
                  this->OnlyColor[0], this->OnlyColor[1], this->OnlyColor[2]));

#ifdef DBG_Un_Decimation_Conquest
          std::cout << __func__ << "  marker #4"
                    << "  loop " << dbg_loop_counter << std::endl;
          std::cout << __func__
                    << "  vertex=" << vertex_to_string(target(g, _pMesh), _pm)
                    << std::endl;
          std::cout << __func__ << "  Color=" << this->OnlyColor[0] << " "
                    << this->OnlyColor[1] << " " << this->OnlyColor[2]
                    << std::endl;
#endif
        }
      }
    }
    // In case of border edge.
    else if((valence == 8) || (valence == 9))
    {

      type = Find_Type(_pMesh, h, valence - 5);

      halfedge_descriptor pass = h;

      Vector normal = Normal_Patch(pass, valence - 5, _pMesh, _pm);
      Vector T2 = Vector(0.0, 0.0, 0.0);
      Vector T1 = Calculate_T1_T2(h, normal, T2, _pMesh, _pm);

      if(normal == Vector(0.0, 0.0, 0.0))
      {
        T1 = Vector(1, 0, 0);
        T2 = Vector(0, 1, 0);
        normal = Vector(0, 0, 1);
      }

      pass = h;
      Point_Int Frenet;
      if(check_alpha == false)
      {
        Frenet.x = Decoder.decode(alpha);
        Frenet.y = Decoder.decode(alpha);
      }
      else if(check_alpha == true)
      {
        Frenet.x = 0;
        Frenet.y = 0;
      }

      if(check_gamma == false)
        Frenet.z = Decoder.decode(gamma);
      else if(check_gamma == true)
        Frenet.z = 0;

      Frenet.x -= Alpha_offset;
      Frenet.y -= Alpha_offset;
      Frenet.z -= Gamma_offset;

      Point_Int Diff;
      if(this->Is_Bijection_Enabled)
        Diff = Inverse_Frenet_Rotation(Frenet, T1, T2, normal);
      else
        Diff = Frenet;

      // Point_Int Center = BC + Diff;

      // Point_Int Diff = Inverse_Frenet_Rotation(Frenet,T1,T2,normal);
      halfedge_descriptor g = h;

      Color_Unit Predicted_color;
      if((this->IsColored) && (!this->IsOneColor))
      {
        Predicted_color.c0 =
            Decoder.decode(this->Color_0_Model) + this->Smallest_C0;
        Predicted_color.c1 =
            Decoder.decode(this->Color_1_Model) + this->Smallest_C1;
        Predicted_color.c2 =
            Decoder.decode(this->Color_2_Model) + this->Smallest_C2;
      }

      // border edge with valence == 3
      if(valence == 8)
      {
        Point3d Barycenter =
            Barycenter_Patch_After_Removal(pass, 3, _pMesh, _pm);
        Point_Int BC = Change_Real_Int(Barycenter, Component_ID);

        Point_Int Center = BC + Diff;
        Point3d Center_vertex = this->Change_Int_Real(Center, Component_ID);

        //#ifdef PREDICTION_METHOD
        Color_Unit Average_color;
        if((this->IsColored) && (!this->IsOneColor))
        {
          Average_color = Get_Average_Vertex_Color_After_Removal(g, 3, _pMesh);
        }
        //#endif

        std::vector< halfedge_descriptor > Border_edges;
        int Number_jump = 0;

        if(CGAL::is_border_edge(next(g, _pMesh), _pMesh))
          Number_jump = 0;
        if(CGAL::is_border_edge(prev(g, _pMesh), _pMesh))
        {
          Number_jump = 1;
          g = next(g, _pMesh);
        }

        Border_edges.push_back(opposite(g, _pMesh));

        g = prev(g, _pMesh);
        Border_edges.push_back(opposite(g, _pMesh));

        CGAL::Euler::remove_face(g, _pMesh);

        halfedge_descriptor Prev_edge =
            prev(opposite(Border_edges[1], _pMesh), _pMesh);

        // g points the new vertex
        g = fixed_CGAL_Euler_add_vertex_and_face_to_border(
            Prev_edge, opposite(Border_edges[1], _pMesh), _pMesh);
        init_vertex_attributes(target(g, _pMesh));
        put(*_pm, target(g, _pMesh), Center_vertex);
#ifdef DBG_Un_Decimation_Conquest
        std::cout << __func__ << "  marker PM#2"
                  << "  loop " << dbg_loop_counter << std::endl;
        std::cout << __func__
                  << "  vertex=" << vertex_to_string(target(g, _pMesh), _pm)
                  << std::endl;
#endif
        put(this->Vertex_Flag, target(g, _pMesh), CONQUERED);

        //#ifdef PREDICTION_METHOD
        Color_Unit CV;
        if((this->IsColored) && (!this->IsOneColor))
        {
          CV = Average_color + Predicted_color;

          put(this->vertex_color_int,
              target(g, _pMesh),
              Color_Unit(CV.c0, CV.c1, CV.c2));

          float LAB[3];
          LAB[0] = this->C0_Min + CV.c0 * this->Color_Quantization_Step;
          LAB[1] = this->C1_Min + CV.c1 * this->Color_Quantization_Step;
          LAB[2] = this->C2_Min + CV.c2 * this->Color_Quantization_Step;

          float RGB[3];
          LAB_To_RGB(LAB[0], LAB[1], LAB[2], RGB);

          typedef typename boost::property_traits< VertexColorMap >::value_type
              Color;
          put(*_v_cm, target(g, _pMesh), Color(RGB[0], RGB[1], RGB[2]));

#ifdef DBG_Un_Decimation_Conquest
          std::cout << __func__ << "  marker #5"
                    << "  loop " << dbg_loop_counter << std::endl;
          std::cout << __func__
                    << "  vertex=" << vertex_to_string(target(g, _pMesh), _pm)
                    << std::endl;
          std::cout << __func__ << "  Color=" << RGB[0] << " " << RGB[1] << " "
                    << RGB[2] << std::endl;
#endif
        }
        //#endif

        if((this->IsColored) && (this->IsOneColor))
        {
          typedef typename boost::property_traits< VertexColorMap >::value_type
              Color;
          put(*_v_cm,
              target(g, _pMesh),
              Color(
                  this->OnlyColor[0], this->OnlyColor[1], this->OnlyColor[2]));

#ifdef DBG_Un_Decimation_Conquest
          std::cout << __func__ << "  marker #6"
                    << "  loop " << dbg_loop_counter << std::endl;
          std::cout << __func__
                    << "  vertex=" << vertex_to_string(target(g, _pMesh), _pm)
                    << std::endl;
          std::cout << __func__ << "  Color=" << this->OnlyColor[0] << " "
                    << this->OnlyColor[1] << " " << this->OnlyColor[2]
                    << std::endl;
#endif
        }


        Prev_edge = next(Prev_edge, _pMesh);
        CGAL::Euler::add_face_to_border(
            Prev_edge, opposite(Border_edges[0], _pMesh), _pMesh);

        halfedge_descriptor Tag_handle;

        if(Number_jump == 0)
          Tag_handle = Border_edges[1];
        if(Number_jump == 1)
          Tag_handle = opposite(Border_edges[0], _pMesh);
        put(this->Vertex_Flag, target(Tag_handle, _pMesh), CONQUERED);

        if((type == 1) || (type == 2) || (type == 4))
        {
          if(get(this->Vertex_Sign, target(Tag_handle, _pMesh)) == NOSIGN)
            put(this->Vertex_Sign, target(Tag_handle, _pMesh), PLUS);
        }
        else if(type == 3)
        {
          if(get(this->Vertex_Sign, target(Tag_handle, _pMesh)) == NOSIGN)
            put(this->Vertex_Sign, target(Tag_handle, _pMesh), MINUS);
        }
        for(int i = 0; i < 2; i++)
        {
          put(this->Facet_Flag,
              face(opposite(Border_edges[i], _pMesh), _pMesh),
              CONQUERED);
          if(i != Number_jump)
            Halfedges.push(Border_edges[i]);
        }
      }
      // border edge with valence == 4
      else if(valence == 9)
      {
        int Number_jump = -1;
        std::vector< halfedge_descriptor > Border_edges;

        if((type == 5) || (type == 8))
        {
          // jump == 0;
          if(CGAL::is_border_edge(next(g, _pMesh), _pMesh))
          {
            Number_jump = 0;

            Border_edges.push_back(opposite(g, _pMesh));

            Border_edges.push_back(opposite(
                prev(opposite(prev(g, _pMesh), _pMesh), _pMesh), _pMesh));
            Border_edges.push_back(opposite(
                next(opposite(prev(g, _pMesh), _pMesh), _pMesh), _pMesh));

            CGAL::Euler::remove_face(opposite(Border_edges[0], _pMesh), _pMesh);
            CGAL::Euler::remove_face(opposite(Border_edges[1], _pMesh), _pMesh);
          }

          // jump == 1;
          else if(CGAL::is_border_edge(
                      next(opposite(prev(g, _pMesh), _pMesh), _pMesh), _pMesh))
          {
            Number_jump = 1;

            Border_edges.push_back(opposite(next(g, _pMesh), _pMesh));
            Border_edges.push_back(opposite(g, _pMesh));
            Border_edges.push_back(opposite(
                prev(opposite(prev(g, _pMesh), _pMesh), _pMesh), _pMesh));

            CGAL::Euler::remove_face(opposite(Border_edges[2], _pMesh), _pMesh);
            CGAL::Euler::remove_face(opposite(Border_edges[0], _pMesh), _pMesh);
          }

          // jump == 2;
          else
          {
            Number_jump = 2;

            Border_edges.push_back(opposite(
                next(opposite(prev(g, _pMesh), _pMesh), _pMesh), _pMesh));
            Border_edges.push_back(opposite(next(g, _pMesh), _pMesh));
            Border_edges.push_back(opposite(g, _pMesh));

            CGAL::Euler::remove_face(opposite(Border_edges[0], _pMesh), _pMesh);
            CGAL::Euler::remove_face(opposite(Border_edges[1], _pMesh), _pMesh);
          }
        }
        else
        {
          if(CGAL::is_border_edge(prev(g, _pMesh), _pMesh))
          {
            Number_jump = 2;

            Border_edges.push_back(opposite(
                prev(opposite(next(g, _pMesh), _pMesh), _pMesh), _pMesh));
            Border_edges.push_back(opposite(
                next(opposite(next(g, _pMesh), _pMesh), _pMesh), _pMesh));
            Border_edges.push_back(opposite(g, _pMesh));

            CGAL::Euler::remove_face(opposite(Border_edges[2], _pMesh), _pMesh);
            CGAL::Euler::remove_face(opposite(Border_edges[1], _pMesh), _pMesh);
          }

          else if(CGAL::is_border_edge(
                      prev(opposite(next(g, _pMesh), _pMesh), _pMesh), _pMesh))
          {
            Number_jump = 1;

            Border_edges.push_back(opposite(
                next(opposite(next(g, _pMesh), _pMesh), _pMesh), _pMesh));
            Border_edges.push_back(opposite(g, _pMesh));
            Border_edges.push_back(opposite(prev(g, _pMesh), _pMesh));

            CGAL::Euler::remove_face(opposite(Border_edges[0], _pMesh), _pMesh);
            CGAL::Euler::remove_face(opposite(Border_edges[1], _pMesh), _pMesh);
          }
          else
          {
            Number_jump = 0;

            Border_edges.push_back(opposite(g, _pMesh));
            Border_edges.push_back(opposite(prev(g, _pMesh), _pMesh));
            Border_edges.push_back(opposite(
                prev(opposite(next(g, _pMesh), _pMesh), _pMesh), _pMesh));

            CGAL::Euler::remove_face(opposite(Border_edges[2], _pMesh), _pMesh);
            CGAL::Euler::remove_face(opposite(Border_edges[1], _pMesh), _pMesh);
          }
        }

        g = h;

        Point3d p0 = get(*_pm, target(Border_edges[0], _pMesh));
        Point3d p1 =
            get(*_pm, target(opposite(Border_edges[0], _pMesh), _pMesh));
        Point3d p2 = get(*_pm, target(Border_edges[1], _pMesh));
        Point3d p3 = get(*_pm, target(Border_edges[2], _pMesh));

        Point3d Barycenter((p0[0] + p1[0] + p2[0] + p3[0]) / 4.,
                           (p0[1] + p1[1] + p2[1] + p3[1]) / 4.,
                           (p0[2] + p1[2] + p2[2] + p3[2]) / 4.);

        Point_Int BC = Change_Real_Int(Barycenter, Component_ID);

        Point_Int Center = BC + Diff;

        Point3d Center_vertex = this->Change_Int_Real(Center, Component_ID);

        // to create the new facets
        halfedge_descriptor Prev_edge =
            prev(opposite(Border_edges[2], _pMesh), _pMesh);
        g = fixed_CGAL_Euler_add_vertex_and_face_to_border(
            Prev_edge, opposite(Border_edges[2], _pMesh), _pMesh);
        init_vertex_attributes(target(g, _pMesh));
        put(*_pm, target(g, _pMesh), Center_vertex);
#ifdef DBG_Un_Decimation_Conquest
        std::cout << __func__ << "  marker PM#3"
                  << "  loop " << dbg_loop_counter << std::endl;
        std::cout << __func__
                  << "  vertex=" << vertex_to_string(target(g, _pMesh), _pm)
                  << std::endl;
#endif

        //#ifdef PREDICTION_METHOD
        Color_Unit CV;
        if((this->IsColored) && (!this->IsOneColor))
        {
          Color_Unit Average_color;

          Color_Unit Color_0 =
              Get_Vertex_Color(opposite(Border_edges[0], _pMesh), _pMesh);
          Color_Unit Color_1 = Get_Vertex_Color(Border_edges[0], _pMesh);
          Color_Unit Color_2 = Get_Vertex_Color(Border_edges[1], _pMesh);
          Color_Unit Color_3 = Get_Vertex_Color(Border_edges[2], _pMesh);

          Average_color.c0 =
              (int)(Color_0.c0 + Color_1.c0 + Color_2.c0 + Color_3.c0) / 4.0;
          Average_color.c1 =
              (int)(Color_0.c1 + Color_1.c1 + Color_2.c1 + Color_3.c1) / 4.0;
          Average_color.c2 =
              (int)(Color_0.c2 + Color_1.c2 + Color_2.c2 + Color_3.c2) / 4.0;

          CV = Average_color + Predicted_color;

          put(this->vertex_color_int,
              target(g, _pMesh),
              Color_Unit(CV.c0, CV.c1, CV.c2));

          float LAB[3];
          LAB[0] = this->C0_Min + CV.c0 * this->Color_Quantization_Step;
          LAB[1] = this->C1_Min + CV.c1 * this->Color_Quantization_Step;
          LAB[2] = this->C2_Min + CV.c2 * this->Color_Quantization_Step;

          float RGB[3];
          LAB_To_RGB(LAB[0], LAB[1], LAB[2], RGB);

          typedef typename boost::property_traits< VertexColorMap >::value_type
              Color;
          put(*_v_cm, target(g, _pMesh), Color(RGB[0], RGB[1], RGB[2]));

#ifdef DBG_Un_Decimation_Conquest
          std::cout << __func__ << "  marker #7"
                    << "  loop " << dbg_loop_counter << std::endl;
          std::cout << __func__
                    << "  vertex=" << vertex_to_string(target(g, _pMesh), _pm)
                    << std::endl;
          std::cout << __func__ << "  Color=" << RGB[0] << " " << RGB[1] << " "
                    << RGB[2] << std::endl;
#endif
        }
        //#endif
        if((this->IsColored) && (this->IsOneColor))
        {
          typedef typename boost::property_traits< VertexColorMap >::value_type
              Color;
          put(*_v_cm,
              target(g, _pMesh),
              Color(
                  this->OnlyColor[0], this->OnlyColor[1], this->OnlyColor[2]));

#ifdef DBG_Un_Decimation_Conquest
          std::cout << __func__ << "  marker #8"
                    << "  loop " << dbg_loop_counter << std::endl;
          std::cout << __func__
                    << "  vertex=" << vertex_to_string(target(g, _pMesh), _pm)
                    << std::endl;
          std::cout << __func__ << "  Color=" << this->OnlyColor[0] << " "
                    << this->OnlyColor[1] << " " << this->OnlyColor[2]
                    << std::endl;
#endif
        }
        Prev_edge = next(Prev_edge, _pMesh);
        CGAL::Euler::add_face_to_border(
            Prev_edge, opposite(Border_edges[1], _pMesh), _pMesh);
        CGAL::Euler::add_face_to_border(
            Prev_edge, opposite(Border_edges[0], _pMesh), _pMesh);


        // vertex_tag
        if(Number_jump == 0)
        {
          if((type == 5) || (type == 8))
          {
            if(get(this->Vertex_Sign, target(Border_edges[2], _pMesh)) ==
               NOSIGN)
              put(this->Vertex_Sign, target(Border_edges[2], _pMesh), PLUS);
            if(get(this->Vertex_Sign,
                   target(opposite(Border_edges[2], _pMesh), _pMesh)) == NOSIGN)
              put(this->Vertex_Sign,
                  target(opposite(Border_edges[2], _pMesh), _pMesh),
                  MINUS);
          }
          else
          {
            if(get(this->Vertex_Sign, target(Border_edges[2], _pMesh)) ==
               NOSIGN)
              put(this->Vertex_Sign, target(Border_edges[2], _pMesh), MINUS);
            if(get(this->Vertex_Sign,
                   target(opposite(Border_edges[2], _pMesh), _pMesh)) == NOSIGN)
              put(this->Vertex_Sign,
                  target(opposite(Border_edges[2], _pMesh), _pMesh),
                  PLUS);
          }
        }

        else if(Number_jump == 1)
        {
          if((type == 5) || (type == 8))
          {
            if(get(this->Vertex_Sign, target(Border_edges[2], _pMesh)) ==
               NOSIGN)
              put(this->Vertex_Sign, target(Border_edges[2], _pMesh), MINUS);
            if(get(this->Vertex_Sign,
                   target(opposite(Border_edges[0], _pMesh), _pMesh)) == NOSIGN)
              put(this->Vertex_Sign,
                  target(opposite(Border_edges[0], _pMesh), _pMesh),
                  PLUS);
          }
          else
          {
            if(get(this->Vertex_Sign, target(Border_edges[2], _pMesh)) ==
               NOSIGN)
              put(this->Vertex_Sign, target(Border_edges[2], _pMesh), PLUS);
            if(get(this->Vertex_Sign,
                   target(opposite(Border_edges[0], _pMesh), _pMesh)) == NOSIGN)
              put(this->Vertex_Sign,
                  target(opposite(Border_edges[0], _pMesh), _pMesh),
                  MINUS);
          }
        }
        else // jump == 2
        {
          if((type == 5) || (type == 8))
          {
            if(get(this->Vertex_Sign, target(Border_edges[0], _pMesh)) ==
               NOSIGN)
              put(this->Vertex_Sign, target(Border_edges[0], _pMesh), PLUS);
            if(get(this->Vertex_Sign,
                   target(opposite(Border_edges[0], _pMesh), _pMesh)) == NOSIGN)
              put(this->Vertex_Sign,
                  target(opposite(Border_edges[0], _pMesh), _pMesh),
                  MINUS);
          }
          else
          {
            if(get(this->Vertex_Sign, target(Border_edges[0], _pMesh)) ==
               NOSIGN)
              put(this->Vertex_Sign, target(Border_edges[0], _pMesh), MINUS);
            if(get(this->Vertex_Sign,
                   target(opposite(Border_edges[0], _pMesh), _pMesh)) == NOSIGN)
              put(this->Vertex_Sign,
                  target(opposite(Border_edges[0], _pMesh), _pMesh),
                  PLUS);
          }
        }

        for(int i = 0; i < 3; i++)
        {
          put(this->Facet_Flag,
              face(opposite(Border_edges[i], _pMesh), _pMesh),
              CONQUERED);

          put(this->Vertex_Flag, target(Border_edges[i], _pMesh), CONQUERED);
          put(this->Vertex_Flag,
              target(opposite(Border_edges[i], _pMesh), _pMesh),
              CONQUERED);

          if(i != Number_jump)
            Halfedges.push(Border_edges[i]);
        }
      }
    }


    //  the symbol == N
    else if(valence == 7)
    {
      // its front face is tagged CONQUERED
      put(this->Facet_Flag, face(h, _pMesh), CONQUERED);
      put(this->Vertex_Flag, target(next(h, _pMesh), _pMesh), CONQUERED);

      if((get(this->Vertex_Sign, target(h, _pMesh)) == PLUS) &&
         (get(this->Vertex_Sign, target(opposite(h, _pMesh), _pMesh)) == MINUS))
      {
        if(get(this->Vertex_Sign, target(next(h, _pMesh), _pMesh)) == NOSIGN)
          put(this->Vertex_Sign, target(next(h, _pMesh), _pMesh), PLUS);
      }
      else if((get(this->Vertex_Sign, target(h, _pMesh)) == MINUS) &&
              (get(this->Vertex_Sign, target(opposite(h, _pMesh), _pMesh)) ==
               PLUS))
      {
        if(get(this->Vertex_Sign, target(next(h, _pMesh), _pMesh)) == NOSIGN)
          put(this->Vertex_Sign, target(next(h, _pMesh), _pMesh), PLUS);
      }
      else if((get(this->Vertex_Sign, target(h, _pMesh)) == PLUS) &&
              (get(this->Vertex_Sign, target(opposite(h, _pMesh), _pMesh)) ==
               PLUS))
      {
        if(get(this->Vertex_Sign, target(next(h, _pMesh), _pMesh)) == NOSIGN)
          put(this->Vertex_Sign, target(next(h, _pMesh), _pMesh), MINUS);
      }
      else if((get(this->Vertex_Sign, target(h, _pMesh)) == MINUS) &&
              (get(this->Vertex_Sign, target(opposite(h, _pMesh), _pMesh)) ==
               MINUS))
      {
        if(get(this->Vertex_Sign, target(next(h, _pMesh), _pMesh)) == NOSIGN)
          put(this->Vertex_Sign, target(next(h, _pMesh), _pMesh), PLUS);
      }

      if(CGAL::is_border_edge(next(h, _pMesh), _pMesh) == false)
        Halfedges.push(opposite(next(h, _pMesh), _pMesh));

      if(CGAL::is_border_edge(prev(h, _pMesh), _pMesh) == false)
        Halfedges.push(opposite(prev(h, _pMesh), _pMesh));
    }
  }
}


// When the last loop doen not remove any vertex, all information are deleted
template< typename HalfedgeGraph, typename PointMap, typename VertexColorMap >
void
Compression_Valence_Component< HalfedgeGraph, PointMap, VertexColorMap >::
    Remove_Last_Phase_Elements(const int &Component_ID)
{
  this->NumberDecimation[Component_ID] -=
      1; // this->NumberDecimation[Component_ID] - 1;
  this->ComponentOperations[Component_ID] -= 1;
  this->ListOperation[Component_ID].pop_front();

  for(int i = 0; i < (this->DumpSymbolDecimation + this->DumpSymbolRegulation);
      i++)
    this->Connectivity[Component_ID].pop_front();

  for(int i = 0; i < 2; i++)
  {
    this->NumberSymbol[Component_ID].pop_front();
    this->NumberVertices[Component_ID].pop_front();
  }
}


// To store information needed for the reconstruction of the base mesh.
template< typename HalfedgeGraph, typename PointMap, typename VertexColorMap >
void
Compression_Valence_Component< HalfedgeGraph, PointMap, VertexColorMap >::
    Write_Base_Mesh(const HalfedgeGraph &_pMesh,
                    Arithmetic_Codec &enc,
                    unsigned &Connectivity_size,
                    unsigned &Color_size,
                    const int &Num_color_base_mesh,
                    const PointMap *_pm)
{

  unsigned int Max_Qbit = 0;
  for(int i = 0; i < this->NumberComponents; i++)
  {
    enc.put_bits(this->ComponentOperations[i], 8); // number of operations < 256
    enc.put_bits(this->Qbit[i] - 4, 4);            // quantization bit < 16
    enc.put_bits(this->NumberChangeQuantization[i],
                 4); // number of decrease of quantization resolution
    enc.put_bits(this->NumberColorQuantization[i], 3);

    if(this->Qbit[i] > Max_Qbit)
      Max_Qbit = this->Qbit[i];
  }

  enc.put_bits(FEVV::size_of_vertices(_pMesh),
               15); // number of vertices of base mesh < 4096
  enc.put_bits(FEVV::size_of_faces(_pMesh),
               16); // number of facets of base mesh < 8192

  // int Base_color_index_bit = 0;

  int Basemesh_vertex_number = 0;
  std::vector< int > Seed_edges(2 * this->NumberComponents, -1);


  // Encoding of vertex information of base mesh //
  for(vertex_iterator pVertex = vertices(_pMesh).first;
      pVertex != vertices(_pMesh).second;
      Basemesh_vertex_number++, pVertex++)
  {
    put(this->Vertex_Number, *pVertex, Basemesh_vertex_number);

    // int SEDD = pVertex->Seed_Edge;

    if(get(this->vertex_Seed_Edge, *pVertex) != OTHER_COORDINATE)
      Seed_edges[get(this->vertex_Seed_Edge, *pVertex)] =
          Basemesh_vertex_number;

    // int cid = pVertex->Component_Number;

    Point_Int Vertex = Change_Real_Int(
        get(*_pm, *pVertex), get(this->vertex_Component_Number, *pVertex));


    enc.put_bits(Vertex.x, Max_Qbit + 1);
    enc.put_bits(Vertex.y, Max_Qbit + 1);
    enc.put_bits(Vertex.z, Max_Qbit + 1);

    if((this->IsColored) && (!this->IsOneColor))
    {
      //#ifdef PREDICTION_METHOD
      Color_Unit color = get(this->vertex_color_int, *pVertex);
      int C0 = color.c0;
      int C1 = color.c1;
      int C2 = color.c2;

      enc.put_bits(C0, C0_QUANTIZATION);
      enc.put_bits(C1, C1_QUANTIZATION);
      enc.put_bits(C2, C2_QUANTIZATION);

      Color_size += 3 * C0_QUANTIZATION;
      //#endif
    }
  }

  // Bits needed for each edge.
  int Facet_index_bit = (int)ceil(
      log((double)(FEVV::size_of_vertices(_pMesh) + 1)) / log((double)2));

  int Count_facet_index = 0;
  for(face_iterator pFacet = faces(_pMesh).first;
      pFacet != faces(_pMesh).second;
      pFacet++)
  {
    Count_facet_index++;

    halfedge_descriptor pHalfedge = halfedge(*pFacet, _pMesh);
    do
    {
      enc.put_bits(get(this->Vertex_Number, target(pHalfedge, _pMesh)),
                   Facet_index_bit);
      pHalfedge = next(pHalfedge, _pMesh);


    } while(pHalfedge != halfedge(*pFacet, _pMesh));
  }

  // Store seed edge information.
  for(int i = 0; i < (int)Seed_edges.size(); i++) // MT
    enc.put_bits(Seed_edges[i], Facet_index_bit);

  Connectivity_size +=
      Facet_index_bit * (Count_facet_index + Seed_edges.size());
}


template< typename HalfedgeGraph, typename PointMap, typename VertexColorMap >
void
Compression_Valence_Component< HalfedgeGraph, PointMap, VertexColorMap >::
    Calculate_Geometry_Color_Offset_Range(void)
{
  for(int Component_ID = 0; Component_ID < this->NumberComponents;
      Component_ID++)
  {
    std::list< int >::iterator Number_iterator =
        this->NumberVertices[Component_ID].begin();
    std::list< Point_Int >::iterator Vertex_iterator =
        this->Geometry[Component_ID].begin();

    unsigned Number_phases = this->NumberVertices[Component_ID].size();


    for(unsigned i = 0; i < Number_phases; i++)
    {
      int Number_vertices_layer = *Number_iterator;
      Number_iterator++;

      int alpha_max = -10000;
      int alpha_min = 10000;
      int gamma_max = -10000;
      int gamma_min = 10000;
      int alpha = 0, beta = 0, gamma = 0;


      if(Number_vertices_layer != 0)
      {
        for(int j = 0; j < Number_vertices_layer; j++)
        {
          alpha = Vertex_iterator->x;
          if(alpha > alpha_max)
            alpha_max = alpha;
          if(alpha < alpha_min)
            alpha_min = alpha;

          beta = Vertex_iterator->y;
          if(beta > alpha_max)
            alpha_max = beta;
          if(beta < alpha_min)
            alpha_min = beta;

          gamma = Vertex_iterator->z;
          if(gamma > gamma_max)
            gamma_max = gamma;
          if(gamma < gamma_min)
            gamma_min = gamma;

          Vertex_iterator++;
        }

        this->AlphaRange[Component_ID].push_back(alpha_max - alpha_min + 1);
        this->AlphaOffset[Component_ID].push_back(-alpha_min);

        this->GammaRange[Component_ID].push_back(gamma_max - gamma_min + 1);
        this->GammaOffset[Component_ID].push_back(-gamma_min);
      }
      else
      {
        this->AlphaRange[Component_ID].push_back(0);
        this->AlphaOffset[Component_ID].push_back(0);

        this->GammaRange[Component_ID].push_back(0);
        this->GammaOffset[Component_ID].push_back(0);
      }
    }
  }


  /* Calculate alpha_min and gamma_min of all coefficients
   * in order to prevent negative symbols.Compression is not possible. */
  std::list< int >::iterator it_gamma, it_alpha;
  for(int Component_ID = 0; Component_ID < this->NumberComponents;
      Component_ID++)
  {
    for(it_alpha = this->AlphaOffset[Component_ID].begin();
        it_alpha != this->AlphaOffset[Component_ID].end();
        it_alpha++)
    {
      if(*it_alpha < this->Smallest_Alpha)
        this->Smallest_Alpha = *it_alpha;
    }
    for(it_gamma = this->GammaOffset[Component_ID].begin();
        it_gamma != this->GammaOffset[Component_ID].end();
        it_gamma++)
    {
      if(*it_gamma < this->Smallest_Gamma)
        this->Smallest_Gamma = *it_gamma;
    }
  }


  if(this->IsColored)
  {

    int C0_min = 50000, C1_min = 50000, C2_min = 50000;
    int C0_max = -50000, C1_max = -50000, C2_max = -50000;

    for(int Component_ID = 0; Component_ID < this->NumberComponents;
        Component_ID++)
    {
      //#ifdef PREDICTION_METHOD
      std::list< Color_Unit >::iterator Vertex_color_iterator;
      for(Vertex_color_iterator = this->VertexColor[Component_ID].begin();
          Vertex_color_iterator != this->VertexColor[Component_ID].end();
          Vertex_color_iterator++)
      //#endif
      {
        if(Vertex_color_iterator->c0 < C0_min)
          C0_min = Vertex_color_iterator->c0;
        if(Vertex_color_iterator->c0 > C0_max)
          C0_max = Vertex_color_iterator->c0;

        if(Vertex_color_iterator->c1 < C1_min)
          C1_min = Vertex_color_iterator->c1;
        if(Vertex_color_iterator->c1 > C1_max)
          C1_max = Vertex_color_iterator->c1;

        if(Vertex_color_iterator->c2 < C2_min)
          C2_min = Vertex_color_iterator->c2;
        if(Vertex_color_iterator->c2 > C2_max)
          C2_max = Vertex_color_iterator->c2;
      }
    }


    this->C0_Range = C0_max - C0_min + 1;
    this->C1_Range = C1_max - C1_min + 1;
    this->C2_Range = C2_max - C2_min + 1;

    if(this->C0_Range <= 1)
      this->C0_Range = 2;
    if(this->C1_Range <= 1)
      this->C1_Range = 2;
    if(this->C2_Range <= 1)
      this->C2_Range = 2;

    this->Smallest_C0 = C0_min;
    this->Smallest_C1 = C1_min;
    this->Smallest_C2 = C2_min;
  }
}


//#define DBG_Simplification

template< typename HalfedgeGraph, typename PointMap, typename VertexColorMap >
void
Compression_Valence_Component< HalfedgeGraph, PointMap, VertexColorMap >::
    Simplification(HalfedgeGraph &_pMesh,
                   const PointMap *_pm,
                   const int &NVertices,
                   const bool Normal_flipping,
                   const bool Use_metric,
                   const float &Metric_thread,
                   const bool Use_forget_metric,
                   const int &Forget_value)
{


  // bool Is_any_vertex_removed = true;

  unsigned Last_Number = 0;
  unsigned Current_Number = FEVV::size_of_vertices(_pMesh);

  // int Operation_choice = -1;

  do
  {
    Last_Number = Current_Number;

    // Simplify component by component
    for(int Component_ID = 0; Component_ID < this->NumberComponents;
        Component_ID++)
    {
      // if the ith component did not remove any vertex in last loop, it is not
      // necessary to simplifiy it.
      if(this->ComponentOperations[Component_ID] == this->GlobalCountOperation)
      {
        unsigned Initial_number_vertices = FEVV::size_of_vertices(_pMesh);

#ifdef DBG_Simplification
        std::cout << "in " << __func__
                  << ", before Decimation_Conquest, num_vertices(_pMesh) = "
                  << FEVV::size_of_vertices(_pMesh) << std::endl;
#endif

        this->Decimation_Conquest(_pMesh,
                                  _pm,
                                  Normal_flipping,
                                  Use_metric,
                                  Metric_thread,
                                  Use_forget_metric,
                                  Forget_value,
                                  Component_ID);

#ifdef DBG_Simplification
        std::cout << "in " << __func__
                  << ", after Decimation_Conquest, num_vertices(_pMesh) = "
                  << FEVV::size_of_vertices(_pMesh) << std::endl;
#endif

        this->Regulation(_pMesh,
                         Normal_flipping,
                         Use_metric,
                         Metric_thread,
                         Use_forget_metric,
                         Forget_value,
                         Component_ID,
                         _pm);

        int Diff_number_vertices =
            FEVV::size_of_vertices(_pMesh) - Initial_number_vertices;

        this->ComponentOperations[Component_ID] += 1;
        this->NumberDecimation[Component_ID] += 1;

        this->ListOperation[Component_ID].push_front(0);

        if(Diff_number_vertices == 0)
          this->Remove_Last_Phase_Elements(Component_ID);
      }
    }

    Current_Number = FEVV::size_of_vertices(_pMesh);

#ifdef DBG_Simplification
    // DBG-ELO+beg
    {
      // count real vertices number
      vertex_iterator v_begin = vertices(_pMesh).first;
      vertex_iterator v_end = vertices(_pMesh).second;
      unsigned int real_vertex_count = 0;
      for(vertex_iterator vi = v_begin; vi != v_end; ++vi)
        real_vertex_count++;
      assert((Current_Number == real_vertex_count));
    }
    // DBG-ELO+end
#endif

    if(Current_Number != Last_Number)
      this->GlobalCountOperation++;

    if(Current_Number < (unsigned)NVertices) // MT
      break;

  } while((Current_Number != Last_Number));

  compute_normals(_pMesh, _pm);
}


template< typename HalfedgeGraph, typename PointMap, typename VertexColorMap >
void
Compression_Valence_Component< HalfedgeGraph, PointMap, VertexColorMap >::
    Compression(HalfedgeGraph &_pMesh,
                const char *File_Name,
                const int &_Qbit,
                unsigned &Connectivity_size,
                unsigned &Color_size,
                unsigned &Total_size,
                const PointMap *_pm)
// const unsigned & Initial_file_size)
{
  // Calculate offset and range for the compression.
  this->Calculate_Geometry_Color_Offset_Range();

  FILE *fp =
      fopen(File_Name, "wb"); // Main FILE to save compression information.

  int res;

  res = fwrite(&this->Smallest_Alpha,
               sizeof(int),
               1,
               fp); // smallest value of alpha (to save the negative value)
  res = fwrite(&this->Smallest_Gamma,
               sizeof(int),
               1,
               fp); // smallest value of gamma (to save the negative value)
  res = fwrite(
      &this->Initial_file_size,
      sizeof(unsigned),
      1,
      fp); // Intial size of the input file (To visualize during decompression)
  res = fwrite(&this->NumberComponents, sizeof(int), 1, fp);

  for(int i = 0; i < this->NumberComponents; i++)
  {
    res = fwrite(&this->Quantization_Step[i],
                 sizeof(float),
                 1,
                 fp); // Quantization_Step(step of quantization)
    res = fwrite(&this->xmin[i], sizeof(float), 1, fp); // xmin value
    res = fwrite(&this->ymin[i], sizeof(float), 1, fp); // ymin value
    res = fwrite(&this->zmin[i], sizeof(float), 1, fp); // zmin value
  }


  // Type of mesh
  char Colored;
  char OneColor;

  if(this->IsColored)
  {
    Colored = 1;
    if(this->IsOneColor)
      OneColor = 1;
    else
      OneColor = 0;
  }
  else
    Colored = 0;

  res = fwrite(&Colored, sizeof(char), 1, fp);
  if(this->IsColored)
    res = fwrite(&OneColor, sizeof(char), 1, fp);


  int Num_color_base_mesh = 0;
  if((this->IsColored) && (!this->IsOneColor))
  {

    res = fwrite(&this->Color_Quantization_Step, sizeof(float), 1, fp);

    // En-tete pour la couleur
    res = fwrite(&this->C0_Min, sizeof(float), 1, fp); // smallest value of c0
    res = fwrite(&this->C1_Min, sizeof(float), 1, fp); // smallest value of c1
    res = fwrite(&this->C2_Min, sizeof(float), 1, fp); // smallest value of c2

    res = fwrite(&this->Smallest_C0, sizeof(int), 1, fp);
    res = fwrite(&this->Smallest_C1, sizeof(int), 1, fp);
    res = fwrite(&this->Smallest_C2, sizeof(int), 1, fp);

    Color_size += sizeof(int) * 8 * 9;
  }

  if((this->IsColored) && (this->IsOneColor))
  {

    res = fwrite(
        &this->OnlyColor[0], sizeof(float), 1, fp); // smallest value of c0
    res = fwrite(
        &this->OnlyColor[1], sizeof(float), 1, fp); // smallest value of c1
    res = fwrite(
        &this->OnlyColor[2], sizeof(float), 1, fp); // smallest value of c2
  }

  res += 0; // just to avoid gcc 4.6 warning

  // Declaration du codeur.
  Arithmetic_Codec enc(AC_BUFFER);
  enc.start_encoder();

  // To calculate connectivity rate
  Arithmetic_Codec Connectivity_encoder(AC_BUFFER);
  Connectivity_encoder.start_encoder();

  for(int i = 0; i < this->NumberComponents; i++)
  {
    if(this->IsClosed[i])
      enc.put_bits(0, 1);

    else
      enc.put_bits(1, 1);
  }

  if(this->Is_Bijection_Enabled)
    enc.put_bits(1, 1);
  else
    enc.put_bits(0, 1);

  /*	Write information of base mesh.
          geometry + connectivity + color information(if the mesh is colored).
   */
  this->Write_Base_Mesh(
      _pMesh, enc, Connectivity_size, Color_size, Num_color_base_mesh, _pm);

  // To calculate color rate
  Arithmetic_Codec Color_enc(AC_BUFFER);
  Color_enc.start_encoder();


  //#ifdef PREDICTION_METHOD

  Adaptive_Data_Model C0_model;
  Adaptive_Data_Model C1_model;
  Adaptive_Data_Model C2_model;

  // To calculate color rate
  Adaptive_Data_Model PC_C0_model;
  Adaptive_Data_Model PC_C1_model;
  Adaptive_Data_Model PC_C2_model;

  if((this->IsColored) && (!this->IsOneColor))
  {
    enc.put_bits(this->C0_Range, C0_QUANTIZATION + 1);
    enc.put_bits(this->C1_Range, C1_QUANTIZATION + 1);
    enc.put_bits(this->C2_Range, C2_QUANTIZATION + 1);

    C0_model.set_alphabet(this->C0_Range);
    C1_model.set_alphabet(this->C1_Range);
    C2_model.set_alphabet(this->C2_Range);

    // To calculate color rate
    Color_enc.put_bits(this->C0_Range, C0_QUANTIZATION + 1);
    Color_enc.put_bits(this->C1_Range, C1_QUANTIZATION + 1);
    Color_enc.put_bits(this->C2_Range, C2_QUANTIZATION + 1);

    PC_C0_model.set_alphabet(this->C0_Range);
    PC_C1_model.set_alphabet(this->C1_Range);
    PC_C2_model.set_alphabet(this->C2_Range);
  }

  this->DM_JCW_MOVE_ERROR.set_alphabet(3);

  // Main loop of compression //
  for(int i = 0; i < this->GlobalCountOperation; i++)
  {

    for(int Component_ID = 0; Component_ID < this->NumberComponents;
        Component_ID++)
    {
      if(i < this->ComponentOperations[Component_ID])
      {
        int Number_connectivity_symbols;

        if(this->IsClosed[Component_ID])
          Number_connectivity_symbols = 5;
        else
          Number_connectivity_symbols = 7;

        int Type_operation = this->ListOperation[Component_ID].front();
        this->ListOperation[Component_ID].pop_front();

        // decimation is chosen to be applied.
        if(Type_operation == 0)
        {
          enc.put_bits(0, 2);

          for(int j = 0; j < 2; j++) // Decimation and regulation
          {
            Adaptive_Data_Model Connectivity;

            if(j == 0)
              Connectivity.set_alphabet(2);
            else
              Connectivity.set_alphabet(Number_connectivity_symbols);

            Adaptive_Data_Model Temp_connectivity;
            if(j == 0)
              Temp_connectivity.set_alphabet(2);
            else
              Temp_connectivity.set_alphabet(Number_connectivity_symbols);

            int Alpha_range = this->AlphaRange[Component_ID].front();
            this->AlphaRange[Component_ID].pop_front();
            int Alpha_offset = this->AlphaOffset[Component_ID].front();
            this->AlphaOffset[Component_ID].pop_front();

            int Gamma_range = this->GammaRange[Component_ID].front();
            this->GammaRange[Component_ID].pop_front();
            int Gamma_offset = this->GammaOffset[Component_ID].front();
            this->GammaOffset[Component_ID].pop_front();

            enc.put_bits(Alpha_range, _Qbit + 1);
            if(this->Smallest_Alpha < 0)
              enc.put_bits(Alpha_offset - this->Smallest_Alpha, _Qbit + 1);
            else
              enc.put_bits(Alpha_offset, _Qbit + 1);

            enc.put_bits(Gamma_range, _Qbit + 1);
            if(this->Smallest_Gamma < 0)
              enc.put_bits(Gamma_offset - this->Smallest_Gamma, _Qbit + 1);
            else
              enc.put_bits(Gamma_offset, _Qbit + 1);

            bool check_alpha = false;
            bool check_gamma = false;

            if((Alpha_range == 0) || (Alpha_range == 1))
            {
              check_alpha = true;
              Alpha_range = 2;
            }
            if((Gamma_range == 0) || (Gamma_range == 1))
            {
              check_gamma = true;
              Gamma_range = 2;
            }

            Adaptive_Data_Model alpha(Alpha_range);
            Adaptive_Data_Model gamma(Gamma_range);

            int Number_symbols = this->NumberSymbol[Component_ID].front();
            this->NumberSymbol[Component_ID].pop_front();

            for(unsigned k = 0; k < (unsigned)Number_symbols; k++) // MT
            {
              unsigned symbol = this->Connectivity[Component_ID].front();
              this->Connectivity[Component_ID].pop_front();

              enc.encode(symbol, Connectivity);

              // To calculare connectivity rate
              Connectivity_encoder.encode(symbol, Temp_connectivity);

              if(((j == 0) && (symbol != 1)) || ((j == 1) && (symbol != 4)))
              {
                Point_Int Coeff = this->Geometry[Component_ID].front();
                this->Geometry[Component_ID].pop_front();

                int x = Coeff.x + Alpha_offset;
                int y = Coeff.y + Alpha_offset;
                int z = Coeff.z + Gamma_offset;

                if(check_alpha == false)
                {
                  enc.encode(x, alpha);
                  enc.encode(y, alpha);
                }
                if(check_gamma == false)
                {
                  enc.encode(z, gamma);
                }


                //#ifdef PREDICTION_METHOD
                if((this->IsColored) && (!this->IsOneColor))
                {
                  Color_Unit VC = this->VertexColor[Component_ID].front();
                  this->VertexColor[Component_ID].pop_front();

                  enc.encode(VC.c0 - this->Smallest_C0, C0_model);
                  enc.encode(VC.c1 - this->Smallest_C1, C1_model);
                  enc.encode(VC.c2 - this->Smallest_C2, C2_model);

                  // To calculate color rate
                  Color_enc.encode(VC.c0 - this->Smallest_C0, PC_C0_model);
                  Color_enc.encode(VC.c1 - this->Smallest_C1, PC_C1_model);
                  Color_enc.encode(VC.c2 - this->Smallest_C2, PC_C2_model);
                }
                //#endif
              }
            }
            alpha.reset();
            gamma.reset();
          }

          if(!this->m_N_Errors.empty())
          {
            int tt_count = this->m_N_Errors.front();
            this->m_N_Errors.pop_front();

            // if(!this->m_JCW_Move_Error.empty())
            for(int tt_i = 0; tt_i < tt_count; tt_i++)
            {
              std::vector< int > Temp_error = this->m_JCW_Move_Error.front();
              this->m_JCW_Move_Error.pop_front();

              for(int m = 0; m < 3; m++)
              {
                int Te = Temp_error[m];

                if(Te == -1)
                  Te = 2;

                enc.encode(Te, this->DM_JCW_MOVE_ERROR);
              }
            }
          }
        }

        // Decrease of geometry quantization resolution.
        else if(Type_operation == 1)
        {
          enc.put_bits(1, 2);

          int Number_vertices =
              this->NumberQuantizationLayer[Component_ID].front();
          this->NumberQuantizationLayer[Component_ID].pop_front();

          Adaptive_Data_Model Under_quantization_model(8);

          for(int i = 0; i < Number_vertices; i++)
          {
            int Under_quantization_coeff =
                this->QuantizationCorrectVector[Component_ID].front();
            this->QuantizationCorrectVector[Component_ID].pop_front();

            enc.encode(Under_quantization_coeff, Under_quantization_model);
          }
        }

        // Decrease of color quantization resolution.
        else
        {
          enc.put_bits(2, 2);
          Arithmetic_Codec Color_enc(AC_BUFFER);
          Color_enc.start_encoder();

          int Number_vertices =
              this->NumberProcessedVertices[Component_ID].front();
          this->NumberProcessedVertices[Component_ID].pop_front();

          Adaptive_Data_Model *Color_quantization_model =
              new Adaptive_Data_Model[COLOR_NUMBER];
          for(int i = 0; i < COLOR_NUMBER; i++)
            Color_quantization_model[i].set_alphabet(8);

          // To measure color rates
          Adaptive_Data_Model *Temp_quantization_model =
              new Adaptive_Data_Model[COLOR_NUMBER];
          for(int i = 0; i < COLOR_NUMBER; i++)
            Temp_quantization_model[i].set_alphabet(8);

          for(int i = 0; i < Number_vertices; i++)
          {
            int Color_index = this->ColorEncoderIndex[Component_ID].front();
            this->ColorEncoderIndex[Component_ID].pop_front();

            int Childcell_index =
                this->ColorChildcellIndex[Component_ID].front();
            this->ColorChildcellIndex[Component_ID].pop_front();

            enc.encode(Childcell_index, Color_quantization_model[Color_index]);

            Color_enc.encode(Childcell_index,
                             Temp_quantization_model[Color_index]);
          }
          Color_size += Color_enc.stop_encoder() * 8;

          delete[] Color_quantization_model;
          delete[] Temp_quantization_model;
        }
      }
    }
  }

  Connectivity_size += Connectivity_encoder.stop_encoder() * 8;
  Color_size += Color_enc.stop_encoder() * 8;

  enc.write_to_file(fp);
  fclose(fp);

  FILE *f_size = fopen(File_Name, "rb");
  fseek(f_size, 0, SEEK_END);
  Total_size = ftell(f_size);
}


//#define DBG_Calculate_Edge_Color_Difference
//#define DBG_Calculate_Edge_Color_Difference_VERBOSE

template< typename HalfedgeGraph, typename PointMap, typename VertexColorMap >
void
Compression_Valence_Component< HalfedgeGraph, PointMap, VertexColorMap >::
    Calculate_Edge_Color_Difference(HalfedgeGraph &_pMesh,
                                    const int &_Component_ID,
                                    double &_Max_color,
                                    double &_Mean_color,
                                    int &Number_of_vertices)
{
  typedef boost::graph_traits< HalfedgeGraph > GraphTraits;
  typedef typename GraphTraits::vertex_iterator vertex_iterator;
  typedef typename GraphTraits::halfedge_iterator halfedge_iterator;

#ifdef DBG_Calculate_Edge_Color_Difference
  static unsigned int Calculate_Edge_Color_Difference_call_cnt = 0;
  Calculate_Edge_Color_Difference_call_cnt++;
  std::string dbg_header =
      std::string("in ") + std::string(__func__) + std::string(", call_cnt=") +
      std::to_string(Calculate_Edge_Color_Difference_call_cnt);
  std::cout << dbg_header << ", this->C0_Min = " << this->C0_Min << std::endl;
  std::cout << dbg_header << ", this->C1_Min = " << this->C1_Min << std::endl;
  std::cout << dbg_header << ", this->C2_Min = " << this->C2_Min << std::endl;
  vertex_iterator vi_beg_dbg = vertices(_pMesh).first;
  vertex_iterator vi_end_dbg = vertices(_pMesh).second;
  for(vertex_iterator vi = vi_beg_dbg; vi != vi_end_dbg; ++vi)
  {
    std::cout << dbg_header
              << ", vertex_Seed_Edge = " << get(this->vertex_Seed_Edge, *vi)
              << std::endl;
    Color_Unit color_int;
    color_int = get(this->vertex_color_int, *vi);
    std::cout << dbg_header << ", vertex_color_int = " << color_int.c0 << ","
              << color_int.c1 << "," << color_int.c2 << std::endl;
  }
  PointMap pm = get(boost::vertex_point, _pMesh);
  PointMap *_pm = &pm;
#endif

  float C0_min = this->C0_Min;
  float C1_min = this->C1_Min;
  float C2_min = this->C2_Min;

  float Color_small_step = 0.0;

  if(this->NumberColorQuantization[_Component_ID] == 0)
  {
    Color_small_step = this->Color_Quantization_Step;
  }

  else
  {
    Color_small_step =
        this->Color_Quantization_Step *
        std::pow(2.0, this->NumberColorQuantization[_Component_ID]);
  }
#ifdef DBG_Calculate_Edge_Color_Difference
  std::cout << dbg_header << ", Color_small_step = " << Color_small_step
            << std::endl;
#endif

  // To find first points to start the conquest.
  auto halfedge_iterator_pair = halfedges(_pMesh);
  halfedge_iterator hi = halfedge_iterator_pair.first;

  while(
      (get(this->vertex_Seed_Edge, target(*hi, _pMesh)) != 2 * _Component_ID) ||
      (get(this->vertex_Seed_Edge, target(opposite(*hi, _pMesh), _pMesh)) !=
       2 * _Component_ID + 1))
  {
    hi++;
  }
#ifdef DBG_Calculate_Edge_Color_Difference
  print_halfedge(
      "in Calculate_Edge_Color_Difference, *hi = ", *hi, _pMesh, _pm);
#endif

  // Vertex_Flag est donnee free a tous les sommets
  auto vertex_iterator_pair = vertices(_pMesh);
  vertex_iterator vi_beg = vertex_iterator_pair.first;
  vertex_iterator vi_end = vertex_iterator_pair.second;
  for(vertex_iterator pVert = vi_beg; pVert != vi_end; pVert++)
  {
    put(this->Vertex_Flag, *pVert, FREE);
    put(this->Vertex_Number, *pVert, -1);
  }

  std::queue< vertex_descriptor > vertices;

  // push a input gate to begin loop
  // in order to treat all vertices.
  vertices.push(target(*hi, _pMesh));
  vertices.push(target(opposite(*hi, _pMesh), _pMesh));

  // To count number of vertices;;
  int Count_vertices = 0;
  double Max_color = -5000.;
  double Mean_color = 0.;
  double L_min = 5000, A_min = 5000, B_min = 5000;
  double L_max = -5000, A_max = -5000, B_max = -5000;

  unsigned int while_loop_cnt = 0;
  while(!vertices.empty())
  {
#ifdef DBG_Calculate_Edge_Color_Difference_VERBOSE
    std::cout << dbg_header << ", while_loop_cnt = " << ++while_loop_cnt
              << std::endl;
    std::cout << dbg_header << ", vertices.size = " << vertices.size()
              << std::endl;
#endif
    vertex_descriptor v = vertices.front();
    vertices.pop();
#ifdef DBG_Calculate_Edge_Color_Difference_VERBOSE
    print_vertex("in Calculate_Edge_Color_Difference, pop v = ", v, _pm);
#endif


    if(get(this->Vertex_Flag, v) == CONQUERED)
    {
#ifdef DBG_Calculate_Edge_Color_Difference_VERBOSE
      // std::cout << "in " << __func__ << "_mark_#1" << std::endl;
#endif
      continue;
    }
    else
    {
#ifdef DBG_Calculate_Edge_Color_Difference_VERBOSE
      // std::cout << "in " << __func__ << "_mark_#2" << std::endl;
      std::vector< std::string > dbg_diff_per_he;
#endif
      put(this->Vertex_Flag, v, CONQUERED);
      put(this->Vertex_Number, v, Count_vertices);

      CGAL::Halfedge_around_target_circulator< HalfedgeGraph > hvc(v, _pMesh);
      CGAL::Halfedge_around_target_circulator< HalfedgeGraph > phvc(hvc);

      int neighbour_cnt = 0;
      double Mean = 0.0;
      int Count = 0;
      CGAL_For_all(hvc, phvc)
      {
#ifdef DBG_Calculate_Edge_Color_Difference_VERBOSE
        // print_halfedge("in Calculate_Edge_Color_Difference, *hvc = ", *hvc,
        // _pMesh, _pm); if( vertex_to_string(v, _pm) == "(0.000483196,
        // 0.0101471, 0.714355)")
        //{
        // std::cout << "in " << __func__ << ", v=" << vertex_to_string(v, _pm)
        // << "  neighbour=" << ++neighbour_cnt << "  v0=" <<
        // vertex_to_string(target(*hvc, _pMesh), _pm) << "  v1=" <<
        // vertex_to_string(target(opposite(*hvc, _pMesh), _pMesh), _pm) << "
        // Vertex_Flag=" << (get(this->Vertex_Flag, target(opposite(*hvc,
        // _pMesh), _pMesh)) == FREE) << std::endl;
        //}
#endif

        if(get(this->Vertex_Flag, target(opposite(*hvc, _pMesh), _pMesh)) ==
           FREE)
        {
#ifdef DBG_Calculate_Edge_Color_Difference_VERBOSE
          // std::cout << "in " << __func__ << "_mark_#3" << std::endl;
#endif
          Color_Unit Color_0, Color_1;
          Color_0 = get(this->vertex_color_int, target(*hvc, _pMesh));

          Color_1 = get(this->vertex_color_int,
                        target(opposite(*hvc, _pMesh), _pMesh));

          double LAB_0[3], LAB_1[3];

          LAB_0[0] = C0_min + Color_0.c0 * Color_small_step;
          LAB_0[1] = C1_min + Color_0.c1 * Color_small_step;
          LAB_0[2] = C2_min + Color_0.c2 * Color_small_step;

          LAB_1[0] = C0_min + Color_1.c0 * Color_small_step;
          LAB_1[1] = C1_min + Color_1.c1 * Color_small_step;
          LAB_1[2] = C2_min + Color_1.c2 * Color_small_step;

          if(LAB_0[0] > L_max)
            L_max = LAB_0[0];
          if(LAB_0[1] > A_max)
            A_max = LAB_0[1];
          if(LAB_0[2] > B_max)
            B_max = LAB_0[2];

          if(LAB_1[0] > L_max)
            L_max = LAB_1[0];
          if(LAB_1[1] > A_max)
            A_max = LAB_1[1];
          if(LAB_1[2] > B_max)
            B_max = LAB_1[2];

          if(LAB_0[0] < L_min)
            L_min = LAB_0[0];
          if(LAB_0[1] < A_min)
            A_min = LAB_0[1];
          if(LAB_0[2] < B_min)
            B_min = LAB_0[2];

          if(LAB_1[0] < L_min)
            L_min = LAB_1[0];
          if(LAB_1[1] < A_min)
            A_min = LAB_1[1];
          if(LAB_1[2] < B_min)
            B_min = LAB_1[2];
          double diff = 0.0;
          for(int i = 0; i < 3; i++)
          {
            diff += (LAB_0[i] - LAB_1[i]) * (LAB_0[i] - LAB_1[i]);
          }

          diff = std::sqrt(diff);
#ifdef DBG_Calculate_Edge_Color_Difference_VERBOSE
          // std::cout << "in " << __func__ << ", v0 = " <<
          // vertex_to_string(target(*hvc, _pMesh), _pm) << "  v1 = " <<
          // vertex_to_string(target(opposite(*hvc, _pMesh), _pMesh), _pm) << "
          // diff = " << diff << std::endl;
          char dbg_buffer[20];
          snprintf(dbg_buffer, sizeof(dbg_buffer), "%f", diff);
          dbg_diff_per_he.push_back(edge_to_string(*hvc, _pMesh, _pm) +
                                    "  diff=" + dbg_buffer);
#endif

          Mean += diff;
          Count++;

#ifdef DBG_Calculate_Edge_Color_Difference_VERBOSE
          // if( vertex_to_string(v, _pm) == "(0.000483196, 0.0101471,
          // 0.714355)")
          //{
          // std::cout << "in " << __func__ << ", v=" << vertex_to_string(v,
          // _pm) << "  v0=" << vertex_to_string(target(*hvc, _pMesh), _pm) << "
          // v1=" << vertex_to_string(target(opposite(*hvc, _pMesh), _pMesh),
          // _pm) << "  diff=" << diff << "  Mean=" << Mean << "  Count=" <<
          // Count << std::endl;
          //}
#endif
        }
      }
#ifdef DBG_Calculate_Edge_Color_Difference_VERBOSE
      // std::cout << "in " << __func__ << ", Mean = " << Mean << std::endl;
      // std::cout << "in " << __func__ << ", Count = " << Count << std::endl;
      for(int i = 0; i < dbg_diff_per_he.size(); i++)
        std::cout << dbg_header << "  " << dbg_diff_per_he[i]
                  << " Count=" << Count << std::endl;
#endif

      if(Count != 0)
      {
        Mean /= Count;
        Mean_color += Mean;
      }
#ifdef DBG_Calculate_Edge_Color_Difference_VERBOSE
      std::cout << dbg_header << ", v = " << vertex_to_string(v, _pm)
                << "  Mean = " << Mean << std::endl;
      std::cout << dbg_header << ", v = " << vertex_to_string(v, _pm)
                << "  Mean_color = " << Mean_color << std::endl;
#endif

      CGAL::Halfedge_around_target_circulator< HalfedgeGraph > h(v, _pMesh);
      CGAL::Halfedge_around_target_circulator< HalfedgeGraph > h2(h);
      CGAL_For_all(h, h2)
      {
#ifdef DBG_Calculate_Edge_Color_Difference_VERBOSE
        // print_halfedge("in Calculate_Edge_Color_Difference, *h = ", *h,
        // _pMesh, _pm);
#endif
        if(get(this->Vertex_Flag, target(opposite(*h, _pMesh), _pMesh)) == FREE)
        {
#ifdef DBG_Calculate_Edge_Color_Difference
          std::cout << dbg_header << ", v = " << vertex_to_string(v, _pm)
                    << "  push vertex "
                    << vertex_to_string(target(opposite(*h, _pMesh), _pMesh),
                                        _pm)
                    << std::endl;
#endif
          vertices.push(target(opposite(*h, _pMesh), _pMesh));
        }
      }

      // increment number of vertices;
      Count_vertices++;
    }
  }
#ifdef DBG_Calculate_Edge_Color_Difference
  std::cout << dbg_header << ", after while loop Mean_color = " << Mean_color
            << std::endl;
#endif

  Max_color = (L_max - L_min) * (L_max - L_min) +
              (A_max - A_min) * (A_max - A_min) +
              (B_max - B_min) * (B_max - B_min);

  Max_color = std::sqrt(Max_color);
  _Max_color = Max_color;
  _Mean_color = 3 * Mean_color / Count_vertices;
  Number_of_vertices = Count_vertices;

#ifdef DBG_Calculate_Edge_Color_Difference
  for(vertex_iterator vi = vi_beg_dbg; vi != vi_end_dbg; ++vi)
  {
    std::cout << dbg_header << ", Vertex_Flag = " << get(this->Vertex_Flag, *vi)
              << std::endl;
    std::cout << dbg_header
              << ", Vertex_Number = " << get(this->Vertex_Number, *vi)
              << std::endl;
    // note: Vertex_Number is the order in which the vertices are processed, so
    // it's different in Mepp1 and Mepp2
  }
  std::cout << dbg_header << ", _Mean_color = " << _Mean_color << std::endl;
  std::cout << dbg_header << ", _Max_color = " << _Max_color << std::endl;
  std::cout << dbg_header << ", Number_of_vertices = " << Number_of_vertices
            << std::endl;
#endif
}

// Differentes histogrammes pour chaque couleur.
template< typename HalfedgeGraph, typename PointMap, typename VertexColorMap >
void
Compression_Valence_Component< HalfedgeGraph, PointMap, VertexColorMap >::
    Diminush_Color_Quantization_Precision(HalfedgeGraph &_pMesh,
                                          const int _Component_ID,
                                          VertexColorMap *_v_cm)
{
  typedef typename boost::property_traits< VertexColorMap >::value_type Color;

  // int N = 16; // Number of neighboring colors to use;
  // float Color_diff_seuil = 30.0;

  /*float C0_min = this->C0_Min;
  float C1_min = this->C1_Min;
  float C2_min = this->C2_Min;*/

  float Color_small_step = 0.0;

  if(this->NumberColorQuantization[_Component_ID] == 0)
    Color_small_step = this->Color_Quantization_Step;

  else
    Color_small_step =
        this->Color_Quantization_Step *
        std::pow(2.0, this->NumberColorQuantization[_Component_ID]);

  double Color_large_step = Color_small_step * 2;

  // To find first points to start the conquest.
  auto halfedge_iterator_pair = halfedges(_pMesh);
  halfedge_iterator hi = halfedge_iterator_pair.first;

  while(
      (get(this->vertex_Seed_Edge, target(*hi, _pMesh)) != 2 * _Component_ID) ||
      (get(this->vertex_Seed_Edge, target(opposite(*hi, _pMesh), _pMesh)) !=
       2 * _Component_ID + 1))
    hi++;

  // Vertex_Flag est donnee free a tous les sommets
  for(vertex_iterator pVert = vertices(_pMesh).first;
      pVert != vertices(_pMesh).second;
      pVert++)
  {
    put(this->Vertex_Flag, *pVert, FREE);
    put(this->Vertex_Number, *pVert, -1);
  }

  // premiere passe pour sous quantifier et donne l'indice de symbol a chaque
  // sommet
  for(vertex_iterator pVert = vertices(_pMesh).first;
      pVert != vertices(_pMesh).second;
      pVert++)
  {
    if(get(vertex_Component_Number, *pVert) == _Component_ID)
    {
      Color_Unit color = get(this->vertex_color_int, *pVert);
      int LAB_init_C0 = color.c0;
      int LAB_init_C1 = color.c1;
      int LAB_init_C2 = color.c2;

      // Les coordonnees apres under-quantization
      int LAB_final_C0 = LAB_init_C0 / 2;
      int LAB_final_C1 = LAB_init_C1 / 2;
      int LAB_final_C2 = LAB_init_C2 / 2;

      int i = LAB_init_C0 % 2;
      int j = LAB_init_C1 % 2;
      int k = LAB_init_C2 % 2;

      int Q_index = Get_Correct_Vector(i, j, k);

      put(this->vertex_color_int,
          *pVert,
          Color_Unit(LAB_final_C0, LAB_final_C1, LAB_final_C2));
      put(this->vertex_Q_Index, *pVert, Q_index);

      float RGB[3];

      // To re-colorify the vertex with decreased color precision
      std::vector< float > LAB;

      LAB.push_back(this->C0_Min + LAB_final_C0 * Color_large_step);
      LAB.push_back(this->C1_Min + LAB_final_C1 * Color_large_step);
      LAB.push_back(this->C2_Min + LAB_final_C2 * Color_large_step);

      LAB_To_RGB(LAB[0], LAB[1], LAB[2], RGB);

      put(*_v_cm, *pVert, Color(RGB[0], RGB[1], RGB[2]));
    }
  }

  std::queue< vertex_descriptor > vertices;

  std::list< int > Inter_color_quantization;
  std::list< int > Inter_index_per_color;

  // push a input gate to begin loop
  // in order to treat all vertices.
  vertices.push(target(*hi, _pMesh));
  vertices.push(target(opposite(*hi, _pMesh), _pMesh));

  int Vertex_index = 0;

  // Color table which stocks all present colors
  std::vector< std::vector< int > > Color_table;

  while(!vertices.empty())
  {
    vertex_descriptor v = vertices.front();
    vertices.pop();

    if(get(this->Vertex_Flag, v) == CONQUERED)
      continue;

    else
    {
      put(this->Vertex_Flag, v, CONQUERED);
      put(this->Vertex_Number, v, Vertex_index);

      int Q_index = get(this->vertex_Q_Index, v);
      std::vector< int > LAB;

      LAB.push_back((get(this->vertex_color_int, v)).c0);
      LAB.push_back((get(this->vertex_color_int, v)).c1);
      LAB.push_back((get(this->vertex_color_int, v)).c2);

      bool Is_existed_color = false;

      // C'est une couleur deja existante?
      for(unsigned i = 0; i < Color_table.size(); i++)
      {
        if((LAB[0] == Color_table[i][0]) && (LAB[1] == Color_table[i][1]) &&
           (LAB[2] == Color_table[i][2]))
        {
          Is_existed_color = true;
          Inter_color_quantization.push_front(Q_index);
          Inter_index_per_color.push_front(i);
          break;
        }
      }

      // Sinon, agrandir le tabeau.
      if(Is_existed_color == false)
      {
        Inter_color_quantization.push_front(Q_index);
        Inter_index_per_color.push_front(Color_table.size());

        Color_table.push_back(LAB);
      }
      CGAL::Halfedge_around_target_circulator< HalfedgeGraph > h;

      // First_vertex -> To find the seed edge
      if(Vertex_index == 0)
      {
        h = CGAL::Halfedge_around_target_circulator< HalfedgeGraph >(v, _pMesh);
        do
        {
          h++;
        } while((*h) != (*hi));
      }

      // To find an deterministic order = a given vertex and a vertex
      // with the highest value of vertex number
      else
      {
        int Comp_number = -2;

        h = CGAL::Halfedge_around_target_circulator< HalfedgeGraph >(v, _pMesh);

        CGAL::Halfedge_around_target_circulator< HalfedgeGraph > h2 = h;
        CGAL_For_all(h, h2)
        {
          if(get(this->Vertex_Number, target(opposite(*h, _pMesh), _pMesh)) >
             Comp_number)
            Comp_number =
                get(this->Vertex_Number, target(opposite(*h, _pMesh), _pMesh));
        }

        h = h2;
        CGAL_For_all(h, h2)
        {
          if(get(this->Vertex_Number, target(opposite(*h, _pMesh), _pMesh)) ==
             Comp_number)
            break;
        }
      }

      CGAL::Halfedge_around_target_circulator< HalfedgeGraph > h2 = h;
      CGAL_For_all(h, h2)
      {
        if(get(this->Vertex_Flag, target(opposite(*h, _pMesh), _pMesh)) == FREE)
          vertices.push(target(opposite(*h, _pMesh), _pMesh));
      }

      Vertex_index++;
    }
  }


  int ssss = Inter_color_quantization.size();
  this->NumberProcessedVertices[_Component_ID].push_front(ssss);

  while(!Inter_color_quantization.empty())
  {
    int index = Inter_color_quantization.front();
    Inter_color_quantization.pop_front();

    this->ColorChildcellIndex[_Component_ID].push_front(index);
  }
  while(!Inter_index_per_color.empty())
  {
    int Color = Inter_index_per_color.front();
    Inter_index_per_color.pop_front();

    this->ColorEncoderIndex[_Component_ID].push_front(Color);
  }
}


template< typename HalfedgeGraph, typename PointMap, typename VertexColorMap >
void
Compression_Valence_Component< HalfedgeGraph, PointMap, VertexColorMap >::
    Recalculate_Component_Area(HalfedgeGraph &_pMesh,
                               const PointMap *_pm,
                               const int &_Component_ID,
                               int &Number_facets)
{
  typedef boost::graph_traits< HalfedgeGraph > GraphTraits;
  typedef typename GraphTraits::face_iterator face_iterator;

  double Area = 0.0;
  Number_facets = 0;
  auto face_iterator_pair = faces(_pMesh);
  face_iterator pFacet_beg = face_iterator_pair.first;
  face_iterator pFacet_end = face_iterator_pair.second;
  for(face_iterator pFacet = pFacet_beg; pFacet != pFacet_end; ++pFacet)
  {
    if(get(this->facet_Component_Number, *pFacet) == _Component_ID)
    {
      Number_facets += 1;
      Area += Area_Facet_Triangle(halfedge(*pFacet, _pMesh), _pMesh, _pm);
    }
  }
  Area *= std::pow(
      ((double)10.0 / (double)this->HighestLengthBB[_Component_ID]), 2.0);
  this->ComponentArea[_Component_ID] = Area;
}


template< typename HalfedgeGraph, typename PointMap, typename VertexColorMap >
void
Compression_Valence_Component< HalfedgeGraph, PointMap, VertexColorMap >::
    Augment_Geometry_Quantization_Precision(HalfedgeGraph &_pMesh,
                                            Arithmetic_Codec &Decoder,
                                            const int &Component_ID,
                                            PointMap *_pm)
{
  FEVV::Geometry_traits< HalfedgeGraph > gt(_pMesh);

  Adaptive_Data_Model Under_quantization_model(8);

  // Premier_Pas == distance d'un Quantization_Step de grill de quantification
  // (Q)
  double Small_step = 0.0;

  if(this->NumberChangeQuantization[Component_ID] == 0)
    Small_step = this->Quantization_Step[Component_ID];

  else
    Small_step =
        this->Quantization_Step[Component_ID] *
        std::pow(2.0, this->NumberChangeQuantization[Component_ID] - 1);

  // double Large_step = Small_step * 2;

  // To find first points to start the conquest.
  halfedge_iterator hi = halfedges(_pMesh).first;
  while(
      (get(this->vertex_Seed_Edge, target(*hi, _pMesh)) != 2 * Component_ID) ||
      (get(this->vertex_Seed_Edge, target(opposite(*hi, _pMesh), _pMesh)) !=
       2 * Component_ID + 1))
    hi++;

  // Vertex_Flag est donnee free a tous les sommets
  for(vertex_iterator pVert = vertices(_pMesh).first;
      pVert != vertices(_pMesh).second;
      pVert++)
  {
    put(this->Vertex_Flag, *pVert, FREE);
    put(this->Vertex_Number, *pVert, -1);
    put(this->vertex_Component_Number, *pVert, -1);
  }
  //_pMesh.compute_normals();

  std::queue< vertex_descriptor > verticesq;
  std::list< int > Layer_symbols;

  // push first vertex to begin loop
  // to treat all vertices.
  verticesq.push(target(*hi, _pMesh));
  verticesq.push(target(opposite(*hi, _pMesh), _pMesh));

  int Vertex_index = 0;

  while(!verticesq.empty())
  {
    vertex_descriptor v = verticesq.front();
    verticesq.pop();

    if(get(this->Vertex_Flag, v) == CONQUERED)
      continue;

    else
    {
      put(this->Vertex_Flag, v, CONQUERED);
      put(this->Vertex_Number, v, Vertex_index);
      put(this->vertex_Component_Number, v, Component_ID);

      unsigned valence = out_degree(v, _pMesh);

      Point3d *Neighbors = new Point3d[valence];
      CGAL::Halfedge_around_target_circulator< HalfedgeGraph > hvc(v, _pMesh);
      CGAL::Halfedge_around_target_circulator< HalfedgeGraph > hvc2 = hvc;
      unsigned count_neighbor = 0;
      CGAL_For_all(hvc, hvc2)
      {
        Neighbors[count_neighbor] =
            get(*_pm, target(opposite(*hvc, _pMesh), _pMesh));
        count_neighbor++;
      }

      Point3d Center = get(*_pm, v);
      std::vector< int > D(6, 0);

      for(unsigned i = 0; i < valence; i++)
      {
        Vector Vect =
            FEVV::Math::Vector::sub< FEVV::Geometry_traits< HalfedgeGraph > >(
                Neighbors[i], Center, gt);

        if(Vect[0] <= 0)
        {
          int I_dist = Vect[0] / Small_step / 2;
          D[0] += std::abs(I_dist);
        }
        else
        {
          int I_dist = Vect[0] / Small_step / 2;
          D[1] += std::abs(I_dist);
        }

        if(Vect[1] <= 0)
        {
          int I_dist = Vect[1] / Small_step / 2;
          D[2] += std::abs(I_dist);
        }
        else
        {
          int I_dist = Vect[1] / Small_step / 2;
          D[3] += std::abs(I_dist);
        }

        if(Vect[2] <= 0)
        {
          int I_dist = Vect[2] / Small_step / 2;
          D[4] += std::abs(I_dist);
        }
        else
        {
          int I_dist = Vect[2] / Small_step / 2;
          D[5] += std::abs(I_dist);
        }
      }

      std::vector< double > U(3, 0);
      U[0] = std::abs((double)D[0] / (double)(D[0] + D[1]) - 0.5);
      U[1] = std::abs((double)D[2] / (double)(D[2] + D[3]) - 0.5);
      U[2] = std::abs((double)D[4] / (double)(D[4] + D[5]) - 0.5);

      std::multimap< double, int > U_order;
      U_order.insert(std::pair< double, int >(U[0], 0));
      U_order.insert(std::pair< double, int >(U[1], 1));
      U_order.insert(std::pair< double, int >(U[2], 2));

      std::multimap< double, int >::iterator it;

      std::vector< int > Weight_axe(3, 0);
      int Temp_weight = 1;
      for(it = U_order.begin(); it != U_order.end(); it++)
      {
        Weight_axe[it->second] = Temp_weight;
        Temp_weight++;
      }

      std::vector< int > Priority(8, 0);
      Priority[0] =
          Weight_axe[0] * D[0] + Weight_axe[1] * D[2] + Weight_axe[2] * D[4];
      Priority[1] =
          Weight_axe[0] * D[1] + Weight_axe[1] * D[2] + Weight_axe[2] * D[4];
      Priority[2] =
          Weight_axe[0] * D[0] + Weight_axe[1] * D[3] + Weight_axe[2] * D[4];
      Priority[3] =
          Weight_axe[0] * D[1] + Weight_axe[1] * D[3] + Weight_axe[2] * D[4];
      Priority[4] =
          Weight_axe[0] * D[0] + Weight_axe[1] * D[2] + Weight_axe[2] * D[5];
      Priority[5] =
          Weight_axe[0] * D[1] + Weight_axe[1] * D[2] + Weight_axe[2] * D[5];
      Priority[6] =
          Weight_axe[0] * D[0] + Weight_axe[1] * D[3] + Weight_axe[2] * D[5];
      Priority[7] =
          Weight_axe[0] * D[1] + Weight_axe[1] * D[3] + Weight_axe[2] * D[5];

      std::multimap< int, int > Priority_map;
      std::multimap< int, int >::iterator it_reorder;
      for(int i = 0; i < 8; i++)
        Priority_map.insert(std::pair< int, int >(Priority[i], i));

      std::vector< int > Priority_reorder(8, 0);

      int Temp_priority = 0;
      ///////////////////
      bool Is_same_priority_value = false;

      it_reorder = Priority_map.begin();
      for(int i = 0; i < 7; i++)
      {
        int P0 = it_reorder->first;
        it_reorder++;
        int P1 = it_reorder->first;

        if(P0 == P1)
          Is_same_priority_value = true;
      }

      if(!Is_same_priority_value)
      {
        for(it_reorder = Priority_map.begin(); it_reorder != Priority_map.end();
            it_reorder++)
        {
          Priority_reorder[it_reorder->second] = 7 - Temp_priority;
          Temp_priority++;
        }
      }
      else
      {
        for(int i = 0; i < 8; i++)
        {
          Priority_reorder[i] = i;
        }
      }

      int Reordered_number = Decoder.decode(Under_quantization_model);

      for(int i = 0; i < 8; i++)
      {
        if(Priority_reorder[i] == Reordered_number)
        {
          put(this->vertex_Q_Index, v, i);
          break;
        }
      }

      CGAL::Halfedge_around_target_circulator< HalfedgeGraph > h;
      if(Vertex_index == 0)
      {
        h = CGAL::Halfedge_around_target_circulator< HalfedgeGraph >(v, _pMesh);
        do
        {
          h++;
        } while((*h) != (*hi));
      }
      else
      {
        int Comp_number = -2;
        h = CGAL::Halfedge_around_target_circulator< HalfedgeGraph >(v, _pMesh);
        CGAL::Halfedge_around_target_circulator< HalfedgeGraph > h2 = h;
        CGAL_For_all(h, h2)
        {
          if(get(this->Vertex_Number, target(opposite(*h, _pMesh), _pMesh)) >
             Comp_number)
          {
            Comp_number =
                get(this->Vertex_Number, target(opposite(*h, _pMesh), _pMesh));
          }
        }

        h = h2;
        CGAL_For_all(h, h2)
        {
          if(get(this->Vertex_Number, target(opposite(*h, _pMesh), _pMesh)) ==
             Comp_number)
            break;
        }
      }

      CGAL::Halfedge_around_target_circulator< HalfedgeGraph > h2 = h;
      CGAL_For_all(h, h2)
      {
        if(get(this->Vertex_Flag, target(opposite(*h, _pMesh), _pMesh)) == FREE)
        {
          verticesq.push(target(opposite(*h, _pMesh), _pMesh));
        }
      }
      Vertex_index++;

      delete[] Neighbors;
    }
  }
  vertex_iterator pVertex;
  for(pVertex = vertices(_pMesh).first; pVertex != vertices(_pMesh).second;
      pVertex++)
  {
    if(get(this->vertex_Component_Number, *pVertex) == Component_ID)
    {
      int Coeff[3];
      Get_Coefficient_Up_Quantization(get(this->vertex_Q_Index, *pVertex),
                                      Coeff);
      Point3d pt = get(*_pm, *pVertex);
      Point3d New_position(pt[0] + Coeff[0] * Small_step / 2,
                           pt[1] + Coeff[1] * Small_step / 2,
                           pt[2] + Coeff[2] * Small_step / 2);
      put(*_pm, *pVertex, New_position);
    }
  }
  this->Qbit[Component_ID]++;
  this->NumberChangeQuantization[Component_ID]--;
}


/*  Description : ADAPTIVE_QUANTIZATION
Decreasing of quantization resolution based on the prediction of PENG.
Opposite function is up_quantization. */
template< typename HalfedgeGraph, typename PointMap, typename VertexColorMap >
void
Compression_Valence_Component< HalfedgeGraph, PointMap, VertexColorMap >::
    Diminush_Geometry_Quantization_Precision(HalfedgeGraph &_pMesh,
                                             const int &Component_ID,
                                             PointMap *_pm)
{
  FEVV::Geometry_traits< HalfedgeGraph > gt(_pMesh);

  // stock three mins for the reconstruction.
  float _xmin = this->xmin[Component_ID];
  float _ymin = this->ymin[Component_ID];
  float _zmin = this->zmin[Component_ID];

  // Premier_Pas == distance d'un Quantization_Step de grill de quantification
  // (Q)
  double Small_step = 0.0;

  if(this->NumberChangeQuantization[Component_ID] == 0)
    Small_step = this->Quantization_Step[Component_ID];

  else
    Small_step = this->Quantization_Step[Component_ID] *
                 std::pow(2.0, this->NumberChangeQuantization[Component_ID]);

  // Large_step == distance d'un Quantization_Step de grille de quantification(Q
  // - 1)
  double Large_step = Small_step * 2;

  // To find first points to start the conquest.
  auto halfedge_iterator_pair = halfedges(_pMesh);
  halfedge_iterator hi = halfedge_iterator_pair.first;

  while(
      (get(this->vertex_Seed_Edge, target(*hi, _pMesh)) != 2 * Component_ID) ||
      (get(this->vertex_Seed_Edge, target(opposite(*hi, _pMesh), _pMesh)) !=
       2 * Component_ID + 1))
    hi++;

  // Vertex_Flag est donnee free a tous les sommets
  for(vertex_iterator pVert = vertices(_pMesh).first;
      pVert != vertices(_pMesh).second;
      pVert++)
  {
    put(this->Vertex_Flag, *pVert, FREE);
    put(this->Vertex_Number, *pVert, -1);
  }


  // premiere passe pour sous quantifie et donne l'indice de symbol a chaque
  // sommet
  for(vertex_iterator pVert = vertices(_pMesh).first;
      pVert != vertices(_pMesh).second;
      pVert++)
  {
    if(get(this->vertex_Component_Number, *pVert) == Component_ID)
    {
      Point3d pt = get(*_pm, *pVert);
      double x = pt[0];
      double y = pt[1];
      double z = pt[2];

      // Qx, Qy, Qz sont des coordonnees avant under-Quantization
      int Qx = (int)(ceil((x - _xmin) / Small_step)) - 1;
      if(Qx == -1)
        Qx = 0;
      int Qy = (int)(ceil((y - _ymin) / Small_step)) - 1;
      if(Qy == -1)
        Qy = 0;
      int Qz = (int)(ceil((z - _zmin) / Small_step)) - 1;
      if(Qz == -1)
        Qz = 0;

      // Les coordonnees apres under-quantization
      int Second_Qx = Qx / 2;
      int Second_Qy = Qy / 2;
      int Second_Qz = Qz / 2;

      int i = Qx % 2;
      int j = Qy % 2;
      int k = Qz % 2;

      int Q_index = Get_Correct_Vector(i, j, k);
      put(this->vertex_Q_Index, *pVert, Q_index);

      // Les sommets sont deplaces vers le centre de cellule mere.
      put(*_pm,
          *pVert,
          Point3d(_xmin + (Second_Qx + 0.5) * Large_step,
                  _ymin + (Second_Qy + 0.5) * Large_step,
                  _zmin + (Second_Qz + 0.5) * Large_step));
    }
  }

  /////// appliquer calcul des normales?? ou Quantization_Step???
  //_pMesh.compute_normals();

  std::queue< vertex_descriptor > vertices;

  std::list< int > Layer_symbols;

  // push a input gate to begin loop
  // in order to treat all vertices.

  vertices.push(target(*hi, _pMesh));
  vertices.push(target(opposite(*hi, _pMesh), _pMesh));

  int Vertex_index = 0;


  while(!vertices.empty())
  {
    vertex_descriptor v = vertices.front();
    vertices.pop();

    if(get(this->Vertex_Flag, v) == CONQUERED)
      continue;

    else
    {
      put(this->Vertex_Flag, v, CONQUERED);
      put(this->Vertex_Number, v, Vertex_index);

      int Q_index = get(this->vertex_Q_Index, v);
      unsigned Count_neighbor = 0;
      unsigned valence = out_degree(v, _pMesh);

      Point3d *Neighbors = new Point3d[valence];
      CGAL::Halfedge_around_target_circulator< HalfedgeGraph > hvc =
          CGAL::Halfedge_around_target_circulator< HalfedgeGraph >(v, _pMesh);
      CGAL::Halfedge_around_target_circulator< HalfedgeGraph > phvc = hvc;

      CGAL_For_all(hvc, phvc)
      {
        Neighbors[Count_neighbor] =
            get(*_pm, target(opposite(*hvc, _pMesh), _pMesh));
        Count_neighbor++;
      }

      Point3d Center = get(*_pm, v);

      std::vector< int > D(6, 0);

      for(unsigned i = 0; i < valence; i++)
      {
        Vector Vect =
            FEVV::Math::Vector::sub< FEVV::Geometry_traits< HalfedgeGraph > >(
                Neighbors[i], Center, gt);

        if(Vect[0] <= 0)
        {
          int I_dist = Vect[0] / Small_step / 2;
          D[0] += std::abs(I_dist);
        }
        else
        {
          int I_dist = Vect[0] / Small_step / 2;
          D[1] += std::abs(I_dist);
        }

        if(Vect[1] <= 0)
        {
          int I_dist = Vect[1] / Small_step / 2;
          D[2] += std::abs(I_dist);
        }
        else
        {
          int I_dist = Vect[1] / Small_step / 2;
          D[3] += std::abs(I_dist);
        }

        if(Vect[2] <= 0)
        {
          int I_dist = Vect[2] / Small_step / 2;
          D[4] += std::abs(I_dist);
        }
        else
        {
          int I_dist = Vect[2] / Small_step / 2;
          D[5] += std::abs(I_dist);
        }
      }

      std::vector< double > U(3, 0);
      U[0] = std::abs((double)D[0] / (double)(D[0] + D[1]) - 0.5);
      U[1] = std::abs((double)D[2] / (double)(D[2] + D[3]) - 0.5);
      U[2] = std::abs((double)D[4] / (double)(D[4] + D[5]) - 0.5);

      std::multimap< double, int > U_order;
      U_order.insert(std::pair< double, int >(U[0], 0));
      U_order.insert(std::pair< double, int >(U[1], 1));
      U_order.insert(std::pair< double, int >(U[2], 2));

      std::multimap< double, int >::iterator it;

      std::vector< int > Weight_axe(3, 0);
      int Temp_weight = 1;
      for(it = U_order.begin(); it != U_order.end(); it++)
      {
        Weight_axe[it->second] = Temp_weight;
        Temp_weight++;
      }

      std::vector< int > Priority(8, 0);
      Priority[0] =
          Weight_axe[0] * D[0] + Weight_axe[1] * D[2] + Weight_axe[2] * D[4];
      Priority[1] =
          Weight_axe[0] * D[1] + Weight_axe[1] * D[2] + Weight_axe[2] * D[4];
      Priority[2] =
          Weight_axe[0] * D[0] + Weight_axe[1] * D[3] + Weight_axe[2] * D[4];
      Priority[3] =
          Weight_axe[0] * D[1] + Weight_axe[1] * D[3] + Weight_axe[2] * D[4];
      Priority[4] =
          Weight_axe[0] * D[0] + Weight_axe[1] * D[2] + Weight_axe[2] * D[5];
      Priority[5] =
          Weight_axe[0] * D[1] + Weight_axe[1] * D[2] + Weight_axe[2] * D[5];
      Priority[6] =
          Weight_axe[0] * D[0] + Weight_axe[1] * D[3] + Weight_axe[2] * D[5];
      Priority[7] =
          Weight_axe[0] * D[1] + Weight_axe[1] * D[3] + Weight_axe[2] * D[5];

      std::multimap< int, int > Priority_map;
      std::multimap< int, int >::iterator it_reorder;
      for(int i = 0; i < 8; i++)
        Priority_map.insert(std::pair< int, int >(Priority[i], i));

      std::vector< int > Priority_reorder(8, 0);

      int Temp_priority = 0;
      bool Is_same_priority_value = false;

      it_reorder = Priority_map.begin();
      for(int i = 0; i < 7; i++)
      {
        int P0 = it_reorder->first;
        it_reorder++;
        int P1 = it_reorder->first;

        if(P0 == P1)
          Is_same_priority_value = true;
      }

      if(!Is_same_priority_value)
      {
        for(it_reorder = Priority_map.begin(); it_reorder != Priority_map.end();
            it_reorder++)
        {
          Priority_reorder[it_reorder->second] = 7 - Temp_priority;
          Temp_priority++;
        }
      }
      else
      {
        for(int i = 0; i < 8; i++)
        {
          Priority_reorder[i] = i;
        }
      }

      int Reordered_number = Priority_reorder[Q_index];

      CGAL::Halfedge_around_target_circulator< HalfedgeGraph > h;

      // First_vertex -> To find the seed edge
      if(Vertex_index == 0)
      {
        h = CGAL::Halfedge_around_target_circulator< HalfedgeGraph >(v, _pMesh);
        do
        {
          h++;
        } while((*h) != (*hi));
      }

      else
      {
        int Comp_number = -2;
        h = CGAL::Halfedge_around_target_circulator< HalfedgeGraph >(v, _pMesh);
        CGAL::Halfedge_around_target_circulator< HalfedgeGraph > h2 = h;
        CGAL_For_all(h, h2)
        {
          if(get(this->Vertex_Number, target(opposite(*h, _pMesh), _pMesh)) >
             Comp_number)
            Comp_number =
                get(this->Vertex_Number, target(opposite(*h, _pMesh), _pMesh));
        }

        h = h2;
        CGAL_For_all(h, h2)
        {
          if(get(this->Vertex_Number, target(opposite(*h, _pMesh), _pMesh)) ==
             Comp_number)
            break;
        }
      }

      CGAL::Halfedge_around_target_circulator< HalfedgeGraph > h2 = h;
      CGAL_For_all(h, h2)
      {
        if(get(this->Vertex_Flag, target(opposite(*h, _pMesh), _pMesh)) == FREE)
          vertices.push(target(opposite(*h, _pMesh), _pMesh));
      }
      Vertex_index++;

      delete[] Neighbors;
      Layer_symbols.push_front(Reordered_number);
    }
  }

  this->NumberQuantizationLayer[Component_ID].push_front(Layer_symbols.size());
  while(!Layer_symbols.empty())
  {
    int index = Layer_symbols.front();
    Layer_symbols.pop_front();

    this->QuantizationCorrectVector[Component_ID].push_front(index);
  }
}


template< typename HalfedgeGraph, typename PointMap, typename VertexColorMap >
void
Compression_Valence_Component< HalfedgeGraph, PointMap, VertexColorMap >::
    Decompress_Init(HalfedgeGraph &_pMesh,
                    PointMap *_pm,
                    VertexColorMap *_v_cm,
                    const std::string &Input_File_Name)
{
  this->File_name = Input_File_Name;

  if(FILE *file2 = fopen(this->File_name.c_str(), "r"))
  {
    fseek(file2, 0, SEEK_END);
    this->Compressed_file_size = ftell(file2);
    fclose(file2);
  }
  else
  {
    throw std::runtime_error("Error: failed to open file '" + this->File_name +
                             "'.");
  }

  CGAL::clear(_pMesh);

  this->IsClosed.clear();
  this->Qbit.clear();
  this->xmin.clear();
  this->ymin.clear();
  this->zmin.clear();
  this->Quantization_Step.clear();
  this->NumberColorQuantization.clear();
  this->NumberChangeQuantization.clear();
  this->ComponentOperations.clear();

  this->Decompress_count = 0;

  // ELO note:
  // the .p3d input file is opened here
  // and remains open at the end of the function
  // for further reading
  FILE *fp = fopen(this->File_name.c_str(), "rb");

  this->DM_JCW_MOVE_ERROR.set_alphabet(3);

  int res;

  res = fread(&this->Smallest_Alpha, sizeof(int), 1, fp);
  res = fread(&this->Smallest_Gamma, sizeof(int), 1, fp);
  res = fread(&this->Initial_file_size, sizeof(unsigned), 1, fp);
  res = fread(&this->NumberComponents, sizeof(int), 1, fp);

  float Qpas;
  float t_xmin, t_ymin, t_zmin;

  // Read geometry information for each component;
  for(int i = 0; i < this->NumberComponents; i++)
  {
    res = fread(&Qpas, sizeof(float), 1, fp);   // quantization step
    res = fread(&t_xmin, sizeof(float), 1, fp); // x_min
    res = fread(&t_ymin, sizeof(float), 1, fp); // y_min
    res = fread(&t_zmin, sizeof(float), 1, fp); // z_min

    this->Quantization_Step.push_back(Qpas);
    this->xmin.push_back(t_xmin);
    this->ymin.push_back(t_ymin);
    this->zmin.push_back(t_zmin);
  }

  char Colored, One_color;
  res = fread(&Colored, sizeof(char), 1, fp);

  this->IsOneColor = false;
  if(Colored == 1)
  {
    this->IsColored = true;
    res = fread(&One_color, sizeof(char), 1, fp);
    if(One_color == 1)
      this->IsOneColor = true;
  }
  else
    this->IsColored = false;

  // read color information for each component
  if((this->IsColored) && (!this->IsOneColor))
  {
    res = fread(&this->Color_Quantization_Step, sizeof(float), 1, fp);

    res = fread(&this->C0_Min,
                sizeof(float),
                1,
                fp); // smallest absolute position of c0
    res = fread(&this->C1_Min,
                sizeof(float),
                1,
                fp); // smallest absolute position of c1
    res = fread(&this->C2_Min,
                sizeof(float),
                1,
                fp); // smallest absolute position of c2

    res = fread(
        &this->Smallest_C0, sizeof(int), 1, fp); // smallest quantized positions
    res = fread(&this->Smallest_C1, sizeof(int), 1, fp);
    res = fread(&this->Smallest_C2, sizeof(int), 1, fp);
  }

  if((this->IsColored) && (this->IsOneColor))
  {
    res = fread(&this->OnlyColor[0],
                sizeof(float),
                1,
                fp); // smallest absolute position of c0
    res = fread(
        &this->OnlyColor[1], sizeof(float), 1, fp); // smallest value of c1
    res = fread(
        &this->OnlyColor[2], sizeof(float), 1, fp); // smallest value of c2
  }
  this->Decoder.set_buffer(AC_BUFFER);
  this->Decoder.read_from_file(fp);

  res += 0; // just to avoid gcc 4.6 warning

  // To know if each component is colored or not, and closed of not.
  for(int i = 0; i < this->NumberComponents; i++)
  {
    if(Decoder.get_bits(1) == 0)
      this->IsClosed.push_back(true);
    else
      this->IsClosed.push_back(false);
  }

  if(Decoder.get_bits(1) == 1)
    this->Is_Bijection_Enabled = true;
  else
    this->Is_Bijection_Enabled = false;

  this->GlobalCountOperation = -1;

  unsigned Max_Qbit = 0;
  for(int i = 0; i < this->NumberComponents; i++)
  {
    int Number_operation = Decoder.get_bits(8); // Number of total operations.
    this->ComponentOperations.push_back(Number_operation);

    if(Number_operation > this->GlobalCountOperation)
      this->GlobalCountOperation = Number_operation;

    int Qbit = Decoder.get_bits(4); // Initial quantization bit of geometry
    Qbit += 4;
    this->Qbit.push_back(Qbit);

    int NCQ = Decoder.get_bits(4); // Number of change of geometry quantization
    this->NumberChangeQuantization.push_back(NCQ);

    int Number_color_quantization_change =
        Decoder.get_bits(3); // Number of change of color quantization
    this->NumberColorQuantization.push_back(Number_color_quantization_change);

    if(this->Qbit[i] > Max_Qbit)
      Max_Qbit = this->Qbit[i];
  }

  int Number_basemesh_vertex =
      Decoder.get_bits(15); // Number of vertices of base mesh
  int Number_basemesh_facet =
      Decoder.get_bits(16); // Number of facets of base mesh

  // Vectors for generation of base mesh
  std::vector< Point3d > vlist;
  std::vector< int > flist;
  std::vector< float > clist;
  std::vector< int > Color_index_list;

  for(int i = 0; i < Number_basemesh_vertex; i++)
  {
    Point_Int Pt_int;
    Pt_int.x = Decoder.get_bits(Max_Qbit + 1); // Read geometry info
    Pt_int.y = Decoder.get_bits(Max_Qbit + 1);
    Pt_int.z = Decoder.get_bits(Max_Qbit + 1);

    // All vertices have quantization precision of component 0
    // That'll be corrected below.
    Point3d Pt_real = Change_Int_Real(Pt_int, 0);
    vlist.push_back(Pt_real);

    if((this->IsColored) && (!this->IsOneColor))
    {
      //#ifdef PREDICTION_METHOD

      Color_Unit TC;
      TC.c0 = Decoder.get_bits(C0_QUANTIZATION); // Read color info.
      TC.c1 = Decoder.get_bits(C1_QUANTIZATION);
      TC.c2 = Decoder.get_bits(C2_QUANTIZATION);

      float L = this->C0_Min + TC.c0 * this->Color_Quantization_Step;
      float a = this->C1_Min + TC.c1 * this->Color_Quantization_Step;
      float b = this->C2_Min + TC.c2 * this->Color_Quantization_Step;

      clist.push_back(L);
      clist.push_back(a);
      clist.push_back(b);
      //	#endif
    }
  }

  // Read connectivity information
  int Facet_index_bit =
      (int)ceil(log((double)(Number_basemesh_vertex + 1)) / log((double)2));

  for(int i = 0; i < (Number_basemesh_facet * 3); i++)
  {
    int v = Decoder.get_bits(Facet_index_bit);
    flist.push_back(v);
  }

  // Generation of base mesh
  build_mesh(_pMesh, _pm, _v_cm, vlist, flist, clist, Color_index_list);

  compute_normals(_pMesh, _pm);

  // Seed Edges;
  std::map< int, int > Seed_Edges;

  for(int i = 0; i < 2 * this->NumberComponents; i++)
  {
    int Vertex_number = Decoder.get_bits(Facet_index_bit);
    Seed_Edges.insert(std::pair< int, int >(Vertex_number, i));
  }

  int Basemesh_vertex_number = 0;

  vertex_iterator pVertex;

  std::map< int, int >::iterator Seed_edge_iterator = Seed_Edges.begin();

  int Count_detected_vertices = 0;
  for(pVertex = vertices(_pMesh).first; pVertex != vertices(_pMesh).second;
      Basemesh_vertex_number++, pVertex++)
  {
    if(Count_detected_vertices < this->NumberComponents * 2)
    {
      if(Basemesh_vertex_number == Seed_edge_iterator->first)
      {
        put(this->vertex_Seed_Edge, *pVertex, Seed_edge_iterator->second);
        Seed_edge_iterator++;
        Count_detected_vertices++;
      }
    }
    else
      put(this->vertex_Seed_Edge, *pVertex, OTHER_COORDINATE);

    put(this->vertex_Component_Number, *pVertex, -1);
  }

  /*float C0_min = this->C0_Min;
  float C1_min = this->C1_Min;
  float C2_min = this->C2_Min;*/

  float Color_small_step = 0.0;

  if((this->IsColored) && (!this->IsOneColor))
  {
    // ELO  useless ?
    for(pVertex = vertices(_pMesh).first; pVertex != vertices(_pMesh).second;
        pVertex++)
    {
      /*float LAB[3];
      LAB[0] = pVertex->color(0);
      LAB[1] = pVertex->color(1);
      LAB[2] = pVertex->color(2);

      float RGB[3];
      LAB_To_RGB(LAB[0], LAB[1], LAB[2], RGB);
      pVertex->color(RGB[0], RGB[1], RGB[2]);		*/
    }

    this->C0_Range = Decoder.get_bits(C0_QUANTIZATION + 1);
    this->Color_0_Model.set_alphabet(this->C0_Range);
    this->C1_Range = Decoder.get_bits(C1_QUANTIZATION + 1);
    this->Color_1_Model.set_alphabet(this->C1_Range);
    this->C2_Range = Decoder.get_bits(C2_QUANTIZATION + 1);
    this->Color_2_Model.set_alphabet(this->C2_Range);
  }
  if((this->IsColored) && (this->IsOneColor))
  {
    for(pVertex = vertices(_pMesh).first; pVertex != vertices(_pMesh).second;
        pVertex++)
    {
      typedef
          typename boost::property_traits< VertexColorMap >::value_type Color;
      put(*_v_cm,
          *pVertex,
          Color(this->OnlyColor[0], this->OnlyColor[1], this->OnlyColor[2]));
    }
  }

  // To get number of vertices of each component and restore the real position
  // of vertices
  // if (this->NumberComponents != 1)
  {
    for(face_iterator pFacet = faces(_pMesh).first;
        pFacet != faces(_pMesh).second;
        ++pFacet)
    {
      put(this->facet_tag, *pFacet, -1);
    }
    int Component_index = 0;

    for(int Component_number = 0; Component_number < this->NumberComponents;
        Component_number++)
    {
      halfedge_iterator hi = halfedges(_pMesh).first;

      while(
          (get(this->vertex_Seed_Edge, target(*hi, _pMesh)) !=
           2 * Component_number) ||
          (get(this->vertex_Seed_Edge, target(opposite(*hi, _pMesh), _pMesh)) !=
           2 * Component_number + 1))
        hi++;

      face_descriptor fh = face(*hi, _pMesh);
      put(this->facet_tag, fh, Component_index);

      std::list< face_descriptor > facets;
      facets.push_front(fh);

      while(!facets.empty())
      {
        face_descriptor F = facets.front();
        facets.pop_front();

        put(this->facet_tag, F, Component_index);

        CGAL::Halfedge_around_face_circulator< HalfedgeGraph > pHalfedge =
            CGAL::Halfedge_around_face_circulator< HalfedgeGraph >(
                halfedge(F, _pMesh), _pMesh);
        CGAL::Halfedge_around_face_circulator< HalfedgeGraph > end = pHalfedge;

        CGAL_For_all(pHalfedge, end)
        {
          // tag the vertex to its corresponding component number
          if(get(this->vertex_Component_Number, target(*pHalfedge, _pMesh)) ==
             -1)
          {
            put(this->vertex_Component_Number,
                target(*pHalfedge, _pMesh),
                Component_index);
            Point3d Wrong_position = get(*_pm, target(*pHalfedge, _pMesh));

            // The correct position of vertex is restored.
            Point_Int Temp_pos = Change_Real_Int(Wrong_position, 0);
            Point3d Real_position = Change_Int_Real(Temp_pos, Component_number);

            put(*_pm, target(*pHalfedge, _pMesh), Real_position);

            float Wrong_LAB[3];
            Wrong_LAB[0] = get(*_v_cm, target(*pHalfedge, _pMesh))[0];
            Wrong_LAB[1] = get(*_v_cm, target(*pHalfedge, _pMesh))[1];
            Wrong_LAB[2] = get(*_v_cm, target(*pHalfedge, _pMesh))[2];

            int Original_LAB[3];
            Original_LAB[0] = (int)floor((Wrong_LAB[0] - this->C0_Min) /
                                             this->Color_Quantization_Step +
                                         0.5);
            Original_LAB[1] = (int)floor((Wrong_LAB[1] - this->C1_Min) /
                                             this->Color_Quantization_Step +
                                         0.5);
            Original_LAB[2] = (int)floor((Wrong_LAB[2] - this->C2_Min) /
                                             this->Color_Quantization_Step +
                                         0.5);

            put(this->vertex_color_int,
                target(*pHalfedge, _pMesh),
                Color_Unit(Original_LAB[0], Original_LAB[1], Original_LAB[2]));

            if(this->NumberColorQuantization[Component_number] == 0)
              Color_small_step = this->Color_Quantization_Step;
            else
              Color_small_step =
                  this->Color_Quantization_Step *
                  std::pow(2.0,
                           this->NumberColorQuantization[Component_number]);

            float RGB[3];
            float LAB[3];
            LAB[0] = this->C0_Min + Original_LAB[0] * Color_small_step;
            LAB[1] = this->C1_Min + Original_LAB[1] * Color_small_step;
            LAB[2] = this->C2_Min + Original_LAB[2] * Color_small_step;

            LAB_To_RGB(LAB[0], LAB[1], LAB[2], RGB);
            typedef
                typename boost::property_traits< VertexColorMap >::value_type
                    Color;
            put(*_v_cm,
                target(*pHalfedge, _pMesh),
                Color(RGB[0], RGB[1], RGB[2]));
          }

          face_descriptor pNFacet = face(opposite(*pHalfedge, _pMesh), _pMesh);
          if(pNFacet != boost::graph_traits< HalfedgeGraph >::null_face() &&
             get(this->facet_tag, pNFacet) == -1)
          {
            facets.push_front(pNFacet);
            put(this->facet_tag, pNFacet, Component_index);
          }
        }
      }
      // this->ComponentNumberVertices.push_back(Number_vertices);
      Component_index++;
    }
  }

  this->IsDecompress = true;
  this->Current_level = 0;
  // this->Initial_file_size = _Initial_file_size;
  this->Total_layer = this->GlobalCountOperation;


  float prog = (float)this->Calculate_Current_File_Size() /
               this->Compressed_file_size * 100;
  float ratio = 0;
  if(this->Initial_file_size != 0)
    ratio = 1 / ((float)this->Calculate_Current_File_Size() /
                 this->Initial_file_size);

  std::ostringstream infos_tmp;
  infos_tmp << "Number of all levels : ";
  infos_tmp << int(this->Total_layer);
  infos_tmp << "   |   ";
  infos_tmp << std::setw(3) << prog;
  infos_tmp << "   |   ";
  infos_tmp << std::setw(3) << ratio;
  std::string infos(infos_tmp.str());

  this->Prog.clear();
  this->Prog.push_back(prog);
  this->Ratio.clear();
  this->Ratio.push_back(ratio);

  std::cout << infos << std::endl;
}


//#define DBG_Decompress_Each_Step

// Description : To decode step by step - show intermediate meshes
template< typename HalfedgeGraph, typename PointMap, typename VertexColorMap >
int
Compression_Valence_Component< HalfedgeGraph, PointMap, VertexColorMap >::
    Decompress_Each_Step(HalfedgeGraph &_pMesh,
                         PointMap *_pm,
                         VertexColorMap *_v_cm)
{
  if(this->Decompress_count < this->GlobalCountOperation)
  {
    for(int Component_ID = 0; Component_ID < this->NumberComponents;
        Component_ID++)
    {
      if(this->Decompress_count < this->ComponentOperations[Component_ID])
      {
        int Operation = Decoder.get_bits(2);
        if(Operation == 0)
        {
#ifdef DBG_Decompress_Each_Step
          std::cout << "DBG in " << __func__
                    << "  Decompress_count=" << this->Decompress_count
                    << "  Component_ID=" << Component_ID
                    << "  Operation=" << Operation << std::endl;
          DBG_print_mesh_geometry(_pMesh, _pm, "before Un_Regulation");
          DBG_print_mesh_vertexcolor(_pMesh, _v_cm, "before Un_Regulation");
#endif

          this->Un_Regulation(_pMesh, Decoder, Component_ID, _pm, _v_cm);

#ifdef DBG_Decompress_Each_Step
          DBG_print_mesh_geometry(_pMesh, _pm, "after Un_Regulation");
          DBG_print_mesh_vertexcolor(_pMesh, _v_cm, "after Un_Regulation");
#endif

          this->Un_Decimation_Conquest(
              _pMesh, Decoder, Component_ID, _pm, _v_cm);

#ifdef DBG_Decompress_Each_Step
          DBG_print_mesh_geometry(_pMesh, _pm, "after Un_Decimation_Conquest");
          DBG_print_mesh_vertexcolor(
              _pMesh, _v_cm, "after Un_Decimation_Conquest");
#endif
        }
        else if(Operation == 1)
        {
          this->Augment_Geometry_Quantization_Precision(
              _pMesh, Decoder, Component_ID, _pm);

#ifdef DBG_Decompress_Each_Step
          std::cout << "DBG in " << __func__
                    << "  Decompress_count=" << this->Decompress_count
                    << "  Component_ID=" << Component_ID
                    << "  Operation=" << Operation << std::endl;
          DBG_print_mesh_geometry(
              _pMesh, _pm, "after Augment_Geometry_Quantization_Precision");
          DBG_print_mesh_vertexcolor(
              _pMesh, _v_cm, "after Augment_Geometry_Quantization_Precision");
#endif
        }
        else if(Operation == 2)
        {
          this->Augment_Color_Quantization_Precision(
              _pMesh, Decoder, Component_ID, _v_cm);

#ifdef DBG_Decompress_Each_Step
          std::cout << "DBG in " << __func__
                    << "  Decompress_count=" << this->Decompress_count
                    << "  Component_ID=" << Component_ID
                    << "  Operation=" << Operation << std::endl;
          DBG_print_mesh_geometry(
              _pMesh, _pm, "after Augment_Color_Quantization_Precision");
          DBG_print_mesh_vertexcolor(
              _pMesh, _v_cm, "after Augment_Color_Quantization_Precision");
#endif
        }
      }
    }
  }

  this->Decompress_count++;
  compute_normals(_pMesh, _pm);

  return this->Decompress_count;
  //
  /////////////////////////////////////////////////
  //                                             //
  //        WIP       WIP   WIP   WIPWIP         //
  //        WIP       WIP         WIP WIP        //
  //        WIP       WIP   WIP   WIPWIP         //
  //         WIP WIP WIP    WIP   WIP            //
  //           WIP WIP      WIP   WIP            //
  //                                             //
  /////////////////////////////////////////////////
  //
  //
  //
  // std::cout << "fixme here: " << __FILE__ << ":" << __LINE__ << " in " <<
  // __func__ << std::endl;
}


template< typename HalfedgeGraph, typename PointMap, typename VertexColorMap >
void
Compression_Valence_Component< HalfedgeGraph, PointMap, VertexColorMap >::
    Augment_Color_Quantization_Precision(HalfedgeGraph &_pMesh,
                                         Arithmetic_Codec &Decoder,
                                         const int &_Component_ID,
                                         VertexColorMap *_v_cm)
{
  Adaptive_Data_Model *Color_quantization_model =
      new Adaptive_Data_Model[COLOR_NUMBER];
  for(int i = 0; i < COLOR_NUMBER; i++)
    Color_quantization_model[i].set_alphabet(8);

  /*float C0_min = this->C0_Min;
  float C1_min = this->C1_Min;
  float C2_min = this->C2_Min;*/

  float Color_small_step = 0.0;
  float Color_large_step = 0.0;

  // define quantization step for (QC = large step) and (QC+1 = small_step)
  if(this->NumberColorQuantization[_Component_ID] == 0)
    Color_large_step = this->Color_Quantization_Step;
  else
    Color_large_step =
        this->Color_Quantization_Step *
        std::pow(2.0, this->NumberColorQuantization[_Component_ID]);

  Color_small_step = Color_large_step / (float)2.0;

  // To find first points to start the conquest.
  halfedge_iterator hi = halfedges(_pMesh).first;

  while(
      (get(this->vertex_Seed_Edge, target(*hi, _pMesh)) != 2 * _Component_ID) ||
      (get(this->vertex_Seed_Edge, target(opposite(*hi, _pMesh), _pMesh)) !=
       2 * _Component_ID + 1))
    hi++;

  // Vertex_Flag est donnee free a tous les sommets
  for(vertex_iterator pVert = vertices(_pMesh).first;
      pVert != vertices(_pMesh).second;
      pVert++)
  {
    put(this->Vertex_Flag, *pVert, FREE);
    put(this->Vertex_Number, *pVert, -1);
    put(this->vertex_Component_Number, *pVert, -1);
  }

  // Generation of color table
  std::vector< std::vector< int > > Color_table;

  std::queue< vertex_descriptor > verticesq;
  // push a input gate to begin loop
  // in order to treat all vertices.
  verticesq.push(target(*hi, _pMesh));
  verticesq.push(target(opposite(*hi, _pMesh), _pMesh));

  int Vertex_index = 0;

  while(!verticesq.empty())
  {
    vertex_descriptor v = verticesq.front();
    verticesq.pop();

    if(get(this->Vertex_Flag, v) == CONQUERED)
      continue;

    else
    {
      put(this->Vertex_Flag, v, CONQUERED);
      put(this->Vertex_Number, v, Vertex_index);
      put(this->vertex_Component_Number, v, _Component_ID);

      std::vector< int > LAB;

      Color_Unit color = get(this->vertex_color_int, v);
      LAB.push_back(color.c0);
      LAB.push_back(color.c1);
      LAB.push_back(color.c2);

      bool Is_existed_color = false;

      int Color_index = -1;

      // Color table creation
      for(int i = 0; i < (int)Color_table.size(); i++)
      {
        if((LAB[0] == Color_table[i][0]) && (LAB[1] == Color_table[i][1]) &&
           (LAB[2] == Color_table[i][2]))
        {
          Is_existed_color = true;
          Color_index = Decoder.decode(Color_quantization_model[i]);

          break;
        }
      }

      if(Is_existed_color == false)
      {
        Color_index =
            Decoder.decode(Color_quantization_model[Color_table.size()]);
        Color_table.push_back(LAB);
      }

      put(this->vertex_Q_Index, v, Color_index);

      // unsigned Count_neighbor = 0;
      // unsigned valence = v->vertex_degree();

      CGAL::Halfedge_around_target_circulator< HalfedgeGraph > h;

      // First_vertex -> To find the seed edge
      if(Vertex_index == 0)
      {
        h = CGAL::Halfedge_around_target_circulator< HalfedgeGraph >(v, _pMesh);
        do
        {
          h++;
        } while((*h) != (*hi));
      }

      // Synchronization proble = find an deterministic order using the vertex
      // with the highest vertex number
      else
      {
        int Comp_number = -2;
        h = CGAL::Halfedge_around_target_circulator< HalfedgeGraph >(v, _pMesh);
        CGAL::Halfedge_around_target_circulator< HalfedgeGraph > h2 = h;
        CGAL_For_all(h, h2)
        {
          if(get(this->Vertex_Number, target(opposite(*h, _pMesh), _pMesh)) >
             Comp_number)
            Comp_number =
                get(this->Vertex_Number, target(opposite(*h, _pMesh), _pMesh));
        }

        h = h2;
        CGAL_For_all(h, h2)
        {
          if(get(this->Vertex_Number, target(opposite(*h, _pMesh), _pMesh)) ==
             Comp_number)
            break;
        }
      }

      CGAL::Halfedge_around_target_circulator< HalfedgeGraph > h2 = h;
      CGAL_For_all(h, h2)
      {
        if(get(this->Vertex_Flag, target(opposite(*h, _pMesh), _pMesh)) == FREE)
          verticesq.push(target(opposite(*h, _pMesh), _pMesh));
      }
      Vertex_index++;
    }
  }

  // Here, we increase a bit of quantization precision of corresponding
  // component
  for(vertex_iterator pVert = vertices(_pMesh).first;
      pVert != vertices(_pMesh).second;
      pVert++)
  {
    if(get(this->vertex_Component_Number, *pVert) == _Component_ID)
    {
      int Coeff[3];
      Get_Coefficient_Up_Quantization(get(this->vertex_Q_Index, *pVert), Coeff);

      Color_Unit color = get(this->vertex_color_int, *pVert);
      int LAB_init_C0 = color.c0;
      int LAB_init_C1 = color.c1;
      int LAB_init_C2 = color.c2;

      // Each coordinate is multiplied by 2
      int LAB_final_C0 = LAB_init_C0 * 2;
      int LAB_final_C1 = LAB_init_C1 * 2;
      int LAB_final_C2 = LAB_init_C2 * 2;

      // We sum the childcell index to restore true color
      if(Coeff[0] == 1)
        LAB_final_C0 += 1;
      if(Coeff[1] == 1)
        LAB_final_C1 += 1;
      if(Coeff[2] == 1)
        LAB_final_C2 += 1;

      put(this->vertex_color_int,
          *pVert,
          Color_Unit(LAB_final_C0, LAB_final_C1, LAB_final_C2));

      float Inter_color[3];
      float RGB[3];

      Inter_color[0] = this->C0_Min + LAB_final_C0 * Color_small_step;
      Inter_color[1] = this->C1_Min + LAB_final_C1 * Color_small_step;
      Inter_color[2] = this->C2_Min + LAB_final_C2 * Color_small_step;

      LAB_To_RGB(Inter_color[0], Inter_color[1], Inter_color[2], RGB);

      typedef
          typename boost::property_traits< VertexColorMap >::value_type Color;
      put(*_v_cm, *pVert, Color(RGB[0], RGB[1], RGB[2]));
    }
  }

  // this->Q_color++;
  this->NumberColorQuantization[_Component_ID]--;

  delete[] Color_quantization_model;
}


// Description : Change a point coordinates in real to integer coordinates
template< typename HalfedgeGraph, typename PointMap, typename VertexColorMap >
Point_Int
Compression_Valence_Component< HalfedgeGraph, PointMap, VertexColorMap >::
    Change_Real_Int(const Point3d &pt, const int &Component_ID)
{
  Point_Int Point;

  double Quantization_step = 0.0;

  // If the quantization resolution is decreased,
  // we increase the step of quantization by a power of two.
  if(this->NumberChangeQuantization[Component_ID] == 0)
    Quantization_step = this->Quantization_Step[Component_ID];
  else
    Quantization_step =
        this->Quantization_Step[Component_ID] *
        std::pow(2.0, (int)this->NumberChangeQuantization[Component_ID]);

  float xmin = this->xmin[Component_ID];
  float ymin = this->ymin[Component_ID];
  float zmin = this->zmin[Component_ID];

  double x = pt[0];
  double y = pt[1];
  double z = pt[2];

  Point.x = (int)(ceil((x - xmin) / Quantization_step)) - 1;
  if(Point.x == -1)
    Point.x = 0;

  Point.y = (int)(ceil((y - ymin) / Quantization_step)) - 1;
  if(Point.y == -1)
    Point.y = 0;

  Point.z = (int)(ceil((z - zmin) / Quantization_step)) - 1;
  if(Point.z == -1)
    Point.z = 0;

  return Point;
}


// Description : Change a point coordinates in integer to real coordinates
template< typename HalfedgeGraph, typename PointMap, typename VertexColorMap >
typename Compression_Valence_Component< HalfedgeGraph,
                                        PointMap,
                                        VertexColorMap >::Point3d
Compression_Valence_Component< HalfedgeGraph, PointMap, VertexColorMap >::
    Change_Int_Real(const Point_Int &Point, const int &Component_ID)
{
  float Quantization_step = 0;

  // If the quantization resolution is decreased,
  // we increase the step of quantization by a power of two.
  if(this->NumberChangeQuantization[Component_ID] == 0)
    Quantization_step = this->Quantization_Step[Component_ID];
  else
    Quantization_step =
        this->Quantization_Step[Component_ID] *
        std::pow(2.0, (int)this->NumberChangeQuantization[Component_ID]);

  float xmin = this->xmin[Component_ID];
  float ymin = this->ymin[Component_ID];
  float zmin = this->zmin[Component_ID];

  double x = xmin + (Point.x + 0.5) * Quantization_step;
  double y = ymin + (Point.y + 0.5) * Quantization_step;
  double z = zmin + (Point.z + 0.5) * Quantization_step;

  Point3d pt(x, y, z);

  return pt;
}


//#define DBG_Error_Projected_Surface

// Error metric which measures importance of color and geometry for each vertex.
// Used to prevent removal of the visually important vertex.
template< typename HalfedgeGraph, typename PointMap, typename VertexColorMap >
bool
Compression_Valence_Component< HalfedgeGraph, PointMap, VertexColorMap >::
    Error_Projected_Surface(const HalfedgeGraph &_pMesh,
                            const PointMap *_pm,
                            const halfedge_descriptor &_h,
                            const int &_Component_ID,
                            const double &Mean_color,
                            const double &Mean_area)
{
  halfedge_descriptor g = _h;
  int Valence = (int)out_degree(target(next(g, _pMesh), _pMesh), _pMesh);
  int Type = Find_Type(_pMesh, g, Valence);

  float C0_min = this->C0_Min;
  float C1_min = this->C1_Min;
  float C2_min = this->C2_Min;

  float Color_step = 0.0;

  if(this->NumberColorQuantization[_Component_ID] == 0)
    Color_step = this->Color_Quantization_Step;
  else
    Color_step =
        this->Color_Quantization_Step *
        std::pow(2.0, (double)this->NumberColorQuantization[_Component_ID]);

  std::vector< float > Center_color;
  Center_color.push_back(
      C0_min +
      (get(this->vertex_color_int, target(next(g, _pMesh), _pMesh))).c0 *
          Color_step);
  Center_color.push_back(
      C1_min +
      (get(this->vertex_color_int, target(next(g, _pMesh), _pMesh))).c1 *
          Color_step);
  Center_color.push_back(
      C2_min +
      (get(this->vertex_color_int, target(next(g, _pMesh), _pMesh))).c2 *
          Color_step);

  std::vector< Color_Unit > Neighbors_color;

  std::vector< float > Projected_color;

  double Patch_area = 0.0;

  // Calculate area of patch
  g = next(g, _pMesh);
  for(int i = 0; i < Valence; i++)
  {
    g = opposite(g, _pMesh);
    Color_Unit c;

    c = get(this->vertex_color_int, target(g, _pMesh));
    Neighbors_color.push_back(c);

    Patch_area += Area_Facet_Triangle(g, _pMesh, _pm);

    g = prev(g, _pMesh);
  }

  Color_Unit Average_color;

  if(Valence == 3)
  {
    for(int i = 0; i < 3; i++)
      Average_color = Average_color + Neighbors_color[i];
  }

  if(Valence == 4)
  {
    if((Type == 5) || (Type == 8))
    {
      Average_color = Average_color + Neighbors_color[1];
      Average_color = Average_color + Neighbors_color[3];
    }

    if((Type == 6) || (Type == 7))
    {
      Average_color = Average_color + Neighbors_color[0];
      Average_color = Average_color + Neighbors_color[2];
    }
  }

  if(Valence == 5)
  {
    if((Type == 9) || (Type == 12))
    {
      Average_color = Average_color + Neighbors_color[1];
      Average_color = Average_color + Neighbors_color[3];
      Average_color = Average_color + Neighbors_color[4];
    }

    if(Type == 10)
    {
      Average_color = Average_color + Neighbors_color[0];
      Average_color = Average_color + Neighbors_color[1];
      Average_color = Average_color + Neighbors_color[3];
    }

    if(Type == 11)
    {
      Average_color = Average_color + Neighbors_color[0];
      Average_color = Average_color + Neighbors_color[2];
      Average_color = Average_color + Neighbors_color[4];
    }
  }


  if(Valence == 6)
  {
    if((Type == 13) || (Type == 16))
    {
      Average_color = Average_color + Neighbors_color[1];
      Average_color = Average_color + Neighbors_color[3];
      Average_color = Average_color + Neighbors_color[5];
    }
    else
    {
      Average_color = Average_color + Neighbors_color[0];
      Average_color = Average_color + Neighbors_color[2];
      Average_color = Average_color + Neighbors_color[4];
    }
  }

  float Divised_c0, Divised_c1, Divised_c2;
  if(Valence != 4)
  {
    Divised_c0 = (float)Average_color.c0 / 3.0;
    Divised_c1 = (float)Average_color.c1 / 3.0;
    Divised_c2 = (float)Average_color.c2 / 3.0;
  }
  else
  {
    Divised_c0 = (float)Average_color.c0 / 2.0;
    Divised_c1 = (float)Average_color.c1 / 2.0;
    Divised_c2 = (float)Average_color.c2 / 2.0;
  }

  Projected_color.push_back(C0_min + Divised_c0 * Color_step);
  Projected_color.push_back(C1_min + Divised_c1 * Color_step);
  Projected_color.push_back(C2_min + Divised_c2 * Color_step);


  // Color distance between the original color and estimated color
  double Color_distance = (Projected_color[0] - Center_color[0]) *
                              (Projected_color[0] - Center_color[0]) +
                          (Projected_color[1] - Center_color[1]) *
                              (Projected_color[1] - Center_color[1]) +
                          (Projected_color[2] - Center_color[2]) *
                              (Projected_color[2] - Center_color[2]);

  Color_distance = std::sqrt(Color_distance);

  // Mean color is obtained using number of vertices.
  // X 3.0 cause we should use the number of edges.
  double Relative_color_distance = Color_distance / Mean_color * 3.0;
#ifdef DBG_Error_Projected_Surface
  std::cout << "Color_distance  = " << Color_distance << std::endl;
  std::cout << "Mean_color  = " << Mean_color << std::endl;
  std::cout << "Relative_color_distance = " << Relative_color_distance
            << std::endl;
#endif

  // Averaged area of triangles of patch
  double Area_per_triangle = Patch_area / double(Valence);
  Area_per_triangle *= std::pow(
      ((double)10.0 / (double)this->HighestLengthBB[_Component_ID]), 2.0);

  double Relative_geo_distance = Area_per_triangle / Mean_area;
#ifdef DBG_Error_Projected_Surface
  std::cout << "Area_per_triangle  = " << Area_per_triangle << std::endl;
  std::cout << "Mean_area  = " << Mean_area << std::endl;
  std::cout << "Relative_geo_distance = " << Relative_geo_distance << std::endl;
#endif

  double Global_distance = (Relative_color_distance) * (Relative_geo_distance);
#ifdef DBG_Error_Projected_Surface
  std::cout << "Global_distance = " << Global_distance << std::endl;
#endif

  if(Global_distance > 0.5)
    return true;
  else
    return false;
}


template< typename HalfedgeGraph, typename PointMap, typename VertexColorMap >
void
Compression_Valence_Component< HalfedgeGraph, PointMap, VertexColorMap >::
    write_intermediate_mesh(/*const*/ HalfedgeGraph &_pMesh,
                            const VertexColorMap *_v_cm)
{
  // build output file name
  std::string outputFilePath = this->File_name;
  outputFilePath += ".level" + std::to_string(this->Current_level) + ".off";

  // create a temporary property map bag in order
  // to be able to call the generic writer
  FEVV::PMapsContainer pmaps_bag;
  put_property_map(FEVV::vertex_color, _pMesh, pmaps_bag, *_v_cm);

  // write the mesh to file
  FEVV::Filters::write_mesh(outputFilePath, _pMesh, pmaps_bag);
}


template< typename HalfedgeGraph, typename PointMap, typename VertexColorMap >
void
Compression_Valence_Component< HalfedgeGraph, PointMap, VertexColorMap >::
    copy_mesh(const HalfedgeGraph &_pMesh,
              const VertexColorMap *_v_cm,
              HalfedgeGraph &mesh_copy,
              VertexColorMap *v_cm_copy /* nullptr allowed */)
{
  // storage for vertices correspondance
  boost::unordered_map< vertex_descriptor, vertex_descriptor > v2v;

  // copy topology and geometry
  CGAL::clear(mesh_copy);
  CGAL::copy_face_graph(_pMesh, mesh_copy, std::inserter(v2v, v2v.end()));

  // copy the vertex-color property map
  if(v_cm_copy)
  {
    typename boost::unordered_map< vertex_descriptor,
                                   vertex_descriptor >::iterator
        it = v2v.begin(),
        end = v2v.end();
    for(; it != end; ++it)
    {
      auto color = get(*_v_cm, it->first);
      put(*v_cm_copy, it->second, color);
    }
  }
}


template< typename HalfedgeGraph, typename PointMap, typename VertexColorMap >
void
Compression_Valence_Component< HalfedgeGraph, PointMap, VertexColorMap >::
    keep_intermediate_mesh(
        const HalfedgeGraph &_pMesh,
        const VertexColorMap *_v_cm,
        std::vector< HalfedgeGraph * > *intermediate_meshes,
        std::vector< VertexColorMap * > *intermediate_vertexColorMaps)
{
  if(!intermediate_meshes)
    return;

  // create empty new mesh
  HalfedgeGraph *m_copy = new HalfedgeGraph;

  // create empty new vertex-color property map
  VertexColorMap *v_cm_copy = nullptr;
  if(intermediate_vertexColorMaps)
    v_cm_copy = new VertexColorMap;

  // copy mesh + geometry + vertex-color
  copy_mesh(_pMesh, _v_cm, *m_copy, v_cm_copy);

  // store copied mesh and vertex-color
  intermediate_meshes->push_back(m_copy);
  if(intermediate_vertexColorMaps)
    intermediate_vertexColorMaps->push_back(v_cm_copy);
}


//#define DBG_Decompression_All_From_File

template< typename HalfedgeGraph, typename PointMap, typename VertexColorMap >
std::string
Compression_Valence_Component< HalfedgeGraph, PointMap, VertexColorMap >::
    Decompression_All_From_File(
        HalfedgeGraph &_pMesh,
        PointMap *_pm,
        VertexColorMap *_v_cm,
        bool do_write_info,
        std::vector< HalfedgeGraph * >
            *intermediate_meshes /* nullptr allowed */,
        std::vector< VertexColorMap * >
            *intermediate_vertexColorMaps /* nullptr allowed */,
        int stop_level /* decompression level to stop at */,
        bool do_write_intermediate_meshes)
{
  // start time measurement
  auto time_start = std::chrono::steady_clock::now();

  if(do_write_info)
  {
    if(this->Process_level == 0)
      Write_Info(_pMesh);
  }

  //ELO+ fix color > 1 issue on level 0 mesh
  truncate_colors(_pMesh, _v_cm);

  // save the level 0 mesh (the simplest one) to file
  if(do_write_intermediate_meshes)
    write_intermediate_mesh(_pMesh, _v_cm);

  // keep the level 0 mesh (the simplest one)
  if(intermediate_meshes)
    keep_intermediate_mesh(
        _pMesh, _v_cm, intermediate_meshes, intermediate_vertexColorMaps);

  // check at which level to stop decompression
  if(stop_level == -1  ||  stop_level > this->Total_layer)
    stop_level = this->Total_layer;

  // loop over decompression levels
  while(this->Current_level != stop_level)
  {
#ifdef DBG_Decompression_All_From_File
    std::cout << "DBG in " << __func__
              << "  Current_level=" << this->Current_level << std::endl;
    DBG_print_mesh_geometry(_pMesh, _pm, "loop begin");
    DBG_print_mesh_vertexcolor(_pMesh, _v_cm, "loop begin");
#endif

    this->Current_level = this->Decompress_Each_Step(_pMesh, _pm, _v_cm);
    this->Process_level++;

    if(do_write_info)
    {
      Write_Info(_pMesh);
    }

    //ELO+ fix color > 1 issue on current level mesh
    truncate_colors(_pMesh, _v_cm);

    // save intermediate meshes into files
    if(do_write_intermediate_meshes)
      write_intermediate_mesh(_pMesh, _v_cm);

    // keep intermediate meshes
    if(intermediate_meshes && this->Current_level < stop_level)
      keep_intermediate_mesh(
          _pMesh, _v_cm, intermediate_meshes, intermediate_vertexColorMaps);


#ifdef DBG_Decompression_All_From_File
    DBG_print_mesh_geometry(_pMesh, _pm, "loop end");
    DBG_print_mesh_vertexcolor(_pMesh, _v_cm, "loop end");
#endif
  }

  // stop time measurement
  auto time_end = std::chrono::steady_clock::now();
  std::chrono::duration< double > time_diff = time_end - time_start;
  std::stringstream msgbuffer;
  msgbuffer << "Number of layers : " << this->Total_layer << std::endl;
  msgbuffer << "Calculation time for decompressing layers 1 to "
            << stop_level << " : " << (float)time_diff.count()
            << " seconds" << std::endl;
  msgbuffer << "Time for loading .p3d and decompressing layer 0 was not taken "
               "into account."
            << std::endl;
  std::cout << msgbuffer.str() << std::endl;

  return msgbuffer.str();
}


template< typename HalfedgeGraph, typename PointMap, typename VertexColorMap >
void
Compression_Valence_Component< HalfedgeGraph, PointMap, VertexColorMap >::
    Write_Info(HalfedgeGraph &_pMesh)
{
  if(this->Process_level == 0)
  {
    this->Dec_File_Info = this->File_name;
    (this->Dec_File_Info).append(".infos.txt");

    this->Dec_Info = fopen(this->Dec_File_Info.c_str(), "w");
  }
  else
    this->Dec_Info = fopen(this->Dec_File_Info.c_str(), "a");

  int CLevel = 0, Number_vertices = 0;

  if(this->Sequence)
  {
    CLevel = this->Visu_level;
    Number_vertices = (int)FEVV::size_of_vertices(_pMesh);
  }
  else
  {
    CLevel = this->Current_level;
    Number_vertices = (int)FEVV::size_of_vertices(_pMesh);
  }

  unsigned Current_file_size = this->Calculate_Current_File_Size();

  float prog = (float)Current_file_size / this->Compressed_file_size * 100;

  if(this->Process_level == this->Total_layer)
    prog = 100.0;

  float ratio = 1 / ((float)Current_file_size / this->Initial_file_size);

  fprintf(this->Dec_Info,
          "Level %2d   #v : %8d      %6u bytes     Prog : %7.3f %%    Ratio : "
          "%9.3f\n",
          CLevel,
          Number_vertices,
          Current_file_size,
          prog,
          ratio);
  fclose(this->Dec_Info);
}


//-----------------------------------------------------------
//    Functions moved from Compression_Valence_Common.h
//-----------------------------------------------------------

template< typename HalfedgeGraph, typename PointMap, typename VertexColorMap >
bool
Compression_Valence_Component< HalfedgeGraph, PointMap, VertexColorMap >::
    Is_Border_Vertex(const halfedge_descriptor &h, const HalfedgeGraph &mesh)
{
  bool check = false;

  CGAL::Halfedge_around_target_circulator< HalfedgeGraph > hvc(target(h, mesh),
                                                               mesh);
  CGAL::Halfedge_around_target_circulator< HalfedgeGraph > hvc_end(hvc);

  do
  {
    if(CGAL::is_border_edge(*hvc, mesh))
    {
      check = true;
      break;
    }
  } while(++hvc != hvc_end);

  return check;
}


template< typename HalfedgeGraph, typename PointMap, typename VertexColorMap >
int
Compression_Valence_Component< HalfedgeGraph, PointMap, VertexColorMap >::
    Find_Type(const HalfedgeGraph &_pMesh,
              const halfedge_descriptor &h,
              const int &valence)
{
  int type = 0;

  if(valence == 3)
  {
    if((get(this->Vertex_Sign, target(h, _pMesh)) == MINUS) &&
       (get(this->Vertex_Sign, target(opposite(h, _pMesh), _pMesh)) == PLUS))
      type = 1;
    else if((get(this->Vertex_Sign, target(h, _pMesh)) == PLUS) &&
            (get(this->Vertex_Sign, target(opposite(h, _pMesh), _pMesh)) ==
             MINUS))
      type = 2;
    else if((get(this->Vertex_Sign, target(h, _pMesh)) == PLUS) &&
            (get(this->Vertex_Sign, target(opposite(h, _pMesh), _pMesh)) ==
             PLUS))
      type = 3;
    else if((get(this->Vertex_Sign, target(h, _pMesh)) == MINUS) &&
            (get(this->Vertex_Sign, target(opposite(h, _pMesh), _pMesh)) ==
             MINUS))
      type = 4;
  }

  else if(valence == 4)
  {
    if((get(this->Vertex_Sign, target(h, _pMesh)) == MINUS) &&
       (get(this->Vertex_Sign, target(opposite(h, _pMesh), _pMesh)) == PLUS))
      type = 5;
    else if((get(this->Vertex_Sign, target(h, _pMesh)) == PLUS) &&
            (get(this->Vertex_Sign, target(opposite(h, _pMesh), _pMesh)) ==
             MINUS))
      type = 6;
    else if((get(this->Vertex_Sign, target(h, _pMesh)) == PLUS) &&
            (get(this->Vertex_Sign, target(opposite(h, _pMesh), _pMesh)) ==
             PLUS))
      type = 7;
    else if((get(this->Vertex_Sign, target(h, _pMesh)) == MINUS) &&
            (get(this->Vertex_Sign, target(opposite(h, _pMesh), _pMesh)) ==
             MINUS))
      type = 8;
  }

  else if(valence == 5)
  {
    if((get(this->Vertex_Sign, target(h, _pMesh)) == MINUS) &&
       (get(this->Vertex_Sign, target(opposite(h, _pMesh), _pMesh)) == PLUS))
      type = 9;
    else if((get(this->Vertex_Sign, target(h, _pMesh)) == PLUS) &&
            (get(this->Vertex_Sign, target(opposite(h, _pMesh), _pMesh)) ==
             MINUS))
      type = 10;
    else if((get(this->Vertex_Sign, target(h, _pMesh)) == PLUS) &&
            (get(this->Vertex_Sign, target(opposite(h, _pMesh), _pMesh)) ==
             PLUS))
      type = 11;
    else if((get(this->Vertex_Sign, target(h, _pMesh)) == MINUS) &&
            (get(this->Vertex_Sign, target(opposite(h, _pMesh), _pMesh)) ==
             MINUS))
      type = 12;
  }

  else if(valence == 6)
  {
    if((get(this->Vertex_Sign, target(h, _pMesh)) == MINUS) &&
       (get(this->Vertex_Sign, target(opposite(h, _pMesh), _pMesh)) == PLUS))
      type = 13;
    else if((get(this->Vertex_Sign, target(h, _pMesh)) == PLUS) &&
            (get(this->Vertex_Sign, target(opposite(h, _pMesh), _pMesh)) ==
             MINUS))
      type = 14;
    else if((get(this->Vertex_Sign, target(h, _pMesh)) == PLUS) &&
            (get(this->Vertex_Sign, target(opposite(h, _pMesh), _pMesh)) ==
             PLUS))
      type = 15;
    else if((get(this->Vertex_Sign, target(h, _pMesh)) == MINUS) &&
            (get(this->Vertex_Sign, target(opposite(h, _pMesh), _pMesh)) ==
             MINUS))
      type = 16;
  }
  return type;
}


template< typename HalfedgeGraph, typename PointMap, typename VertexColorMap >
bool
Compression_Valence_Component< HalfedgeGraph, PointMap, VertexColorMap >::
    Is_Manifold_Property_Violated(const HalfedgeGraph &_pMesh,
                                  const halfedge_descriptor &h,
                                  const int &type,
                                  const int &valence)
{
  bool check = false;
  halfedge_descriptor g = h;
  int *Points_index = new int[valence];

  // if valence is 3, no new edge is inserted, so always safe to remove.
  if(valence == 3)
    return false;

  else
  {
    // Points_index[  ] contains all boundary vertices' indices (ordered in
    // counterclockwise)
    Points_index[0] = get(this->Vertex_Number, target(g, _pMesh));
    g = next(g, _pMesh); // g points center vertex;

    for(int i = 1; i < valence; i++)
    {
      g = prev(opposite(g, _pMesh),
               _pMesh); // around the vertex in the counterclockwise way.
      Points_index[i] =
          get(this->Vertex_Number, target(opposite(g, _pMesh), _pMesh));
    }

    // quadrangle
    if(valence == 4)
    {
      if((type == 5) || (type == 8))
      {
        g = opposite(h, _pMesh);
        CGAL::Halfedge_around_target_circulator< HalfedgeGraph > Hvc(
            target(g, _pMesh), _pMesh);
        CGAL::Halfedge_around_target_circulator< HalfedgeGraph > Hvc_end = Hvc;

        CGAL_For_all(Hvc, Hvc_end)
        {
          if(get(this->Vertex_Number, target(opposite(*Hvc, _pMesh), _pMesh)) ==
             Points_index[1])
            check = true;
        }
      }

      else if((type == 6) || (type == 7))
      {
        g = h;
        CGAL::Halfedge_around_target_circulator< HalfedgeGraph > Hvc(
            target(g, _pMesh), _pMesh);
        CGAL::Halfedge_around_target_circulator< HalfedgeGraph > Hvc_end = Hvc;

        CGAL_For_all(Hvc, Hvc_end)
        {
          if(get(this->Vertex_Number, target(opposite(*Hvc, _pMesh), _pMesh)) ==
             Points_index[2])
            check = true;
          ;
        }
      }
    }

    // pendtagone : 2 edges to verify
    if(valence == 5)
    {
      if((type == 9) || (type == 12))
      {
        g = opposite(h, _pMesh);
        CGAL::Halfedge_around_target_circulator< HalfedgeGraph > Hvc(
            target(g, _pMesh), _pMesh);
        CGAL::Halfedge_around_target_circulator< HalfedgeGraph > Hvc_end = Hvc;

        CGAL_For_all(Hvc, Hvc_end)
        {
          if(get(this->Vertex_Number, target(opposite(*Hvc, _pMesh), _pMesh)) ==
             Points_index[1])
            check = true;
        }

        g = next(opposite(next(h, _pMesh), _pMesh), _pMesh);
        Hvc = CGAL::Halfedge_around_target_circulator< HalfedgeGraph >(
            target(g, _pMesh), _pMesh);
        Hvc_end = Hvc;

        CGAL_For_all(Hvc, Hvc_end)
        {
          if(get(this->Vertex_Number, target(opposite(*Hvc, _pMesh), _pMesh)) ==
             Points_index[3])
            check = true;
        }
      }

      else if(type == 10)
      {
        g = h;
        CGAL::Halfedge_around_target_circulator< HalfedgeGraph > Hvc(
            target(g, _pMesh), _pMesh);
        CGAL::Halfedge_around_target_circulator< HalfedgeGraph > Hvc_end = Hvc;

        CGAL_For_all(Hvc, Hvc_end)
        {
          if(get(this->Vertex_Number, target(opposite(*Hvc, _pMesh), _pMesh)) ==
             Points_index[3])
            check = true;
        }

        g = next(opposite(next(h, _pMesh), _pMesh), _pMesh);
        Hvc = CGAL::Halfedge_around_target_circulator< HalfedgeGraph >(
            target(g, _pMesh), _pMesh);
        Hvc_end = Hvc;

        CGAL_For_all(Hvc, Hvc_end)
        {
          if(get(this->Vertex_Number, target(opposite(*Hvc, _pMesh), _pMesh)) ==
             Points_index[3])
            check = true;
        }
      }

      else if(type == 11)
      {
        g = h;
        CGAL::Halfedge_around_target_circulator< HalfedgeGraph > Hvc(
            target(g, _pMesh), _pMesh);
        CGAL::Halfedge_around_target_circulator< HalfedgeGraph > Hvc_end = Hvc;

        CGAL_For_all(Hvc, Hvc_end)
        {
          if(get(this->Vertex_Number, target(opposite(*Hvc, _pMesh), _pMesh)) ==
             Points_index[2])
            check = true;
        }

        g = opposite(h, _pMesh);
        Hvc = CGAL::Halfedge_around_target_circulator< HalfedgeGraph >(
            target(g, _pMesh), _pMesh);
        Hvc_end = Hvc;

        CGAL_For_all(Hvc, Hvc_end)
        {
          if(get(this->Vertex_Number, target(opposite(*Hvc, _pMesh), _pMesh)) ==
             Points_index[2])
            check = true;
        }
      }
    }

    // hexagone

    if(valence == 6)
    {
      if((type == 13) || (type == 16))
      {
        g = opposite(h, _pMesh);
        CGAL::Halfedge_around_target_circulator< HalfedgeGraph > Hvc(
            target(g, _pMesh), _pMesh);
        CGAL::Halfedge_around_target_circulator< HalfedgeGraph > Hvc_end = Hvc;

        CGAL_For_all(Hvc, Hvc_end)
        {
          if(get(this->Vertex_Number, target(opposite(*Hvc, _pMesh), _pMesh)) ==
             Points_index[1])
            check = true;
        }

        g = opposite(h, _pMesh);
        Hvc = CGAL::Halfedge_around_target_circulator< HalfedgeGraph >(
            target(g, _pMesh), _pMesh);
        Hvc_end = Hvc;

        CGAL_For_all(Hvc, Hvc_end)
        {
          if(get(this->Vertex_Number, target(opposite(*Hvc, _pMesh), _pMesh)) ==
             Points_index[3])
            check = true;
        }

        g = next(opposite(next(h, _pMesh), _pMesh), _pMesh);
        Hvc = CGAL::Halfedge_around_target_circulator< HalfedgeGraph >(
            target(g, _pMesh), _pMesh);
        Hvc_end = Hvc;

        CGAL_For_all(Hvc, Hvc_end)
        {
          if(get(this->Vertex_Number, target(opposite(*Hvc, _pMesh), _pMesh)) ==
             Points_index[3])
            check = true;
        }
      }

      else if((type == 14) || (type == 15))
      {
        g = h;
        CGAL::Halfedge_around_target_circulator< HalfedgeGraph > Hvc(
            target(g, _pMesh), _pMesh);
        CGAL::Halfedge_around_target_circulator< HalfedgeGraph > Hvc_end = Hvc;

        CGAL_For_all(Hvc, Hvc_end)
        {
          if(get(this->Vertex_Number, target(opposite(*Hvc, _pMesh), _pMesh)) ==
             Points_index[2])
            check = true;
        }

        g = h;
        Hvc = CGAL::Halfedge_around_target_circulator< HalfedgeGraph >(
            target(g, _pMesh), _pMesh);
        Hvc_end = Hvc;

        CGAL_For_all(Hvc, Hvc_end)
        {
          if(get(this->Vertex_Number, target(opposite(*Hvc, _pMesh), _pMesh)) ==
             Points_index[4])
            check = true;
        }

        g = next(opposite(prev(h, _pMesh), _pMesh), _pMesh);
        Hvc = CGAL::Halfedge_around_target_circulator< HalfedgeGraph >(
            target(g, _pMesh), _pMesh);
        Hvc_end = Hvc;

        CGAL_For_all(Hvc, Hvc_end)
        {
          if(get(this->Vertex_Number, target(opposite(*Hvc, _pMesh), _pMesh)) ==
             Points_index[2])
            check = true;
        }
      }
    }
  }

  delete[] Points_index;
  return check;
}

template< typename HalfedgeGraph, typename PointMap, typename VertexColorMap >
bool
Compression_Valence_Component< HalfedgeGraph, PointMap, VertexColorMap >::
    Is_Normal_Flipping_Occured(const HalfedgeGraph &_pMesh,
                               const PointMap *_pm,
                               const halfedge_descriptor &h,
                               const unsigned &valence)
{
  typedef typename FEVV::Geometry_traits< HalfedgeGraph >::Point Point3d;
  typedef typename FEVV::Geometry_traits< HalfedgeGraph >::Vector Vector;
  FEVV::Geometry_traits< HalfedgeGraph > gt(_pMesh);

  int type = Find_Type(_pMesh, h, valence);
  bool check = false;
  ////////////////////////////////////////////////////////////////////////////////////////////////////////////////

  halfedge_descriptor g = h;
  g = next(g, _pMesh);

  Vector V_normal = get(this->vertex_normal, target(next(h, _pMesh), _pMesh));

  double length =
      FEVV::Math::Vector::l2_distance< FEVV::Geometry_traits< HalfedgeGraph > >(
          V_normal, gt);
  if(length != 0)
    V_normal = V_normal / length;

  Point3d *Points = new Point3d[valence];
  Vector *Normal = new Vector[valence - 2];


  for(unsigned int i = 0; i < valence; i++)
  {
    Points[i] = get(*_pm, target(opposite(g, _pMesh), _pMesh));
    g = prev(opposite(g, _pMesh), _pMesh);
  }

  for(unsigned int j = 0; j < (valence - 2); j++)
  {
    Normal[j] = Vector(0.0, 0.0, 0.0);
  }

  // if valence = 3 , there is no normal flipping;
  if(valence == 3)
  {
    check = false;
  }

  // quadrangle
  else if((type == 5) || (type == 8)) // 0 1 3 , 1 2 3
  {
    Normal[0] = Triangle_Normal(_pMesh, Points[0], Points[1], Points[3]);
    Normal[1] = Triangle_Normal(_pMesh, Points[1], Points[2], Points[3]);
  }

  else if((type == 6) || (type == 7)) // 0 1 2 , 0 2 3
  {
    Normal[0] = Triangle_Normal(_pMesh, Points[0], Points[1], Points[2]);
    Normal[1] = Triangle_Normal(_pMesh, Points[0], Points[2], Points[3]);
  }

  // pentagone
  else if((type == 9) || (type == 12)) // 0 1 4 , 1 2 3 , 1 3 4
  {
    Normal[0] = Triangle_Normal(_pMesh, Points[0], Points[1], Points[4]);
    Normal[1] = Triangle_Normal(_pMesh, Points[1], Points[2], Points[3]);
    Normal[2] = Triangle_Normal(_pMesh, Points[1], Points[3], Points[4]);
  }


  else if(type == 10) // 0 1 3 , 1 2 3 , 0 3 4
  {
    Normal[0] = Triangle_Normal(_pMesh, Points[0], Points[1], Points[3]);
    Normal[1] = Triangle_Normal(_pMesh, Points[1], Points[2], Points[3]);
    Normal[2] = Triangle_Normal(_pMesh, Points[0], Points[3], Points[4]);
  }

  else if(type == 11) // 0 1 2 , 0 2 4 , 2 3 4
  {
    Normal[0] = Triangle_Normal(_pMesh, Points[0], Points[1], Points[2]);
    Normal[1] = Triangle_Normal(_pMesh, Points[0], Points[2], Points[4]);
    Normal[2] = Triangle_Normal(_pMesh, Points[2], Points[3], Points[4]);
  }

  // Hexagone

  else if((type == 13) || (type == 16)) // 0 1 5,  1 2 3 , 3 4 5, 1 3 5
  {
    Normal[0] = Triangle_Normal(_pMesh, Points[0], Points[1], Points[5]);
    Normal[1] = Triangle_Normal(_pMesh, Points[1], Points[2], Points[3]);
    Normal[2] = Triangle_Normal(_pMesh, Points[3], Points[4], Points[5]);
    Normal[3] = Triangle_Normal(_pMesh, Points[1], Points[3], Points[5]);
  }

  else if((type == 14) || (type == 15)) // 0 1 2 , 2 3 4 , 4 5 0, 0 2 4
  {
    Normal[0] = Triangle_Normal(_pMesh, Points[0], Points[1], Points[2]);
    Normal[1] = Triangle_Normal(_pMesh, Points[2], Points[3], Points[4]);
    Normal[2] = Triangle_Normal(_pMesh, Points[4], Points[5], Points[0]);
    Normal[3] = Triangle_Normal(_pMesh, Points[0], Points[2], Points[4]);
  }

  for(unsigned int i = 0; i < (valence - 2); i++)
  {
    double length_normal = FEVV::Math::Vector::l2_distance<
        FEVV::Geometry_traits< HalfedgeGraph > >(Normal[i], gt);
    if(length_normal != 0)
      Normal[i] = Normal[i] / length_normal;
  }

  for(unsigned int i = 0; i < (valence - 2); i++)
  {
    double cosine = FEVV::Math::Vector::dot_product<
        FEVV::Geometry_traits< HalfedgeGraph > >(V_normal, Normal[i], gt);
    double cosine_rad = std::acos(cosine);

    if(cosine_rad >= PI / 2)
      check = true;
  }
  delete[] Points;
  delete[] Normal;

  return check;
}


template< typename HalfedgeGraph, typename PointMap, typename VertexColorMap >
typename FEVV::Geometry_traits< HalfedgeGraph >::Vector
Compression_Valence_Component< HalfedgeGraph, PointMap, VertexColorMap >::
    Triangle_Normal(const HalfedgeGraph &_pMesh,
                    const PointMap *_pm,
                    const halfedge_descriptor &h)
{
  typedef typename FEVV::Geometry_traits< HalfedgeGraph >::Point Point3d;
  typedef typename FEVV::Geometry_traits< HalfedgeGraph >::Vector Vector;
  FEVV::Geometry_traits< HalfedgeGraph > gt(_pMesh);

  Point3d P = get(*_pm, target(h, _pMesh));
  Point3d Q = get(*_pm, target(next(h, _pMesh), _pMesh));
  Point3d R = get(*_pm, target(next(next(h, _pMesh), _pMesh), _pMesh));

  Vector PQ = Q - P;
  // Vector PR = R-P; // MT
  Vector QR = R - Q;

  Vector normal = FEVV::Math::Vector::cross_product<
      FEVV::Geometry_traits< HalfedgeGraph > >(PQ, QR, gt);
  double length =
      FEVV::Math::Vector::l2_distance< FEVV::Geometry_traits< HalfedgeGraph > >(
          normal, gt);
  if(length != 0.0)
    normal = normal / length;

  return normal;
}


template< typename HalfedgeGraph, typename PointMap, typename VertexColorMap >
typename FEVV::Geometry_traits< HalfedgeGraph >::Vector
Compression_Valence_Component< HalfedgeGraph, PointMap, VertexColorMap >::
    Triangle_Normal(
        const HalfedgeGraph &_pMesh,
        const typename FEVV::Geometry_traits< HalfedgeGraph >::Point &P,
        const typename FEVV::Geometry_traits< HalfedgeGraph >::Point &Q,
        const typename FEVV::Geometry_traits< HalfedgeGraph >::Point &R)
{
  typedef typename FEVV::Geometry_traits< HalfedgeGraph >::Point Point3d;
  typedef typename FEVV::Geometry_traits< HalfedgeGraph >::Vector Vector;
  FEVV::Geometry_traits< HalfedgeGraph > gt(_pMesh);

  Vector PQ = Q - P;
  // Vector PR = R - P; // MT
  Vector QR = R - Q;

  Vector normal = FEVV::Math::Vector::cross_product<
      FEVV::Geometry_traits< HalfedgeGraph > >(PQ, QR, gt);
  double length =
      FEVV::Math::Vector::l2_distance< FEVV::Geometry_traits< HalfedgeGraph > >(
          normal, gt);
  if(length != 0.0)
    normal = normal / length;

  return normal;
}


template< typename HalfedgeGraph, typename PointMap, typename VertexColorMap >
bool
Compression_Valence_Component< HalfedgeGraph, PointMap, VertexColorMap >::
    Is_Geometric_Metric_Violated(const HalfedgeGraph &_pMesh,
                                 const PointMap *_pm,
                                 const halfedge_descriptor &h,
                                 const int &type,
                                 const unsigned int &valence,
                                 const float &Threshold)
{
  FEVV::Geometry_traits< HalfedgeGraph > gt(_pMesh);

  bool check = false;
  double volume = 0;    // volume;;
  double perimeter = 0; // perimeter;;;
  double area = 0;      // area;;
  unsigned int i;

  halfedge_descriptor g = h;
  Point3d *Points = new Point3d[valence + 1];
  g = next(g, _pMesh); // g points front vertex

  // Points[0] = coordinates of the front vertex
  // Points[i] = coordinates of the ith vertex in the counterclockwise way.
  Points[0] = get(*_pm, target(g, _pMesh));
  for(i = 1; i < (valence + 1); i++)
  {
    Points[i] = get(*_pm, target(opposite(g, _pMesh), _pMesh));
    g = prev(opposite(g, _pMesh), _pMesh);
  }

  // caculate perimeter
  Vector *Vectors = new Vector[valence];
  Vectors[0] = Points[1] - Points[valence];
  perimeter =
      FEVV::Math::Vector::l2_distance< FEVV::Geometry_traits< HalfedgeGraph > >(
          Vectors[0], gt);

  for(i = 1; i < valence; i++)
  {
    Vectors[i] = Points[(i + 1)] - Points[i];
    double per = FEVV::Math::Vector::l2_distance<
        FEVV::Geometry_traits< HalfedgeGraph > >(Vectors[i], gt);
    perimeter = perimeter + per;
  }

  // calculate volume and area;;;
  if(valence == 3) //[0 1 2 3]
  {
    volume = Volume(Points[0], Points[1], Points[3], Points[2], _pMesh);
    volume = std::abs(volume);

    area = Area_Facet_Triangle(Points[1], Points[2], Points[3], _pMesh);
    area = std::abs(area);
  }

  else if(valence == 4)
  {
    if((type == 5) || (type == 8)) // [0 1 2 4], [0 2 3 4]
    {
      double vol = Volume(Points[0], Points[1], Points[4], Points[2], _pMesh);
      vol = std::abs(vol);
      volume = volume + vol;

      vol = Volume(Points[0], Points[2], Points[4], Points[3], _pMesh);
      vol = std::abs(vol);
      volume = volume + vol;

      double area_intermediate =
          Area_Facet_Triangle(Points[1], Points[2], Points[4], _pMesh);
      area_intermediate = std::abs(area_intermediate);
      area = area + area_intermediate;

      area_intermediate =
          Area_Facet_Triangle(Points[2], Points[3], Points[4], _pMesh);
      area_intermediate = std::abs(area_intermediate);
      area = area + area_intermediate;
    }

    else if((type == 6) || (type == 7)) // [0 1 2 3], [0 1 3 4]
    {
      double vol = Volume(Points[0], Points[1], Points[3], Points[2], _pMesh);
      vol = std::abs(vol);
      volume = volume + vol;

      vol = Volume(Points[0], Points[1], Points[4], Points[3], _pMesh);
      vol = std::abs(vol);
      volume = volume + vol;

      double area_intermediate =
          Area_Facet_Triangle(Points[1], Points[2], Points[3], _pMesh);
      area_intermediate = std::abs(area_intermediate);
      area = area + area_intermediate;

      area_intermediate =
          Area_Facet_Triangle(Points[1], Points[3], Points[4], _pMesh);
      area_intermediate = std::abs(area_intermediate);
      area = area + area_intermediate;
    }
  }

  else if(valence == 5)
  {
    if((type == 9) || (type == 12)) // [0 1 2 5] , [0 2 4 5], [0 2 3 4]
    {
      double vol = Volume(Points[0], Points[1], Points[5], Points[2], _pMesh);
      vol = std::abs(vol);
      volume = volume + vol;

      vol = Volume(Points[0], Points[2], Points[5], Points[4], _pMesh);
      vol = std::abs(vol);
      volume = volume + vol;

      vol = Volume(Points[0], Points[2], Points[4], Points[3], _pMesh);
      vol = std::abs(vol);
      volume = volume + vol;

      double area_intermediate =
          Area_Facet_Triangle(Points[1], Points[2], Points[5], _pMesh);
      area_intermediate = std::abs(area_intermediate);
      area = area + area_intermediate;

      area_intermediate =
          Area_Facet_Triangle(Points[2], Points[4], Points[5], _pMesh);
      area_intermediate = std::abs(area_intermediate);
      area = area + area_intermediate;

      area_intermediate =
          Area_Facet_Triangle(Points[2], Points[3], Points[4], _pMesh);
      area_intermediate = std::abs(area_intermediate);
      area = area + area_intermediate;
    }

    else if(type == 10) // [0 1 4 5], [0 1 2 4], [0 2 3 4]
    {
      double vol = Volume(Points[0], Points[1], Points[5], Points[4], _pMesh);
      vol = std::abs(vol);
      volume = volume + vol;

      vol = Volume(Points[0], Points[1], Points[4], Points[2], _pMesh);
      vol = std::abs(vol);
      volume = volume + vol;

      vol = Volume(Points[0], Points[2], Points[4], Points[3], _pMesh);
      vol = std::abs(vol);
      volume = volume + vol;

      double area_intermediate =
          Area_Facet_Triangle(Points[1], Points[4], Points[5], _pMesh);
      area_intermediate = std::abs(area_intermediate);
      area = area + area_intermediate;

      area_intermediate =
          Area_Facet_Triangle(Points[1], Points[2], Points[4], _pMesh);
      area_intermediate = std::abs(area_intermediate);
      area = area + area_intermediate;

      area_intermediate =
          Area_Facet_Triangle(Points[2], Points[3], Points[4], _pMesh);
      area_intermediate = std::abs(area_intermediate);
      area = area + area_intermediate;
    }

    else if(type == 11) // [0 1 2 3], [0 1 3 5], [0 3 4 5]
    {
      double vol = Volume(Points[0], Points[1], Points[3], Points[2], _pMesh);
      vol = std::abs(vol);
      volume = volume + vol;

      vol = Volume(Points[0], Points[1], Points[5], Points[3], _pMesh);
      vol = std::abs(vol);
      volume = volume + vol;

      vol = Volume(Points[0], Points[3], Points[5], Points[4], _pMesh);
      vol = std::abs(vol);
      volume = volume + vol;

      double area_intermediate =
          Area_Facet_Triangle(Points[1], Points[2], Points[3], _pMesh);
      area_intermediate = std::abs(area_intermediate);
      area = area + area_intermediate;

      area_intermediate =
          Area_Facet_Triangle(Points[1], Points[3], Points[5], _pMesh);
      area_intermediate = std::abs(area_intermediate);
      area = area + area_intermediate;

      area_intermediate =
          Area_Facet_Triangle(Points[3], Points[4], Points[5], _pMesh);
      area_intermediate = std::abs(area_intermediate);
      area = area + area_intermediate;
    }
  }

  else if(valence == 6)
  {
    if((type == 13) ||
       (type == 16)) //[0 1 2 6], [0 2 3 4] , [0 4 5 6], [0 2 4 6]
    {
      double vol = Volume(Points[0], Points[1], Points[6], Points[2], _pMesh);
      vol = std::abs(vol);
      volume = volume + vol;

      vol = Volume(Points[0], Points[2], Points[4], Points[3], _pMesh);
      vol = std::abs(vol);
      volume = volume + vol;

      vol = Volume(Points[0], Points[4], Points[6], Points[5], _pMesh);
      vol = std::abs(vol);
      volume = volume + vol;

      vol = Volume(Points[0], Points[2], Points[6], Points[4], _pMesh);
      vol = std::abs(vol);
      volume = volume + vol;

      double area_intermediate =
          Area_Facet_Triangle(Points[1], Points[2], Points[6], _pMesh);
      area_intermediate = std::abs(area_intermediate);
      area = area + area_intermediate;

      area_intermediate =
          Area_Facet_Triangle(Points[2], Points[3], Points[4], _pMesh);
      area_intermediate = std::abs(area_intermediate);
      area = area + area_intermediate;

      area_intermediate =
          Area_Facet_Triangle(Points[4], Points[5], Points[6], _pMesh);
      area_intermediate = std::abs(area_intermediate);
      area = area + area_intermediate;

      area_intermediate =
          Area_Facet_Triangle(Points[2], Points[6], Points[4], _pMesh);
      area_intermediate = std::abs(area_intermediate);
      area = area + area_intermediate;
    }

    else if((type == 14) ||
            (type == 15)) // [ 0 1 2 3] [ 0 3 4 5] [ 0 1 5 6] [ 0 1 3 5]
    {
      double vol = Volume(Points[0], Points[1], Points[3], Points[2], _pMesh);
      vol = std::abs(vol);
      volume = volume + vol;

      vol = Volume(Points[0], Points[3], Points[5], Points[4], _pMesh);
      vol = std::abs(vol);
      volume = volume + vol;

      vol = Volume(Points[0], Points[1], Points[6], Points[5], _pMesh);
      vol = std::abs(vol);
      volume = volume + vol;

      vol = Volume(Points[0], Points[1], Points[5], Points[3], _pMesh);
      vol = std::abs(vol);
      volume = volume + vol;

      double area_intermediate =
          Area_Facet_Triangle(Points[1], Points[2], Points[3], _pMesh);
      area_intermediate = std::abs(area_intermediate);
      area = area + area_intermediate;

      area_intermediate =
          Area_Facet_Triangle(Points[3], Points[4], Points[5], _pMesh);
      area_intermediate = std::abs(area_intermediate);
      area = area + area_intermediate;

      area_intermediate =
          Area_Facet_Triangle(Points[1], Points[5], Points[6], _pMesh);
      area_intermediate = std::abs(area_intermediate);
      area = area + area_intermediate;

      area_intermediate =
          Area_Facet_Triangle(Points[1], Points[3], Points[5], _pMesh);
      area_intermediate = std::abs(area_intermediate);
      area = area + area_intermediate;
    }
  }

  double volume_3_root = 0.0;

  double volume_measured = 0.0;

  if(volume == 0)
    volume_3_root = volume;
  else
    volume_3_root = exp(1.0 / 3.0 * (double)log(volume));


  // volume_mea
  volume_measured = volume_3_root / (perimeter / (double)valence);

  // if Volume calculated is greater than threshold value, we should give up
  // removal of this vertex
  if(volume_measured > Threshold)
    check = true;

  // tp free memory.
  delete[] Points;
  delete[] Vectors;

  return check;
}


template< typename HalfedgeGraph, typename PointMap, typename VertexColorMap >
double
Compression_Valence_Component< HalfedgeGraph, PointMap, VertexColorMap >::
    Area_Facet_Triangle(
        const typename boost::graph_traits< HalfedgeGraph >::halfedge_descriptor
            &h,
        const HalfedgeGraph &_pMesh,
        const PointMap *_pm)
{
  typedef typename FEVV::Geometry_traits< HalfedgeGraph >::Point Point3d;
  FEVV::Geometry_traits< HalfedgeGraph > gt(_pMesh);

  Point3d P = get(*_pm, target(h, _pMesh));
  Point3d Q = get(*_pm, target(next(h, _pMesh), _pMesh));
  Point3d R = get(*_pm, target(next(next(h, _pMesh), _pMesh), _pMesh));

  return FEVV::Operators::Geometry::triangle_area<
      FEVV::Geometry_traits< HalfedgeGraph > >(P, Q, R, gt);
}

template< typename HalfedgeGraph, typename PointMap, typename VertexColorMap >
double
Compression_Valence_Component< HalfedgeGraph, PointMap, VertexColorMap >::
    Area_Facet_Triangle(const Point3d &P,
                        const Point3d &Q,
                        const Point3d &R,
                        const HalfedgeGraph &_pMesh)
{
  FEVV::Geometry_traits< HalfedgeGraph > gt(_pMesh);

  Vector PQ = FEVV::Math::Vector::sub< FEVV::Geometry_traits< HalfedgeGraph > >(
      Q, P, gt);
  // Vector PR = R - P; // MT
  Vector QR = FEVV::Math::Vector::sub< FEVV::Geometry_traits< HalfedgeGraph > >(
      R, Q, gt);

  Vector normal = FEVV::Math::Vector::cross_product<
      FEVV::Geometry_traits< HalfedgeGraph > >(PQ, QR, gt);
  double area =
      0.5 *
      FEVV::Math::Vector::l2_distance< FEVV::Geometry_traits< HalfedgeGraph > >(
          normal, gt);

  return area;
}


template< typename HalfedgeGraph, typename PointMap, typename VertexColorMap >
double
Compression_Valence_Component< HalfedgeGraph, PointMap, VertexColorMap >::
    Volume(const Point3d &A,
           const Point3d &B,
           const Point3d &C,
           const Point3d &D,
           const HalfedgeGraph &_pMesh)
{
  FEVV::Geometry_traits< HalfedgeGraph > gt(_pMesh);

  //
  // from https://fr.wikipedia.org/wiki/T%C3%A9tra%C3%A8dre
  //    V = (1/6) * det(AB, AC, AD)
  //
  // from https://fr.wikipedia.org/wiki/D%C3%A9terminant_(math%C3%A9matiques)
  //    det(X,X',X") = x.y'.z" + x'.y".z + x".y.z' - x.y".z' - x".y'.z - x'.y.z"
  //

  Vector AB = FEVV::Math::Vector::sub< FEVV::Geometry_traits< HalfedgeGraph > >(
      B, A, gt);
  Vector AC = FEVV::Math::Vector::sub< FEVV::Geometry_traits< HalfedgeGraph > >(
      C, A, gt);
  Vector AD = FEVV::Math::Vector::sub< FEVV::Geometry_traits< HalfedgeGraph > >(
      D, A, gt);

  double x1 = AB[0];
  double y1 = AB[1];
  double z1 = AB[2];

  double x2 = AC[0];
  double y2 = AC[1];
  double z2 = AC[2];

  double x3 = AD[0];
  double y3 = AD[1];
  double z3 = AD[2];

  double vol = x1 * y2 * z3 + x2 * y3 * z1 + x3 * y1 * z2 - x1 * y3 * z2 -
               x3 * y2 * z1 - x2 * y1 * z3;

  return vol;
}


template< typename HalfedgeGraph, typename PointMap, typename VertexColorMap >
Color_Unit
Compression_Valence_Component< HalfedgeGraph, PointMap, VertexColorMap >::
    Get_Vertex_Color(const halfedge_descriptor &h, const HalfedgeGraph &_pMesh)
{
  return get(this->vertex_color_int, target(h, _pMesh));
}


template< typename HalfedgeGraph, typename PointMap, typename VertexColorMap >
Color_Unit
Compression_Valence_Component< HalfedgeGraph, PointMap, VertexColorMap >::
    Get_Average_Vertex_Color_Lee(const halfedge_descriptor &h,
                                 const int &valence,
                                 const HalfedgeGraph &_pMesh)
{
  halfedge_descriptor g = h;

  std::vector< Color_Unit >
      Neighbor_color; // contain differents colors among neighbors
  int Average_c0 = 0, Average_c1 = 0, Average_c2 = 0;

  g = next(g, _pMesh);

  for(int i = 0; i < valence; i++)
  {
    Color_Unit col = Get_Vertex_Color(opposite(g, _pMesh), _pMesh);

    bool check = false;
    for(unsigned j = 0; j < Neighbor_color.size(); j++)
    {
      if(Neighbor_color[j] == col)
      {
        check = true;
        break;
      }
    }

    if(!check)
      Neighbor_color.push_back(col);

    Average_c0 += col.c0;
    Average_c1 += col.c1;
    Average_c2 += col.c2;

    g = prev(opposite(g, _pMesh), _pMesh);
  }

  Color_Unit Average_color;

  Average_color.c0 = (int)floor((float)(Average_c0 / valence) + 0.5);
  Average_color.c1 = (int)floor((float)(Average_c1 / valence) + 0.5);
  Average_color.c2 = (int)floor((float)(Average_c2 / valence) + 0.5);

  double Critere_c0 = 50000.0, Critere_c1 = 50000.0, Critere_c2 = 50000.0;
  int Index_c0 = -1, Index_c1 = -1, Index_c2 = -1;

  // Find for each component among neighbors' color, the component which is the
  // nearst to the average value.
  for(unsigned i = 0; i < Neighbor_color.size(); i++)
  {
    int Color_component_0 = Neighbor_color[i].c0;
    int Color_component_1 = Neighbor_color[i].c1;
    int Color_component_2 = Neighbor_color[i].c2;

    if(std::abs(Color_component_0 - Average_color.c0) < Critere_c0)
    {
      Critere_c0 = std::abs(Color_component_0 - Average_color.c0);
      Index_c0 = i;
    }
    if(std::abs(Color_component_1 - Average_color.c1) < Critere_c1)
    {
      Critere_c1 = std::abs(Color_component_1 - Average_color.c1);
      Index_c1 = i;
    }
    if(std::abs(Color_component_2 - Average_color.c2) < Critere_c2)
    {
      Critere_c2 = std::abs(Color_component_2 - Average_color.c2);
      Index_c2 = i;
    }
  }

  Color_Unit Resulting_color;
  Resulting_color.c0 = Neighbor_color[Index_c0].c0;
  Resulting_color.c1 = Neighbor_color[Index_c1].c1;
  Resulting_color.c2 = Neighbor_color[Index_c2].c2;

  return Resulting_color;
}


template< typename HalfedgeGraph, typename PointMap, typename VertexColorMap >
typename Compression_Valence_Component< HalfedgeGraph,
                                        PointMap,
                                        VertexColorMap >::Point3d
Compression_Valence_Component< HalfedgeGraph, PointMap, VertexColorMap >::
    Barycenter_Patch_Before_Removal(const halfedge_descriptor &h,
                                    const HalfedgeGraph &_pMesh,
                                    const PointMap *_pm)
{
  halfedge_descriptor g = h;
  int valence = (int)out_degree(target(next(g, _pMesh), _pMesh), _pMesh);

  double x = 0., y = 0., z = 0.;

  CGAL::Halfedge_around_target_circulator< HalfedgeGraph > vh_it(
      target(next(g, _pMesh), _pMesh), _pMesh);
  CGAL::Halfedge_around_target_circulator< HalfedgeGraph > vh_it_end = vh_it;

  CGAL_For_all(vh_it, vh_it_end)
  {
    Point3d pt = get(*_pm, target(opposite(*vh_it, _pMesh), _pMesh));
    x += pt[0];
    y += pt[1];
    z += pt[2];
  }

  x = x / (double)valence;
  y = y / (double)valence;
  z = z / (double)valence;

  return Point3d(x, y, z);
}


template< typename HalfedgeGraph, typename PointMap, typename VertexColorMap >
typename Compression_Valence_Component< HalfedgeGraph,
                                        PointMap,
                                        VertexColorMap >::Point3d
Compression_Valence_Component< HalfedgeGraph, PointMap, VertexColorMap >::
    Barycenter_Patch_After_Removal(const halfedge_descriptor &h,
                                   const int &valence,
                                   const HalfedgeGraph &_pMesh,
                                   const PointMap *_pm)
{
  halfedge_descriptor g = h;
  double x = 0, y = 0, z = 0;
  for(int i = 0; i < valence; i++)
  {
    Point3d pt = get(*_pm, target(g, _pMesh));
    x += pt[0];
    y += pt[1];
    z += pt[2];

    g = next(g, _pMesh);
  }

  x /= valence;
  y /= valence;
  z /= valence;

  return Point3d(x, y, z);
}


//#define DBG_Retriangulation

template< typename HalfedgeGraph, typename PointMap, typename VertexColorMap >
void
Compression_Valence_Component< HalfedgeGraph, PointMap, VertexColorMap >::
    Retriangulation(HalfedgeGraph &pMesh,
                    const halfedge_descriptor &ch,
                    const unsigned &valence,
                    const unsigned &Vertex_number,
                    const int &Component_ID,
                    const PointMap *_pm)
{
  int type = Find_Type(pMesh, ch, valence);
  halfedge_descriptor h;
  h = ch;
  halfedge_descriptor g;
  g = next(h, pMesh);


  // Triangle
  if((type == 1) || (type == 2) || (type == 4))
  {
    put(this->Facet_Flag, face(h, pMesh), TO_BE_REMOVED);
    // h->facet()->Patch_Index = Vertex_number;

    put(this->facet_normal, face(h, pMesh), Triangle_Normal(pMesh, _pm, h));

    put(this->facet_Component_Number, face(h, pMesh), Component_ID);

    if(get(this->Vertex_Sign, target(g, pMesh)) == NOSIGN)
      put(this->Vertex_Sign, target(g, pMesh), PLUS);
  }
  else if(type == 3)
  {
    put(this->Facet_Flag, face(h, pMesh), TO_BE_REMOVED);
    // h->facet()->Patch_Index = Vertex_number;

    put(this->facet_Component_Number, face(h, pMesh), Component_ID);

    put(this->facet_normal, face(h, pMesh), Triangle_Normal(pMesh, _pm, h));

    if(get(this->Vertex_Sign, target(g, pMesh)) == NOSIGN)
      put(this->Vertex_Sign, target(g, pMesh), MINUS);
  }

  // quadrangle
  else if((type == 5) || (type == 8))
  {
    if(get(this->Vertex_Sign, target(g, pMesh)) == NOSIGN)
      put(this->Vertex_Sign, target(g, pMesh), PLUS);
    if(get(this->Vertex_Sign, target(next(g, pMesh), pMesh)) == NOSIGN)
      put(this->Vertex_Sign, target(next(g, pMesh), pMesh), MINUS);

    h = prev(h, pMesh);
#ifdef DBG_Retriangulation
    std::cout << "split face" << std::endl;
#endif
    h = CGAL::Euler::split_face(h, g, pMesh);

    put(this->facet_Component_Number, face(h, pMesh), Component_ID);
    put(this->facet_Component_Number,
        face(opposite(h, pMesh), pMesh),
        Component_ID);

    put(this->facet_normal, face(h, pMesh), Triangle_Normal(pMesh, _pm, h));
    put(this->facet_normal,
        face(opposite(h, pMesh), pMesh),
        Triangle_Normal(pMesh, _pm, opposite(h, pMesh)));


    put(this->Facet_Flag, face(h, pMesh), TO_BE_REMOVED);
    put(this->Facet_Flag, face(opposite(h, pMesh), pMesh), TO_BE_REMOVED);
  }

  else if((type == 6) || (type == 7))
  {
    if(get(this->Vertex_Sign, target(g, pMesh)) == NOSIGN)
      put(this->Vertex_Sign, target(g, pMesh), MINUS);
    if(get(this->Vertex_Sign, target(next(g, pMesh), pMesh)) == NOSIGN)
      put(this->Vertex_Sign, target(next(g, pMesh), pMesh), PLUS);

    g = next(g, pMesh);
#ifdef DBG_Retriangulation
    std::cout << "split face" << std::endl;
#endif
    h = CGAL::Euler::split_face(h, g, pMesh);

    put(this->facet_Component_Number, face(h, pMesh), Component_ID);
    put(this->facet_Component_Number,
        face(opposite(h, pMesh), pMesh),
        Component_ID);

    put(this->facet_normal, face(h, pMesh), Triangle_Normal(pMesh, _pm, h));
    put(this->facet_normal,
        face(opposite(h, pMesh), pMesh),
        Triangle_Normal(pMesh, _pm, opposite(h, pMesh)));

    put(this->Facet_Flag, face(h, pMesh), TO_BE_REMOVED);
    put(this->Facet_Flag, face(opposite(h, pMesh), pMesh), TO_BE_REMOVED);
    // h->opposite()->facet()->Patch_Index = Vertex_number;
  }
  // pentagone
  else if((type == 9) || (type == 12))
  {
    if(get(this->Vertex_Sign, target(g, pMesh)) == NOSIGN)
      put(this->Vertex_Sign, target(g, pMesh), PLUS);
    if(get(this->Vertex_Sign, target(next(g, pMesh), pMesh)) == NOSIGN)
      put(this->Vertex_Sign, target(next(g, pMesh), pMesh), MINUS);
    if(get(this->Vertex_Sign, target(next(next(g, pMesh), pMesh), pMesh)) ==
       NOSIGN)
      put(this->Vertex_Sign, target(next(next(g, pMesh), pMesh), pMesh), PLUS);

    h = prev(h, pMesh);

#ifdef DBG_Retriangulation
    std::cout << "split face" << std::endl;
#endif
    h = CGAL::Euler::split_face(h, g, pMesh);

    put(this->facet_Component_Number, face(h, pMesh), Component_ID);
    put(this->facet_Component_Number,
        face(opposite(h, pMesh), pMesh),
        Component_ID);

    put(this->Facet_Flag, face(opposite(h, pMesh), pMesh), TO_BE_REMOVED);
    // h->opposite()->facet()->Patch_Index = Vertex_number;
    put(this->facet_normal,
        face(opposite(h, pMesh), pMesh),
        Triangle_Normal(pMesh, _pm, opposite(h, pMesh)));

    g = next(next(h, pMesh), pMesh);
#ifdef DBG_Retriangulation
    std::cout << "split face" << std::endl;
#endif
    h = CGAL::Euler::split_face(h, g, pMesh);

    put(this->facet_Component_Number, face(h, pMesh), Component_ID);
    put(this->facet_Component_Number,
        face(opposite(h, pMesh), pMesh),
        Component_ID);

    put(this->facet_normal, face(h, pMesh), Triangle_Normal(pMesh, _pm, h));
    put(this->facet_normal,
        face(opposite(h, pMesh), pMesh),
        Triangle_Normal(pMesh, _pm, opposite(h, pMesh)));

    put(this->Facet_Flag, face(h, pMesh), TO_BE_REMOVED);
    // h->facet()->Patch_Index = Vertex_number;
    put(this->Facet_Flag, face(opposite(h, pMesh), pMesh), TO_BE_REMOVED);
    // h->opposite()->facet()->Patch_Index = Vertex_number;
  }

  else if(type == 10)
  {
    if(get(this->Vertex_Sign, target(g, pMesh)) == NOSIGN)
      put(this->Vertex_Sign, target(g, pMesh), PLUS);
    if(get(this->Vertex_Sign, target(next(g, pMesh), pMesh)) == NOSIGN)
      put(this->Vertex_Sign, target(next(g, pMesh), pMesh), MINUS);
    if(get(this->Vertex_Sign, target(next(next(g, pMesh), pMesh), pMesh)) ==
       NOSIGN)
      put(this->Vertex_Sign, target(next(next(g, pMesh), pMesh), pMesh), PLUS);

    g = h;
    h = prev(prev(h, pMesh), pMesh);

#ifdef DBG_Retriangulation
    std::cout << "split face" << std::endl;
#endif
    h = CGAL::Euler::split_face(h, g, pMesh);

    put(this->facet_Component_Number, face(h, pMesh), Component_ID);
    put(this->facet_Component_Number,
        face(opposite(h, pMesh), pMesh),
        Component_ID);

    put(this->Facet_Flag, face(opposite(h, pMesh), pMesh), TO_BE_REMOVED);
    // h->opposite()->facet()->Patch_Index = Vertex_number;
    put(this->facet_normal,
        face(opposite(h, pMesh), pMesh),
        Triangle_Normal(pMesh, _pm, opposite(h, pMesh)));

    g = next(h, pMesh);
    h = prev(h, pMesh);

#ifdef DBG_Retriangulation
    std::cout << "split face" << std::endl;
#endif
    h = CGAL::Euler::split_face(h, g, pMesh);

    put(this->facet_Component_Number, face(h, pMesh), Component_ID);
    put(this->facet_Component_Number,
        face(opposite(h, pMesh), pMesh),
        Component_ID);

    put(this->facet_normal, face(h, pMesh), Triangle_Normal(pMesh, _pm, h));
    put(this->facet_normal,
        face(opposite(h, pMesh), pMesh),
        Triangle_Normal(pMesh, _pm, opposite(h, pMesh)));

    put(this->Facet_Flag, face(h, pMesh), TO_BE_REMOVED);
    // h->facet()->Patch_Index = Vertex_number;
    put(this->Facet_Flag, face(opposite(h, pMesh), pMesh), TO_BE_REMOVED);
    // h->opposite()->facet()->Patch_Index = Vertex_number;
  }

  else if(type == 11)
  {
    if(get(this->Vertex_Sign, target(g, pMesh)) == NOSIGN)
      put(this->Vertex_Sign, target(g, pMesh), MINUS);
    if(get(this->Vertex_Sign, target(next(g, pMesh), pMesh)) == NOSIGN)
      put(this->Vertex_Sign, target(next(g, pMesh), pMesh), PLUS);
    if(get(this->Vertex_Sign, target(next(next(g, pMesh), pMesh), pMesh)) ==
       NOSIGN)
      put(this->Vertex_Sign, target(next(next(g, pMesh), pMesh), pMesh), MINUS);

    g = next(g, pMesh);
#ifdef DBG_Retriangulation
    std::cout << "split face" << std::endl;
#endif
    h = CGAL::Euler::split_face(h, g, pMesh);

    put(this->facet_Component_Number, face(h, pMesh), Component_ID);
    put(this->facet_Component_Number,
        face(opposite(h, pMesh), pMesh),
        Component_ID);

    put(this->Facet_Flag, face(opposite(h, pMesh), pMesh), TO_BE_REMOVED);
    // h->opposite()->facet()->Patch_Index = Vertex_number;

    put(this->facet_normal,
        face(opposite(h, pMesh), pMesh),
        Triangle_Normal(pMesh, _pm, opposite(h, pMesh)));
    g = next(next(h, pMesh), pMesh);

#ifdef DBG_Retriangulation
    std::cout << "split face" << std::endl;
#endif
    h = CGAL::Euler::split_face(h, g, pMesh);

    put(this->facet_Component_Number, face(h, pMesh), Component_ID);
    put(this->facet_Component_Number,
        face(opposite(h, pMesh), pMesh),
        Component_ID);

    put(this->facet_normal, face(h, pMesh), Triangle_Normal(pMesh, _pm, h));
    put(this->facet_normal,
        face(opposite(h, pMesh), pMesh),
        Triangle_Normal(pMesh, _pm, opposite(h, pMesh)));

    put(this->Facet_Flag, face(h, pMesh), TO_BE_REMOVED);
    // h->facet()->Patch_Index = Vertex_number;
    put(this->Facet_Flag, face(opposite(h, pMesh), pMesh), TO_BE_REMOVED);
    // h->opposite()->facet()->Patch_Index = Vertex_number;
  }

  // Hexagone
  else if((type == 13) || (type == 16))
  {
    if(get(this->Vertex_Sign, target(g, pMesh)) == NOSIGN)
      put(this->Vertex_Sign, target(g, pMesh), PLUS);
    if(get(this->Vertex_Sign, target(next(g, pMesh), pMesh)) == NOSIGN)
      put(this->Vertex_Sign, target(next(g, pMesh), pMesh), MINUS);
    if(get(this->Vertex_Sign, target(next(next(g, pMesh), pMesh), pMesh)) ==
       NOSIGN)
      put(this->Vertex_Sign, target(next(next(g, pMesh), pMesh), pMesh), PLUS);
    if(get(this->Vertex_Sign,
           target(next(next(next(g, pMesh), pMesh), pMesh), pMesh)) == NOSIGN)
      put(this->Vertex_Sign,
          target(next(next(next(g, pMesh), pMesh), pMesh), pMesh),
          MINUS);

    h = prev(h, pMesh);
#ifdef DBG_Retriangulation
    std::cout << "split face" << std::endl;
#endif
    h = CGAL::Euler::split_face(h, g, pMesh);

    put(this->facet_Component_Number, face(h, pMesh), Component_ID);
    put(this->facet_Component_Number,
        face(opposite(h, pMesh), pMesh),
        Component_ID);

    put(this->facet_normal,
        face(opposite(h, pMesh), pMesh),
        Triangle_Normal(pMesh, _pm, opposite(h, pMesh)));
    put(this->Facet_Flag, face(opposite(h, pMesh), pMesh), TO_BE_REMOVED);
    // h->opposite()->facet()->Patch_Index = Vertex_number;


    g = next(next(h, pMesh), pMesh);
#ifdef DBG_Retriangulation
    std::cout << "split face" << std::endl;
#endif
    h = CGAL::Euler::split_face(h, g, pMesh);

    put(this->facet_Component_Number, face(h, pMesh), Component_ID);
    put(this->facet_Component_Number,
        face(opposite(h, pMesh), pMesh),
        Component_ID);

    put(this->Facet_Flag, face(opposite(h, pMesh), pMesh), TO_BE_REMOVED);
    // h->opposite()->facet()->Patch_Index = Vertex_number;

    put(this->facet_normal,
        face(opposite(h, pMesh), pMesh),
        Triangle_Normal(pMesh, _pm, opposite(h, pMesh)));

    g = next(next(h, pMesh), pMesh);
#ifdef DBG_Retriangulation
    std::cout << "split face" << std::endl;
#endif
    h = CGAL::Euler::split_face(h, g, pMesh);

    put(this->facet_Component_Number, face(h, pMesh), Component_ID);
    put(this->facet_Component_Number,
        face(opposite(h, pMesh), pMesh),
        Component_ID);

    put(this->facet_normal, face(h, pMesh), Triangle_Normal(pMesh, _pm, h));
    put(this->facet_normal,
        face(opposite(h, pMesh), pMesh),
        Triangle_Normal(pMesh, _pm, opposite(h, pMesh)));

    put(this->Facet_Flag, face(h, pMesh), TO_BE_REMOVED);
    // h->facet()->Patch_Index = Vertex_number;

    put(this->Facet_Flag, face(opposite(h, pMesh), pMesh), TO_BE_REMOVED);
    // h->opposite()->facet()->Patch_Index = Vertex_number;
  }

  else if((type == 14) || (type == 15))
  {
    if(get(this->Vertex_Sign, target(g, pMesh)) == NOSIGN)
      put(this->Vertex_Sign, target(g, pMesh), MINUS);
    if(get(this->Vertex_Sign, target(next(g, pMesh), pMesh)) == NOSIGN)
      put(this->Vertex_Sign, target(next(g, pMesh), pMesh), PLUS);
    if(get(this->Vertex_Sign, target(next(next(g, pMesh), pMesh), pMesh)) ==
       NOSIGN)
      put(this->Vertex_Sign, target(next(next(g, pMesh), pMesh), pMesh), MINUS);
    if(get(this->Vertex_Sign,
           target(next(next(next(g, pMesh), pMesh), pMesh), pMesh)) == NOSIGN)
      put(this->Vertex_Sign,
          target(next(next(next(g, pMesh), pMesh), pMesh), pMesh),
          PLUS);


    g = next(g, pMesh);
#ifdef DBG_Retriangulation
    std::cout << "split face" << std::endl;
#endif
    h = CGAL::Euler::split_face(h, g, pMesh);

    put(this->facet_Component_Number, face(h, pMesh), Component_ID);
    put(this->facet_Component_Number,
        face(opposite(h, pMesh), pMesh),
        Component_ID);

    put(this->Facet_Flag, face(opposite(h, pMesh), pMesh), TO_BE_REMOVED);
    // h->opposite()->facet()->Patch_Index = Vertex_number;

    put(this->facet_normal,
        face(opposite(h, pMesh), pMesh),
        Triangle_Normal(pMesh, _pm, opposite(h, pMesh)));

    g = next(next(h, pMesh), pMesh);
#ifdef DBG_Retriangulation
    std::cout << "split face" << std::endl;
#endif
    h = CGAL::Euler::split_face(h, g, pMesh);

    put(this->facet_Component_Number, face(h, pMesh), Component_ID);
    put(this->facet_Component_Number,
        face(opposite(h, pMesh), pMesh),
        Component_ID);

    put(this->Facet_Flag, face(opposite(h, pMesh), pMesh), TO_BE_REMOVED);
    // h->opposite()->facet()->Patch_Index = Vertex_number;

    put(this->facet_normal,
        face(opposite(h, pMesh), pMesh),
        Triangle_Normal(pMesh, _pm, opposite(h, pMesh)));

    g = next(next(h, pMesh), pMesh);
#ifdef DBG_Retriangulation
    std::cout << "split face" << std::endl;
#endif
    h = CGAL::Euler::split_face(h, g, pMesh);

    put(this->facet_Component_Number, face(h, pMesh), Component_ID);
    put(this->facet_Component_Number,
        face(opposite(h, pMesh), pMesh),
        Component_ID);

    put(this->facet_normal, face(h, pMesh), Triangle_Normal(pMesh, _pm, h));
    put(this->facet_normal,
        face(opposite(h, pMesh), pMesh),
        Triangle_Normal(pMesh, _pm, opposite(h, pMesh)));

    put(this->Facet_Flag, face(h, pMesh), TO_BE_REMOVED);
    // h->facet()->Patch_Index = Vertex_number;
    put(this->Facet_Flag, face(opposite(h, pMesh), pMesh), TO_BE_REMOVED);
    // h->opposite()->facet()->Patch_Index = Vertex_number;
  }
}


template< typename HalfedgeGraph, typename PointMap, typename VertexColorMap >
typename Compression_Valence_Component< HalfedgeGraph,
                                        PointMap,
                                        VertexColorMap >::Vector
Compression_Valence_Component< HalfedgeGraph, PointMap, VertexColorMap >::
    Normal_Patch(const halfedge_descriptor &const_h,
                 const unsigned int &valence,
                 const HalfedgeGraph &_pMesh,
                 const PointMap *_pm)
{
  FEVV::Geometry_traits< HalfedgeGraph > gt(_pMesh);

  halfedge_descriptor h = const_h;
  int type = Find_Type(_pMesh, h, valence);

  double area[5] = {0, 0, 0, 0, 0};

  Vector *normals = new Vector[5];
  Vector normal = Vector(0.0, 0.0, 0.0);
  for(int i = 0; i < 5; i++)
  {
    normals[i] = Vector(0.0, 0.0, 0.0);
  }

  // Triangle
  if((type == 1) || (type == 2) || (type == 4))
  {
    normals[1] = Triangle_Normal(_pMesh, _pm, h);
    area[1] = Area_Facet_Triangle(h, _pMesh, _pm);
  }
  else if(type == 3)
  {
    normals[1] = Triangle_Normal(_pMesh, _pm, h);
    area[1] = Area_Facet_Triangle(h, _pMesh, _pm);
  }

  // quadrangle
  else if((type == 5) || (type == 8))
  {
    normals[1] = Triangle_Normal(_pMesh, _pm, h);
    area[1] = Area_Facet_Triangle(h, _pMesh, _pm);
    h = opposite(prev(h, _pMesh), _pMesh);
    normals[2] = Triangle_Normal(_pMesh, _pm, h);
    area[2] = Area_Facet_Triangle(h, _pMesh, _pm);
  }

  else if((type == 6) || (type == 7))
  {
    normals[1] = Triangle_Normal(_pMesh, _pm, h);
    area[1] = Area_Facet_Triangle(h, _pMesh, _pm);

    h = opposite(next(h, _pMesh), _pMesh);
    normals[2] = Triangle_Normal(_pMesh, _pm, h);
    area[2] = Area_Facet_Triangle(h, _pMesh, _pm);
  }

  // pentagone
  else if((type == 9) || (type == 12))
  {
    normals[1] = Triangle_Normal(_pMesh, _pm, h);
    area[1] = Area_Facet_Triangle(h, _pMesh, _pm);

    h = opposite(prev(h, _pMesh), _pMesh);
    normals[2] = Triangle_Normal(_pMesh, _pm, h);
    area[2] = Area_Facet_Triangle(h, _pMesh, _pm);

    h = opposite(next(h, _pMesh), _pMesh);
    normals[3] = Triangle_Normal(_pMesh, _pm, h);
    area[3] = Area_Facet_Triangle(h, _pMesh, _pm);
  }


  else if(type == 10)
  {
    normals[1] = Triangle_Normal(_pMesh, _pm, h);
    area[1] = Area_Facet_Triangle(h, _pMesh, _pm);

    h = opposite(next(h, _pMesh), _pMesh);
    normals[2] = Triangle_Normal(_pMesh, _pm, h);
    area[2] = Area_Facet_Triangle(h, _pMesh, _pm);

    h = opposite(prev(h, _pMesh), _pMesh);
    normals[3] = Triangle_Normal(_pMesh, _pm, h);
    area[3] = Area_Facet_Triangle(h, _pMesh, _pm);
  }

  else if(type == 11)
  {
    halfedge_descriptor g = h;
    normals[1] = Triangle_Normal(_pMesh, _pm, h);
    area[1] = Area_Facet_Triangle(h, _pMesh, _pm);

    h = opposite(prev(g, _pMesh), _pMesh);
    normals[2] = Triangle_Normal(_pMesh, _pm, h);
    area[2] = Area_Facet_Triangle(h, _pMesh, _pm);

    h = opposite(next(g, _pMesh), _pMesh);
    normals[3] = Triangle_Normal(_pMesh, _pm, h);
    area[3] = Area_Facet_Triangle(h, _pMesh, _pm);
  }

  // Hexagone
  else if((type == 13) || (type == 16))
  {
    halfedge_descriptor g;
    normals[1] = Triangle_Normal(_pMesh, _pm, h);
    area[1] = Area_Facet_Triangle(h, _pMesh, _pm);

    h = opposite(prev(h, _pMesh), _pMesh);
    g = h;
    normals[2] = Triangle_Normal(_pMesh, _pm, h);
    area[2] = Area_Facet_Triangle(h, _pMesh, _pm);

    h = opposite(prev(g, _pMesh), _pMesh);
    normals[3] = Triangle_Normal(_pMesh, _pm, h);
    area[3] = Area_Facet_Triangle(h, _pMesh, _pm);

    h = opposite(next(g, _pMesh), _pMesh);
    normals[4] = Triangle_Normal(_pMesh, _pm, h);
    area[4] = Area_Facet_Triangle(h, _pMesh, _pm);
  }

  else if((type == 14) || (type == 15))
  {
    halfedge_descriptor g;
    normals[1] = Triangle_Normal(_pMesh, _pm, h);
    area[1] = Area_Facet_Triangle(h, _pMesh, _pm);

    h = opposite(next(h, _pMesh), _pMesh);
    g = h;
    normals[2] = Triangle_Normal(_pMesh, _pm, h);
    area[2] = Area_Facet_Triangle(h, _pMesh, _pm);

    h = opposite(prev(g, _pMesh), _pMesh);
    normals[3] = Triangle_Normal(_pMesh, _pm, h);
    area[3] = Area_Facet_Triangle(h, _pMesh, _pm);

    h = opposite(next(g, _pMesh), _pMesh);
    normals[4] = Triangle_Normal(_pMesh, _pm, h);
    area[4] = Area_Facet_Triangle(h, _pMesh, _pm);
  }

  for(unsigned int i = 0; i < (valence - 2); i++)
    area[0] = area[0] + area[i + 1];

  for(unsigned int i = 0; i < (valence - 2); i++)
    area[i + 1] = area[i + 1] / area[0];

  for(unsigned int i = 0; i < (valence - 2); i++)
    normal = normal + FEVV::Math::Vector::scalar_mult<
                          FEVV::Geometry_traits< HalfedgeGraph > >(
                          normals[i + 1], area[i + 1], gt);

  if(area[0] == 0.0)
    normal = Vector(0.0, 0.0, 0.0);

  double length =
      FEVV::Math::Vector::l2_distance< FEVV::Geometry_traits< HalfedgeGraph > >(
          normal, gt);
  if(length != 0)
    normal = FEVV::Math::Vector::scalar_mult<
        FEVV::Geometry_traits< HalfedgeGraph > >(normal, 1.0 / length, gt);

  delete[] normals;
  return normal;
}


template< typename HalfedgeGraph, typename PointMap, typename VertexColorMap >
typename Compression_Valence_Component< HalfedgeGraph,
                                        PointMap,
                                        VertexColorMap >::Vector
Compression_Valence_Component< HalfedgeGraph, PointMap, VertexColorMap >::
    Calculate_T1_T2(const halfedge_descriptor &h,
                    const Vector &normal,
                    Vector &T2,
                    const HalfedgeGraph &_pMesh,
                    const PointMap *_pm)
{
  FEVV::Geometry_traits< HalfedgeGraph > gt(_pMesh);

  Point3d P = get(*_pm, target(h, _pMesh));
  Point3d Q = get(*_pm, target(opposite(h, _pMesh), _pMesh));

  Vector T1 = Vector(0.0, 0.0, 0.0);
  Vector u = FEVV::Math::Vector::sub< FEVV::Geometry_traits< HalfedgeGraph > >(
      P, Q, gt);

  double length =
      FEVV::Math::Vector::l2_distance< FEVV::Geometry_traits< HalfedgeGraph > >(
          u, gt);
  if(length != 0)
    u = FEVV::Math::Vector::scalar_mult<
        FEVV::Geometry_traits< HalfedgeGraph > >(u, 1.0 / length, gt);

  double product =
      FEVV::Math::Vector::l2_distance< FEVV::Geometry_traits< HalfedgeGraph > >(
          u, gt) *
      FEVV::Math::Vector::l2_distance< FEVV::Geometry_traits< HalfedgeGraph > >(
          normal, gt);

  // cosine
  double dot =
      FEVV::Math::Vector::dot_product< FEVV::Geometry_traits< HalfedgeGraph > >(
          u, normal, gt);
  double cosine = 0;
  if(product != 0)
    cosine = dot / product;

  if(cosine > 1.0)
    cosine = 1.0;
  if(cosine < -1.0)
    cosine = -1.0;

  double cosine_rad = std::acos(cosine);

  double beta_rad = 0;
  if(cosine_rad <= PI / 2)
    beta_rad = PI / 2 - cosine_rad;
  else
    beta_rad = cosine_rad - PI / 2;
  double beta = std::cos(beta_rad);

  if(beta != 0)
  {
    Vector tmp1 = FEVV::Math::Vector::scalar_mult<
        FEVV::Geometry_traits< HalfedgeGraph > >(
        normal, cosine, gt); // = cosine * normal
    Vector tmp2 = u - tmp1;  // = u - cosine * normal
    T1 = FEVV::Math::Vector::scalar_mult<
        FEVV::Geometry_traits< HalfedgeGraph > >(
        tmp2, 1.0 / beta, gt); // = (u - cosine * normal) / beta
  }
  else
    T1 = Vector(0.0, 0.0, 0.0);

  T2 = FEVV::Math::Vector::cross_product<
      FEVV::Geometry_traits< HalfedgeGraph > >(normal, T1, gt);

  return T1;
}


static inline int
Signe(const double &x)
{
  return (x < 0.) ? -1 : 1;
}

template< typename T >
static inline T
square(const T x)
{
  return x * x;
}


template< typename HalfedgeGraph, typename PointMap, typename VertexColorMap >
Point_Int
Compression_Valence_Component< HalfedgeGraph, PointMap, VertexColorMap >::
    Frenet_Rotation(const Point_Int &Dist,
                    const Vector &T1,
                    const Vector &T2,
                    const Vector &normal)
{
  typedef Eigen::Matrix3d Matrix;

  Matrix R;
  R << T1[0], T2[0], normal[0], T1[1], T2[1], normal[1], T1[2], T2[2],
      normal[2];
  Matrix M = R;
  R.transposeInPlace();
  Vector Dist1(Dist.x, Dist.y, Dist.z);
  // Vector Dist2 = Rt * Dist1; // MT
  Matrix D1;
  D1 << 1, 0, 0, 0, 1, 0, 0, 0, 1;
  Matrix D2;
  D2 << 0, 1, 0, 1, 0, 0, 0, 0, 1;
  Matrix D3;
  D3 << -1, 0, 0, 0, -1, 0, 0, 0, 1;
  Matrix D4;
  D4 << 1, 0, 0, 0, 0, 1, 0, 1, 0;
  Matrix D5;
  D5 << 1, 0, 0, 0, -1, 0, 0, 0, -1;
  Matrix *S = new Matrix[16];

  // Verify in order to find the smallest rotation angle.
  if(std::abs(M(0, 2)) > std::abs(M(1, 2)))
    S[0] = D2;
  else
    S[0] = D1;
  M = S[0] * M;

  if(M(1, 2) < 0)
    S[1] = D3;
  else
    S[1] = D1;

  M = S[1] * M;

  /// first rotation angle : phi;
  double phi = -100;

  if(square(M(0, 2)) + square(M(1, 2)) == 0)
    phi = 0;
  else
    phi = Signe(-1 * M(0, 2)) *
          std::acos(M(1, 2) / std::sqrt(square(M(0, 2)) + square(M(1, 2))));

  Matrix R1;
  R1 << std::cos(phi), -std::sin(phi), 0, std::sin(phi), std::cos(phi), 0, 0, 0,
      1;

  S[2] << 1, -std::tan(phi / 2), 0, 0, 1, 0, 0, 0, 1;
  S[3] << 1, 0, 0, std::sin(phi), 1, 0, 0, 0, 1;
  S[4] = S[2];

  Matrix R1inv;
  R1inv << std::cos(phi), std::sin(phi), 0, -std::sin(phi), std::cos(phi), 0, 0,
      0, 1;

  M = R1inv * M;

  if(std::abs(M(1, 2)) > std::abs(M(2, 2)))
    S[5] = D4;
  else
    S[5] = D1;

  M = S[5] * M;

  if(M(2, 2) < 0)
    S[6] = D5;
  else
    S[6] = D1;

  M = S[6] * M;
  double psi = -100;

  /// Second rotation angle psi.
  if(square(M(1, 2)) + square(M(2, 2)) == 0)
    psi = 0;
  else
    psi = Signe(-1 * M(1, 2)) *
          std::acos(M(2, 2) / std::sqrt(square(M(1, 2)) + square(M(2, 2))));

  Matrix R2;
  R2 << 1, 0, 0, 0, std::cos(psi), -std::sin(psi), 0, std::sin(psi),
      std::cos(psi);
  S[7] << 1, 0, 0, 0, 1, -std::tan(psi / 2), 0, 0, 1;
  S[8] << 1, 0, 0, 0, 1, 0, 0, std::sin(psi), 1;
  S[9] = S[7];

  Matrix R2inv;
  R2inv << 1, 0, 0, 0, std::cos(psi), std::sin(psi), 0, -std::sin(psi),
      std::cos(psi);
  M = R2inv * M;

  if(std::abs(M(0, 1)) > std::abs(M(1, 1)))
    S[10] = D2;
  else
    S[10] = D1;
  M = S[10] * M;

  if(M(1, 1) < 0)
    S[11] = D3;
  else
    S[11] = D1;
  M = S[11] * M;
  double theta = -100;

  /// Last rotation angle theta.
  if(square(M(0, 1)) + square(M(1, 1)) == 0)
    theta = 0;
  else
    theta = Signe(-1 * M(0, 1)) *
            std::acos(M(1, 1) / std::sqrt(square(M(0, 1)) + square(M(1, 1))));

  S[12] << 1, -std::tan(theta / 2), 0, 0, 1, 0, 0, 0, 1;
  S[13] << 1, 0, 0, std::sin(theta), 1, 0, 0, 0, 1;
  S[14] = S[12];

  Matrix R3;
  R3 << std::cos(theta), -std::sin(theta), 0, std::sin(theta), std::cos(theta),
      0, 0, 0, 1;
  Matrix R3inv = R3.inverse();
  /*Matrix S16 = */ R3inv *M; // MT //ELO useless-line

  Eigen::Vector3d u;
  u << Dist.x, Dist.y, Dist.z;
  Matrix m_inter;

  // Procedure of the bijection.
  for(int i = 0; i < 15; i++)
  {
    if((i == 0) || (i == 1) || (i == 5) || (i == 6) || (i == 10) || (i == 11))
      m_inter = S[i];
    else
      m_inter << S[i](0, 0), -S[i](0, 1), -S[i](0, 2), -S[i](1, 0), S[i](1, 1),
          -S[i](1, 2), -S[i](2, 0), -S[i](2, 1), S[i](2, 2);
    u = m_inter * u;

    int x = 0, y = 0, z = 0;
    x = ceil(u(0) - 0.5);
    y = ceil(u(1) - 0.5);
    z = ceil(u(2) - 0.5);

    u << x, y, z;
  }

  Point_Int New_Coordinates;
  New_Coordinates.x = (int)u(0);
  New_Coordinates.y = (int)u(1);
  New_Coordinates.z = (int)u(2);

  delete[] S;
  return New_Coordinates;
}

//#define DBG_Inverse_Frenet_Rotation

template< typename HalfedgeGraph, typename PointMap, typename VertexColorMap >
Point_Int
Compression_Valence_Component< HalfedgeGraph, PointMap, VertexColorMap >::
    Inverse_Frenet_Rotation(const Point_Int &Frenet,
                            const Vector &T1,
                            const Vector &T2,
                            const Vector &normal)
{

#ifdef DBG_Inverse_Frenet_Rotation
  static unsigned int dbg_Inverse_Frenet_Rotation_call_cnt = 0;
  std::cout << __func__ << "  call #" << ++dbg_Inverse_Frenet_Rotation_call_cnt
            << std::endl;
#endif

  typedef Eigen::Matrix3d Matrix;

  Matrix R;
  R << T1[0], T2[0], normal[0], T1[1], T2[1], normal[1], T1[2], T2[2],
      normal[2];
  Matrix M = R;

  Matrix D1;
  D1 << 1, 0, 0, 0, 1, 0, 0, 0, 1;
  Matrix D2;
  D2 << 0, 1, 0, 1, 0, 0, 0, 0, 1;
  Matrix D3;
  D3 << -1, 0, 0, 0, -1, 0, 0, 0, 1;
  Matrix D4;
  D4 << 1, 0, 0, 0, 0, 1, 0, 1, 0;
  Matrix D5;
  D5 << 1, 0, 0, 0, -1, 0, 0, 0, -1;
  Matrix *S = new Matrix[16];

  // Verify in order to find the smallest rotation angle.
  if(std::abs(M(0, 2)) > std::abs(M(1, 2)))
    S[0] = D2;
  else
    S[0] = D1;
  M = S[0] * M;

  if(M(1, 2) < 0)
    S[1] = D3;
  else
    S[1] = D1;

  M = S[1] * M;

  /// first rotation angle : phi;
  double phi = -100;

  if(square(M(0, 2)) + square(M(1, 2)) == 0)
    phi = 0;
  else
    phi = Signe(-1 * M(0, 2)) *
          std::acos(M(1, 2) / std::sqrt(square(M(0, 2)) + square(M(1, 2))));

  Matrix R1;
  R1 << std::cos(phi), -std::sin(phi), 0, std::sin(phi), std::cos(phi), 0, 0, 0,
      1;

  S[2] << 1, -std::tan(phi / 2), 0, 0, 1, 0, 0, 0, 1;
  S[3] << 1, 0, 0, std::sin(phi), 1, 0, 0, 0, 1;
  S[4] = S[2];

  // Matrix R1inv = inverse_matrix(R1);
  Matrix R1inv;
  R1inv << std::cos(phi), std::sin(phi), 0, -std::sin(phi), std::cos(phi), 0, 0,
      0, 1;

  M = R1inv * M;

  if(std::abs(M(1, 2)) > std::abs(M(2, 2)))
    S[5] = D4;
  else
    S[5] = D1;

  M = S[5] * M;

  if(M(2, 2) < 0)
    S[6] = D5;
  else
    S[6] = D1;

  M = S[6] * M;
  double psi = -100;

  /// Second rotation angle psi.
  if(square(M(1, 2)) + square(M(2, 2)) == 0)
    psi = 0;
  else
    psi = Signe(-1 * M(1, 2)) *
          std::acos(M(2, 2) / std::sqrt(square(M(1, 2)) + square(M(2, 2))));

  Matrix R2;
  R2 << 1, 0, 0, 0, std::cos(psi), -std::sin(psi), 0, std::sin(psi),
      std::cos(psi);
  S[7] << 1, 0, 0, 0, 1, -std::tan(psi / 2), 0, 0, 1;
  S[8] << 1, 0, 0, 0, 1, 0, 0, std::sin(psi), 1;
  S[9] = S[7];

  // Matrix R2inv = inverse_matrix(R2);
  Matrix R2inv;
  R2inv << 1, 0, 0, 0, std::cos(psi), std::sin(psi), 0, -std::sin(psi),
      std::cos(psi);
  M = R2inv * M;

  if(std::abs(M(0, 1)) > std::abs(M(1, 1)))
    S[10] = D2;
  else
    S[10] = D1;
  M = S[10] * M;

  if(M(1, 1) < 0)
    S[11] = D3;
  else
    S[11] = D1;
  M = S[11] * M;
  double theta = -100;

  /// Last rotation angle theta.
  if(square(M(0, 1)) + square(M(1, 1)) == 0)
    theta = 0;
  else
    theta = Signe(-1 * M(0, 1)) *
            std::acos(M(1, 1) / std::sqrt(square(M(0, 1)) + square(M(1, 1))));

  S[12] << 1, -std::tan(theta / 2), 0, 0, 1, 0, 0, 0, 1;
  S[13] << 1, 0, 0, std::sin(theta), 1, 0, 0, 0, 1;
  S[14] = S[12];

  Matrix R3;
  R3 << std::cos(theta), -std::sin(theta), 0, std::sin(theta), std::cos(theta),
      0, 0, 0, 1;
  Matrix R3inv = R3.inverse();
  /*Matrix S16 = */ R3inv *M; // MT //ELO  useless-line

  Eigen::Vector3d u;
  u << Frenet.x, Frenet.y, Frenet.z;
  Matrix m_inter;

  for(int i = 14; i > -1; i--)
  {
    m_inter = S[i];

#ifdef DBG_Inverse_Frenet_Rotation
    if(dbg_Inverse_Frenet_Rotation_call_cnt == 1596)
    {
      std::cout << __func__ << "  mark #1"
                << "  i=" << i << "  u=" << u(0) << " " << u(1) << " " << u(2)
                << "  m_inter=" << m_inter(0, 0) << " " << m_inter(1, 0) << " "
                << m_inter(2, 0) << " " << m_inter(0, 1) << " " << m_inter(1, 1)
                << " " << m_inter(2, 1) << " " << m_inter(0, 2) << " "
                << m_inter(1, 2) << " " << m_inter(2, 2) << std::endl;
    }
#endif

    // ELO-porting-note:
    //  u =  -1 * m_inter * u;
    //  calls
    //  1) MatrixC33 operator* ( RT const& c, MatrixC33 const& m )
    //     at
    //     Mepp1.github/src/components/Compression/Compression_Valence/src/Matrix3X3.h:253
    //     with c = -1
    //  2) Vector_3 operator* ( MatrixC33 const& m, Vector_3 const& v )
    //     at
    //     Mepp1.github/src/components/Compression/Compression_Valence/src/Matrix3X3.h:263
    u = -1 * m_inter * u;

    int x = 0, y = 0, z = 0;
    x = -std::ceil(u(0) - 0.5);
    y = -std::ceil(u(1) - 0.5);
    z = -std::ceil(u(2) - 0.5);

#ifdef DBG_Inverse_Frenet_Rotation
    if(dbg_Inverse_Frenet_Rotation_call_cnt == 1596)
    {
      std::cout << __func__ << "  mark #2"
                << "  i=" << i << "  u=" << u(0) << " " << u(1) << " " << u(2)
                << "  x=" << x << "  y=" << y << "  z=" << z << std::endl;
    }
#endif

    u << x, y, z;
#ifdef DBG_Inverse_Frenet_Rotation
    if(dbg_Inverse_Frenet_Rotation_call_cnt == 1596)
    {
      std::cout << __func__ << "  mark #3"
                << "  i=" << i << "  u=" << u(0) << " " << u(1) << " " << u(2)
                << std::endl;
    }
#endif
  }

  Point_Int Dist;
  Dist.x = (int)u(0);
  Dist.y = (int)u(1);
  Dist.z = (int)u(2);

#ifdef DBG_Inverse_Frenet_Rotation
  if(dbg_Inverse_Frenet_Rotation_call_cnt == 1596)
  {
    std::cout << __func__ << "  mark #4"
              << "  Dist=" << Dist.x << " " << Dist.y << " " << Dist.z
              << std::endl;
  }
#endif

  delete[] S;
  return Dist;
}


template< typename HalfedgeGraph, typename PointMap, typename VertexColorMap >
Color_Unit
Compression_Valence_Component< HalfedgeGraph, PointMap, VertexColorMap >::
    Get_Average_Vertex_Color_After_Removal(const halfedge_descriptor &h,
                                           const int &valence,
                                           const HalfedgeGraph &_pMesh)
{
  halfedge_descriptor g = h;

  Color_Unit Average_color;
  Average_color.c0 = 0;
  Average_color.c1 = 0;
  Average_color.c2 = 0;

  for(int i = 0; i < valence; i++)
  {
    Color_Unit col = Get_Vertex_Color(g, _pMesh);
    Average_color.c0 += col.c0;
    Average_color.c1 += col.c1;
    Average_color.c2 += col.c2;

    g = next(g, _pMesh);
  }

  Color_Unit Resulting_color;
  Resulting_color.c0 = floor((float)Average_color.c0 / valence + 0.5);
  Resulting_color.c1 = floor((float)Average_color.c1 / valence + 0.5);
  Resulting_color.c2 = floor((float)Average_color.c2 / valence + 0.5);

  return Resulting_color;
}


template< typename HalfedgeGraph, typename PointMap, typename VertexColorMap >
int
Compression_Valence_Component< HalfedgeGraph, PointMap, VertexColorMap >::
    Regulation(HalfedgeGraph &_pMesh,
               const bool Normal_flipping,
               const bool Use_metric,
               const float &Metric_thread,
               const bool Use_forget_metric,
               const int &Forget_value,
               const int &Component_ID,
               const PointMap *_pm)
{
  double Max_color, Mean_color;
  int Temp_NV = 0;
  int Number_facets;
  this->Calculate_Edge_Color_Difference(
      _pMesh, Component_ID, Max_color, Mean_color, Temp_NV);
  this->Recalculate_Component_Area(_pMesh, _pm, Component_ID, Number_facets);
  double Mean_area =
      (double)this->ComponentArea[Component_ID] / (double)Number_facets;

  // Initialization
  Init(_pMesh);


  // Number of removed vertices and number of connectivity symbols
  int Number_vertices = 0;
  int Number_symbol = 0;

  auto halfedge_iterator_pair = halfedges(_pMesh);
  halfedge_iterator hi = halfedge_iterator_pair.first;

  while(
      (get(this->vertex_Seed_Edge, target(*hi, _pMesh)) != 2 * Component_ID) ||
      (get(this->vertex_Seed_Edge, target(opposite(*hi, _pMesh), _pMesh)) !=
       2 * Component_ID + 1))
    hi++;

  halfedge_descriptor First_halfedge = *hi;

  put(this->Vertex_Flag, target(First_halfedge, _pMesh), CONQUERED);
  put(this->Vertex_Flag,
      target(opposite(First_halfedge, _pMesh), _pMesh),
      CONQUERED);

  std::queue< halfedge_descriptor > Halfedges;
  Halfedges.push(First_halfedge);

  halfedge_descriptor h;

  while(!Halfedges.empty())
  {
    h = Halfedges.front();
    Halfedges.pop();

    size_t valence = out_degree(target(next(h, _pMesh), _pMesh), _pMesh);
    // int Component_ID = h->next()->vertex()->Component_Number;

    if((get(this->Facet_Flag, face(h, _pMesh)) == CONQUERED) ||
       (get(this->Facet_Flag, face(h, _pMesh)) == TO_BE_REMOVED))
      continue;

    else if((get(this->Vertex_Flag, target(next(h, _pMesh), _pMesh)) == FREE) &&
            (valence == 3) &&
            (Is_Border_Vertex(next(h, _pMesh), _pMesh) ==
             false)) // if valence is 3, remove the front vertex.
    {
      halfedge_descriptor g = h;
      int type = 1; // ant type of valence 3

      // Check if the manifold property is violated.
      // bool Manifold_condition = Is_Manifold_Property_Violated(h, type,
      // valence);

      // calculate error caused by the removal. This metric decides if the
      // vertex can be removed or not.
      bool Geometric_metric_condition = false;
      if(Use_metric == true)
      {
        if(Use_forget_metric == true)
        {
          if(FEVV::size_of_vertices(_pMesh) > (unsigned)Forget_value)
            Geometric_metric_condition = false;
          else
            Geometric_metric_condition = Is_Geometric_Metric_Violated(
                _pMesh, _pm, h, type, valence, Metric_thread);
        }

        else
          Geometric_metric_condition = Is_Geometric_Metric_Violated(
              _pMesh, _pm, h, type, valence, Metric_thread);
      }
      bool Is_Color_Too_Important = false;

#ifdef USE_COLOR_METRIC
      if(this->IsColored)
      {
        Is_Color_Too_Important = this->Error_Projected_Surface(
            _pMesh, _pm, h, Component_ID, Mean_color, Mean_area);
      }
#endif

      // remove the front vertex if its removal does not viloate the manifold
      // property and some metrics
      bool Check_all_condition = false;

      if((!Geometric_metric_condition) && (!Is_Color_Too_Important))
        Check_all_condition = true;

      // All conditions are goods. Remove the front vertex
      if(Check_all_condition)
      {
        Number_vertices++;
        Number_symbol++;

        // Insertion of a symbol 3 to the list
        this->InterConnectivity.push_front(0);

        Point3d Geo_info = get(*_pm, target(next(h, _pMesh), _pMesh));
        Point_Int Geo = Change_Real_Int(Geo_info, Component_ID);

        Point3d Barycenter = Barycenter_Patch_Before_Removal(g, _pMesh, _pm);
        Point_Int BC = Change_Real_Int(Barycenter, Component_ID);


        if((this->IsColored) && (!this->IsOneColor))
        {
          Color_Unit Removed_vertex_color;

          Removed_vertex_color = Get_Vertex_Color(next(h, _pMesh), _pMesh);

          Color_Unit Average_color;
          Average_color = Get_Average_Vertex_Color_Lee(g, valence, _pMesh);

          Color_Unit Color_diff = Removed_vertex_color - Average_color;

          this->InterVertexColor.push_front(Color_diff);
        }

        g = h;

        // P, Q, R are used to compute the normal of triangle after removal
        Point3d P = get(*_pm, target(g, _pMesh));
        g = next(opposite(next(g, _pMesh), _pMesh), _pMesh);
        Point3d Q = get(*_pm, target(g, _pMesh));
        g = next(opposite(next(g, _pMesh), _pMesh), _pMesh);
        Point3d R = get(*_pm, target(g, _pMesh));
        g = h;

        Vector normal = Triangle_Normal(_pMesh, P, Q, R);
        Vector T2 = Vector(0.0, 0.0, 0.0);
        Vector T1 = Calculate_T1_T2(h, normal, T2, _pMesh, _pm);

        if(T1 == Vector(0.0, 0.0, 0.0))
        {
          T1 = Vector(1, 0, 0);
          T2 = Vector(0, 1, 0);
          normal = Vector(0, 0, 1);
        }

        else if(normal == Vector(0.0, 0.0, 0.0))
        {
          T1 = Vector(1, 0, 0);
          T2 = Vector(0, 1, 0);
          normal = Vector(0, 0, 1);
        }
        else if(T2 == Vector(0.0, 0.0, 0.0))
        {
          T1 = Vector(1, 0, 0);
          T2 = Vector(0, 1, 0);
          normal = Vector(0, 0, 1);
        }

        Point_Int Dist = Geo - BC;
        Point_Int Frenet_Coordinates;

        if(this->Is_Bijection_Enabled)
          Frenet_Coordinates = Frenet_Rotation(Dist, T1, T2, normal);
        else
          Frenet_Coordinates = Dist;

        this->InterGeometry.push_front(Frenet_Coordinates);

        g = next(h, _pMesh);
        put(this->Vertex_Flag, target(g, _pMesh), TO_BE_REMOVED);
        put(this->Facet_Flag, face(g, _pMesh), TO_BE_REMOVED);

        g = prev(opposite(g, _pMesh), _pMesh);
        put(this->Vertex_Flag, target(opposite(g, _pMesh), _pMesh), CONQUERED);
        put(this->Facet_Flag, face(g, _pMesh), TO_BE_REMOVED);

        if(!CGAL::is_border_edge(prev(g, _pMesh), _pMesh))
        {
          halfedge_descriptor h1 = opposite(prev(g, _pMesh), _pMesh);
          put(this->Facet_Flag, face(h1, _pMesh), CONQUERED);
          put(this->Vertex_Flag, target(next(h1, _pMesh), _pMesh), CONQUERED);
          if(!CGAL::is_border_edge(next(h1, _pMesh), _pMesh))
            Halfedges.push(opposite(next(h1, _pMesh), _pMesh));
          if(!CGAL::is_border_edge(prev(h1, _pMesh), _pMesh))
            Halfedges.push(opposite(prev(h1, _pMesh), _pMesh));
        }
        g = prev(opposite(g, _pMesh), _pMesh);
        put(this->Vertex_Flag, target(opposite(g, _pMesh), _pMesh), CONQUERED);
        put(this->Facet_Flag, face(g, _pMesh), TO_BE_REMOVED);
        if(!CGAL::is_border_edge(prev(g, _pMesh), _pMesh))
        {
          halfedge_descriptor h2 = opposite(prev(g, _pMesh), _pMesh);
          put(this->Facet_Flag, face(h2, _pMesh), CONQUERED);
          put(this->Vertex_Flag, target(next(h2, _pMesh), _pMesh), CONQUERED);
          if(!CGAL::is_border_edge(next(h2, _pMesh), _pMesh))
            Halfedges.push(opposite(next(h2, _pMesh), _pMesh));
          if(!CGAL::is_border_edge(prev(h2, _pMesh), _pMesh))
            Halfedges.push(opposite(prev(h2, _pMesh), _pMesh));
        }
      }

      else
      {
        this->InterConnectivity.push_front(1);
        Number_symbol++;
        put(this->Facet_Flag, face(h, _pMesh), CONQUERED);
        put(this->Vertex_Flag, target(next(h, _pMesh), _pMesh), CONQUERED);
        if(!CGAL::is_border_edge(next(h, _pMesh), _pMesh))
          Halfedges.push(opposite(next(h, _pMesh), _pMesh));
        if(!CGAL::is_border_edge(prev(h, _pMesh), _pMesh))
          Halfedges.push(opposite(prev(h, _pMesh), _pMesh));
      }
    }
    else // NULL triangle
    {
      this->InterConnectivity.push_front(1);
      Number_symbol++;
      put(this->Facet_Flag, face(h, _pMesh), CONQUERED);
      put(this->Vertex_Flag, target(next(h, _pMesh), _pMesh), CONQUERED);
      if(!CGAL::is_border_edge(next(h, _pMesh), _pMesh))
        Halfedges.push(opposite(next(h, _pMesh), _pMesh));
      if(!CGAL::is_border_edge(prev(h, _pMesh), _pMesh))
        Halfedges.push(opposite(prev(h, _pMesh), _pMesh));
    }
  }


  // Removal of all vertices with TO_BE_REMOVED flag
  vertex_iterator pVertex;
  for(pVertex = vertices(_pMesh).first; pVertex != vertices(_pMesh).second;)
  {
    vertex_descriptor vh = *pVertex;
    pVertex++;

    if(get(this->Vertex_Flag, vh) == TO_BE_REMOVED)
    {
      halfedge_descriptor temp = halfedge(vh, _pMesh);
      temp = CGAL::Euler::remove_center_vertex(temp, _pMesh);
      put(this->facet_Component_Number, face(temp, _pMesh), Component_ID);
    }
  }

  if((this->IsColored) && (!this->IsOneColor))
  {
    while(!this->InterVertexColor.empty())
    {
      Color_Unit Col = this->InterVertexColor.front();
      this->InterVertexColor.pop_front();
      this->VertexColor[Component_ID].push_front(Col);
    }
  }
  while(!InterConnectivity.empty())
  {
    int Symbol = this->InterConnectivity.front();
    this->InterConnectivity.pop_front();
    this->Connectivity[Component_ID].push_front(Symbol);
  }

  while(!InterGeometry.empty())
  {
    Point_Int Geo = InterGeometry.front();
    InterGeometry.pop_front();
    this->Geometry[Component_ID].push_front(Geo);
  }

  this->NumberSymbol[Component_ID].push_front(Number_symbol);
  this->NumberVertices[Component_ID].push_front(Number_vertices);

  this->DumpSymbolRegulation = Number_symbol;


  return Number_vertices;
}


// function imported from
// Mepp1.github/src/mepp/Polyhedron/polyhedron_shared_items.h
template< typename HalfedgeGraph, typename PointMap, typename VertexColorMap >
void
Compression_Valence_Component< HalfedgeGraph, PointMap, VertexColorMap >::
    compute_normals(const HalfedgeGraph &_pMesh, const PointMap *_pm)
{
  compute_normals_per_facet(_pMesh, _pm);
  compute_normals_per_vertex(_pMesh);
#if 0
	compute_normals_per_halfedge();  // ajout Céline
#endif
}


// normals per facet
// function imported from
// Mepp1.github/src/mepp/Polyhedron/polyhedron_shared_items.h
template< typename HalfedgeGraph, typename PointMap, typename VertexColorMap >
void
Compression_Valence_Component< HalfedgeGraph, PointMap, VertexColorMap >::
    compute_normals_per_facet(const HalfedgeGraph &_pMesh, const PointMap *_pm)
{
  face_iterator fi_begin = faces(_pMesh).first;
  face_iterator fi_end = faces(_pMesh).second;
  for(face_iterator fi = fi_begin; fi != fi_end; ++fi)
    compute_facet_normal(*fi, _pMesh, _pm);
}


// normals per vertex
// function imported from
// Mepp1.github/src/mepp/Polyhedron/polyhedron_shared_items.h
template< typename HalfedgeGraph, typename PointMap, typename VertexColorMap >
void
Compression_Valence_Component< HalfedgeGraph, PointMap, VertexColorMap >::
    compute_normals_per_vertex(const HalfedgeGraph &_pMesh)
{
  vertex_iterator vi_begin = vertices(_pMesh).first;
  vertex_iterator vi_end = vertices(_pMesh).second;
  for(vertex_iterator vi = vi_begin; vi != vi_end; ++vi)
    compute_vertex_normal(*vi, _pMesh);
}


// function imported from
// Mepp1.github/src/mepp/Polyhedron/polyhedron_shared_items.h
template< typename HalfedgeGraph, typename PointMap, typename VertexColorMap >
void
Compression_Valence_Component< HalfedgeGraph, PointMap, VertexColorMap >::
    compute_facet_normal(const face_descriptor &f,
                         const HalfedgeGraph &_pMesh,
                         const PointMap *_pm)
{
  FEVV::Geometry_traits< HalfedgeGraph > gt(_pMesh);

  Vector sum = Vector(0.0, 0.0, 0.0);

  CGAL::Halfedge_around_face_circulator< HalfedgeGraph > facet_begin(
      halfedge(f, _pMesh), _pMesh);
  CGAL::Halfedge_around_face_circulator< HalfedgeGraph > h = facet_begin;

  do
  {
    Point3d p1 = get(*_pm, target(*h, _pMesh));
    Point3d p2 = get(*_pm, target(next(*h, _pMesh), _pMesh));
    Point3d p3 = get(*_pm, target(next(next(*h, _pMesh), _pMesh), _pMesh));

    Vector v1 =
        FEVV::Math::Vector::sub< FEVV::Geometry_traits< HalfedgeGraph > >(
            p2, p1, gt);
    Vector v2 =
        FEVV::Math::Vector::sub< FEVV::Geometry_traits< HalfedgeGraph > >(
            p3, p2, gt);

    Vector normal = FEVV::Math::Vector::cross_product<
        FEVV::Geometry_traits< HalfedgeGraph > >(v1, v2, gt);

    double sqnorm = FEVV::Math::Vector::dot_product<
        FEVV::Geometry_traits< HalfedgeGraph > >(normal, normal, gt);
    if(sqnorm != 0)
    {
      normal = FEVV::Math::Vector::scalar_mult<
          FEVV::Geometry_traits< HalfedgeGraph > >(
          normal, 1.0 / std::sqrt(sqnorm), gt);
    }
    sum = sum + normal;
  } while(++h != facet_begin);

  double sqnorm =
      FEVV::Math::Vector::dot_product< FEVV::Geometry_traits< HalfedgeGraph > >(
          sum, sum, gt);
  if(sqnorm != 0.0)
  {
    put(this->facet_normal,
        f,
        FEVV::Math::Vector::scalar_mult<
            FEVV::Geometry_traits< HalfedgeGraph > >(
            sum, 1.0 / std::sqrt(sqnorm), gt));
  }
  else
  {
    put(this->facet_normal, f, Vector(0.0, 0.0, 0.0));
    // TRACE("degenerate face\n");
  }
};


// function imported from
// Mepp1.github/src/mepp/Polyhedron/polyhedron_shared_items.h
template< typename HalfedgeGraph, typename PointMap, typename VertexColorMap >
void
Compression_Valence_Component< HalfedgeGraph, PointMap, VertexColorMap >::
    compute_vertex_normal(const vertex_descriptor &v,
                          const HalfedgeGraph &_pMesh)
{
  FEVV::Geometry_traits< HalfedgeGraph > gt(_pMesh);

  Vector normal = Vector(0.0, 0.0, 0.0);

  CGAL::Halfedge_around_target_circulator< HalfedgeGraph > pHalfedge(v, _pMesh);
  CGAL::Halfedge_around_target_circulator< HalfedgeGraph > begin = pHalfedge;

  CGAL_For_all(pHalfedge, begin)
  {
    if(!CGAL::is_border_edge(*pHalfedge, _pMesh))
    {
      normal = normal + get(this->facet_normal, face(*pHalfedge, _pMesh));
    }
  }

  double sqnorm =
      FEVV::Math::Vector::dot_product< FEVV::Geometry_traits< HalfedgeGraph > >(
          normal, normal, gt);
  if(sqnorm != 0.0f)
  {
    put(this->vertex_normal,
        v,
        FEVV::Math::Vector::scalar_mult<
            FEVV::Geometry_traits< HalfedgeGraph > >(
            normal, 1.0 / std::sqrt(sqnorm), gt));
  }
  else
  {
    put(this->vertex_normal, v, Vector(0.0, 0.0, 0.0));
  }
}


template< typename HalfedgeGraph, typename PointMap, typename VertexColorMap >
int
Compression_Valence_Component< HalfedgeGraph, PointMap, VertexColorMap >::
    Estimate_Geometry_Quantization(double volume,
                                   double area,
                                   int number_vertices)
{
  double C = (double)volume / (double)area / number_vertices;

  double a = -1.248;
  double b = -0.954;

  int Q = floor(a * log(C) + b + 0.5);

  return Q;
}


template< typename HalfedgeGraph, typename PointMap, typename VertexColorMap >
int
Compression_Valence_Component< HalfedgeGraph, PointMap, VertexColorMap >::
    Get_Correct_Vector(int i, int j, int k)
{
  int Correct_symbol = -1;

  if((i == 0) && (j == 0) && (k == 0))
    Correct_symbol = 0;

  else if((i == 1) && (j == 0) && (k == 0))
    Correct_symbol = 1;

  else if((i == 0) && (j == 1) && (k == 0))
    Correct_symbol = 2;

  else if((i == 1) && (j == 1) && (k == 0))
    Correct_symbol = 3;

  else if((i == 0) && (j == 0) && (k == 1))
    Correct_symbol = 4;

  else if((i == 1) && (j == 0) && (k == 1))
    Correct_symbol = 5;

  else if((i == 0) && (j == 1) && (k == 1))
    Correct_symbol = 6;
  else
    Correct_symbol = 7;

  return Correct_symbol;
}


// ported from
//   Mepp1 Compression_Valence/src/Compression_Valence_Basemesh_builder.h
template< typename HalfedgeGraph, typename PointMap, typename VertexColorMap >
void
Compression_Valence_Component< HalfedgeGraph, PointMap, VertexColorMap >::
    build_mesh(HalfedgeGraph &_pMesh,
               PointMap *_pm,
               VertexColorMap *_v_cm,
               const std::vector< Point3d > &vlist, // points list
               const std::vector< int > &flist,     // faces list
               const std::vector< float > &clist,   // vertex-color list
               const std::vector< int > &Color_index_list)
{
  assert(vlist.size() != 0);
  assert(flist.size() != 0);

  typedef typename boost::property_traits< VertexColorMap >::value_type Color;

  // add vertices

  size_t Number_vertex = vlist.size();
  std::vector< vertex_descriptor > vertex_index_to_descriptor;
  vertex_index_to_descriptor.reserve(Number_vertex);
  for(size_t i = 0; i < Number_vertex; i++)
  {
    // create vertex
    vertex_descriptor v = add_vertex(_pMesh);
    vertex_index_to_descriptor.push_back(v);

    // set vertex geometry
    const Point3d &Pt = vlist[i];
    put(*_pm, v, Pt);

    // set vertex color
    if(!clist.empty())
    {
      Color tmp_color(clist[3 * i + 0], clist[3 * i + 1], clist[3 * i + 2]);
      put(*_v_cm, v, tmp_color);
    }
    if(!Color_index_list.empty())
    {
      // TODO-elo:
      // Vertex_Color_Index seems to be unused in the rest of the code ;
      // remove the Color_index_list parameter and this block when confirmed
      // not ported:
      // v->Vertex_Color_Index = (*Color_index_list)[i];
      throw std::runtime_error("Compression_Valence_Component::build_mesh(): "
                               "Vertex Color Index not supported, fix it!");
    }

    init_vertex_attributes(v);
  }

  // add faces

  size_t Number_facet = flist.size() / 3;

  for(size_t i = 0; i < Number_facet; i++)
  {
    // ELO  flist contains indices of face-vertices, 3 vertices per face
    std::vector< vertex_descriptor > face_vertices;
    face_vertices.push_back(vertex_index_to_descriptor[flist[3 * i + 0]]);
    face_vertices.push_back(vertex_index_to_descriptor[flist[3 * i + 1]]);
    face_vertices.push_back(vertex_index_to_descriptor[flist[3 * i + 2]]);

    face_descriptor currentFace = CGAL::Euler::add_face(face_vertices, _pMesh);
    if(currentFace == GraphTraits::null_face())
    {
      throw std::runtime_error("Compression_Valence_Component::build_mesh(): "
                               "failed to create face!");
    }
  }
}


template< typename HalfedgeGraph, typename PointMap, typename VertexColorMap >
bool
Compression_Valence_Component< HalfedgeGraph, PointMap, VertexColorMap >::
    Remove_Edges(HalfedgeGraph &_pMesh,
                 const halfedge_descriptor &h,
                 const int &type)
{
  bool check = false;
  halfedge_descriptor g = h;

  // triangle

  if((type == 1) || (type == 2) || (type == 4))
  {
    if(get(this->Vertex_Sign, target(next(g, _pMesh), _pMesh)) == NOSIGN)
      put(this->Vertex_Sign, target(next(g, _pMesh), _pMesh), PLUS);
  }

  else if(type == 3)
  {
    if(get(this->Vertex_Sign, target(next(g, _pMesh), _pMesh)) == NOSIGN)
      put(this->Vertex_Sign, target(next(g, _pMesh), _pMesh), MINUS);
  }


  // quadrangle
  else if((type == 5) || (type == 8))
  {
    // verification
    if(get(this->Facet_Flag, face(opposite(prev(g, _pMesh), _pMesh), _pMesh)) !=
       FREE)
      check = true;
    if(check == false)
    {
      g = prev(g, _pMesh);
      CGAL::Euler::join_face(g, _pMesh);

      g = h;
      if(get(this->Vertex_Sign, target(next(g, _pMesh), _pMesh)) == NOSIGN)
        put(this->Vertex_Sign, target(next(g, _pMesh), _pMesh), PLUS);
      if(get(this->Vertex_Sign,
             target(next(next(g, _pMesh), _pMesh), _pMesh)) == NOSIGN)
        put(this->Vertex_Sign,
            target(next(next(g, _pMesh), _pMesh), _pMesh),
            MINUS);
    }
  }

  else if((type == 6) || (type == 7))
  {
    // verification
    if(get(this->Facet_Flag, face(opposite(next(g, _pMesh), _pMesh), _pMesh)) !=
       FREE)
      check = true;
    if(check == false)
    {
      g = next(g, _pMesh);
      CGAL::Euler::join_face(g, _pMesh);

      g = h;
      if(get(this->Vertex_Sign, target(next(g, _pMesh), _pMesh)) == NOSIGN)
        put(this->Vertex_Sign, target(next(g, _pMesh), _pMesh), MINUS);
      if(get(this->Vertex_Sign,
             target(next(next(g, _pMesh), _pMesh), _pMesh)) == NOSIGN)
        put(this->Vertex_Sign,
            target(next(next(g, _pMesh), _pMesh), _pMesh),
            PLUS);
    }
  }

  // pentagone
  else if((type == 9) || (type == 12))
  {
    g = opposite(prev(g, _pMesh), _pMesh);
    if(get(this->Facet_Flag, face(g, _pMesh)) != FREE)
      check = true;
    g = opposite(next(g, _pMesh), _pMesh);
    if(get(this->Facet_Flag, face(g, _pMesh)) != FREE)
      check = true;

    if(check == false)
    {
      g = prev(h, _pMesh);
      g = CGAL::Euler::join_face(g, _pMesh);
      g = next(g, _pMesh);
      g = CGAL::Euler::join_face(g, _pMesh);

      g = h;
      if(get(this->Vertex_Sign, target(next(g, _pMesh), _pMesh)) == NOSIGN)
        put(this->Vertex_Sign, target(next(g, _pMesh), _pMesh), PLUS);
      if(get(this->Vertex_Sign,
             target(next(next(g, _pMesh), _pMesh), _pMesh)) == NOSIGN)
        put(this->Vertex_Sign,
            target(next(next(g, _pMesh), _pMesh), _pMesh),
            MINUS);
      if(get(this->Vertex_Sign,
             target(next(next(next(g, _pMesh), _pMesh), _pMesh), _pMesh)) ==
         NOSIGN)
        put(this->Vertex_Sign,
            target(next(next(next(g, _pMesh), _pMesh), _pMesh), _pMesh),
            PLUS);
    }
  }

  else if(type == 10)
  {
    g = opposite(next(g, _pMesh), _pMesh);
    if(get(this->Facet_Flag, face(g, _pMesh)) != FREE)
      check = true;

    g = opposite(prev(g, _pMesh), _pMesh);
    if(get(this->Facet_Flag, face(g, _pMesh)) != FREE)
      check = true;

    if(check == false)
    {
      g = opposite(next(h, _pMesh), _pMesh);
      g = CGAL::Euler::join_face(g, _pMesh);
      g = CGAL::Euler::join_face(g, _pMesh);


      g = h;
      if(get(this->Vertex_Sign, target(next(g, _pMesh), _pMesh)) == NOSIGN)
        put(this->Vertex_Sign, target(next(g, _pMesh), _pMesh), PLUS);
      if(get(this->Vertex_Sign,
             target(next(next(g, _pMesh), _pMesh), _pMesh)) == NOSIGN)
        put(this->Vertex_Sign,
            target(next(next(g, _pMesh), _pMesh), _pMesh),
            MINUS);
      if(get(this->Vertex_Sign,
             target(next(next(next(g, _pMesh), _pMesh), _pMesh), _pMesh)) ==
         NOSIGN)
        put(this->Vertex_Sign,
            target(next(next(next(g, _pMesh), _pMesh), _pMesh), _pMesh),
            PLUS);
    }
  }

  else if(type == 11)
  {
    if(get(this->Facet_Flag, face(opposite(next(g, _pMesh), _pMesh), _pMesh)) !=
       FREE)
      check = true;
    if(get(this->Facet_Flag, face(opposite(prev(g, _pMesh), _pMesh), _pMesh)) !=
       FREE)
      check = true;

    if(check == false)
    {
      g = next(g, _pMesh);
      g = CGAL::Euler::join_face(g, _pMesh);
      g = prev(g, _pMesh);
      g = CGAL::Euler::join_face(g, _pMesh);

      g = h;

      if(get(this->Vertex_Sign, target(next(g, _pMesh), _pMesh)) == NOSIGN)
        put(this->Vertex_Sign, target(next(g, _pMesh), _pMesh), MINUS);
      if(get(this->Vertex_Sign,
             target(next(next(g, _pMesh), _pMesh), _pMesh)) == NOSIGN)
        put(this->Vertex_Sign,
            target(next(next(g, _pMesh), _pMesh), _pMesh),
            PLUS);
      if(get(this->Vertex_Sign,
             target(next(next(next(g, _pMesh), _pMesh), _pMesh), _pMesh)) ==
         NOSIGN)
        put(this->Vertex_Sign,
            target(next(next(next(g, _pMesh), _pMesh), _pMesh), _pMesh),
            MINUS);
    }
  }

  else if((type == 13) || (type == 16))
  {
    g = opposite(prev(g, _pMesh), _pMesh);
    if(get(this->Facet_Flag, face(g, _pMesh)) != FREE)
      check = true;
    if(get(this->Facet_Flag, face(opposite(next(g, _pMesh), _pMesh), _pMesh)) !=
       FREE)
      check = true;
    if(get(this->Facet_Flag, face(opposite(prev(g, _pMesh), _pMesh), _pMesh)) !=
       FREE)
      check = true;

    if(check == false)
    {
      g = opposite(prev(h, _pMesh), _pMesh);
      g = CGAL::Euler::join_face(g, _pMesh);
      g = CGAL::Euler::join_face(g, _pMesh);
      g = CGAL::Euler::join_face(g, _pMesh);

      g = h;
      if(get(this->Vertex_Sign, target(next(g, _pMesh), _pMesh)) == NOSIGN)
        put(this->Vertex_Sign, target(next(g, _pMesh), _pMesh), PLUS);
      if(get(this->Vertex_Sign,
             target(next(next(g, _pMesh), _pMesh), _pMesh)) == NOSIGN)
        put(this->Vertex_Sign,
            target(next(next(g, _pMesh), _pMesh), _pMesh),
            MINUS);
      if(get(this->Vertex_Sign,
             target(next(next(next(g, _pMesh), _pMesh), _pMesh), _pMesh)) ==
         NOSIGN)
        put(this->Vertex_Sign,
            target(next(next(next(g, _pMesh), _pMesh), _pMesh), _pMesh),
            PLUS);
      if(get(this->Vertex_Sign,
             target(next(next(next(next(g, _pMesh), _pMesh), _pMesh), _pMesh),
                    _pMesh)) == NOSIGN)
        put(this->Vertex_Sign,
            target(next(next(next(next(g, _pMesh), _pMesh), _pMesh), _pMesh),
                   _pMesh),
            MINUS);
    }
  }

  else if((type == 14) || (type == 15))
  {
    g = opposite(next(g, _pMesh), _pMesh);
    if(get(this->Facet_Flag, face(g, _pMesh)) != FREE)
      check = true;
    if(get(this->Facet_Flag, face(opposite(next(g, _pMesh), _pMesh), _pMesh)) !=
       FREE)
      check = true;
    if(get(this->Facet_Flag, face(opposite(prev(g, _pMesh), _pMesh), _pMesh)) !=
       FREE)
      check = true;

    if(check == false)
    {
      g = opposite(next(h, _pMesh), _pMesh);
      g = CGAL::Euler::join_face(g, _pMesh);
      g = CGAL::Euler::join_face(g, _pMesh);
      g = CGAL::Euler::join_face(g, _pMesh);

      g = h;
      if(get(this->Vertex_Sign, target(next(g, _pMesh), _pMesh)) == NOSIGN)
        put(this->Vertex_Sign, target(next(g, _pMesh), _pMesh), MINUS);
      if(get(this->Vertex_Sign,
             target(next(next(g, _pMesh), _pMesh), _pMesh)) == NOSIGN)
        put(this->Vertex_Sign,
            target(next(next(g, _pMesh), _pMesh), _pMesh),
            PLUS);
      if(get(this->Vertex_Sign,
             target(next(next(next(g, _pMesh), _pMesh), _pMesh), _pMesh)) ==
         NOSIGN)
        put(this->Vertex_Sign,
            target(next(next(next(g, _pMesh), _pMesh), _pMesh), _pMesh),
            MINUS);
      if(get(this->Vertex_Sign,
             target(next(next(next(next(g, _pMesh), _pMesh), _pMesh), _pMesh),
                    _pMesh)) == NOSIGN)
        put(this->Vertex_Sign,
            target(next(next(next(next(g, _pMesh), _pMesh), _pMesh), _pMesh),
                   _pMesh),
            PLUS);
    }
  }

  return check;
}


template< typename HalfedgeGraph, typename PointMap, typename VertexColorMap >
void
Compression_Valence_Component< HalfedgeGraph, PointMap, VertexColorMap >::
    Get_Coefficient_Up_Quantization(const int &Correct_symbol, int coeff[3])
{
  if(Correct_symbol == 0)
  {
    coeff[0] = -1;
    coeff[1] = -1;
    coeff[2] = -1;
  }
  else if(Correct_symbol == 1)
  {
    coeff[0] = 1;
    coeff[1] = -1;
    coeff[2] = -1;
  }
  else if(Correct_symbol == 2)
  {
    coeff[0] = -1;
    coeff[1] = 1;
    coeff[2] = -1;
  }
  else if(Correct_symbol == 3)
  {
    coeff[0] = 1;
    coeff[1] = 1;
    coeff[2] = -1;
  }
  else if(Correct_symbol == 4)
  {
    coeff[0] = -1;
    coeff[1] = -1;
    coeff[2] = 1;
  }
  else if(Correct_symbol == 5)
  {
    coeff[0] = 1;
    coeff[1] = -1;
    coeff[2] = 1;
  }
  else if(Correct_symbol == 6)
  {
    coeff[0] = -1;
    coeff[1] = 1;
    coeff[2] = 1;
  }
  else if(Correct_symbol == 7)
  {
    coeff[0] = 1;
    coeff[1] = 1;
    coeff[2] = 1;
  }
}


template< typename HalfedgeGraph, typename PointMap, typename VertexColorMap >
void
Compression_Valence_Component< HalfedgeGraph, PointMap, VertexColorMap >::
    init_vertex_attributes(const vertex_descriptor &v)
{
  // ELO-note:
  // some vertex attribute must be initialized,
  // see
  // Mepp1.github/src/components/Compression/Compression_Valence/src/Compression_Valence_Items.h:87
  //   Seed_Edge = -1;
  //   Region_Number = -1;
  //   Removal_Order = -1;
  //   Valid_Vertex = true;
  //   JCW_Move_Error[0] = 0;
  //   JCW_Move_Error[1] = 0;
  //   JCW_Move_Error[2] = 0;
  // only the 3 first are used for now
  put(this->vertex_Seed_Edge, v, -1);
  put(this->vertex_Region_Number, v, -1);
  put(this->vertex_Removal_Order, v, -1);
}


template< typename HalfedgeGraph, typename PointMap, typename VertexColorMap >
void
Compression_Valence_Component< HalfedgeGraph, PointMap, VertexColorMap >::
    truncate_colors(const HalfedgeGraph &_pMesh, VertexColorMap *_v_cm)
{
  typedef typename boost::property_traits< VertexColorMap >::value_type Color;

  // loop over vertices
  vertex_iterator vi_begin = vertices(_pMesh).first;
  vertex_iterator vi_end = vertices(_pMesh).second;
  for(vertex_iterator vi = vi_begin; vi != vi_end; ++vi)
  {
    Color c = get(*_v_cm, *vi);

    if( c[0] <= 1 && c[1] <= 1 && c[2] <= 1 )
      continue; // no need to truncate

    // truncate to 1
    double c0 = std::min(1.0/*double*/, static_cast<double>(c[0]));
    double c1 = std::min(1.0/*double*/, static_cast<double>(c[1]));
    double c2 = std::min(1.0/*double*/, static_cast<double>(c[2]));

    put(*_v_cm, *vi, Color(c0, c1, c2));
    //DBG put(*_v_cm, *vi, Color(0, 0, 0));
   }
}


template< typename HalfedgeGraph, typename PointMap, typename VertexColorMap >
void
Compression_Valence_Component< HalfedgeGraph, PointMap, VertexColorMap >::
    DBG_print_mesh_geometry(const HalfedgeGraph &_pMesh,
                            const PointMap *_pm,
                            const std::string &header)
{
  std::cout << header << "  mesh-geometry:" << std::endl;

  size_t count = 0;
  auto vertex_iterator_pair = vertices(_pMesh);
  vertex_iterator pVertex = vertex_iterator_pair.first;
  vertex_iterator pVertex_end = vertex_iterator_pair.second;
  for(; pVertex != pVertex_end; ++pVertex)
  {
    count++;
    auto p3d = get(*_pm, *pVertex);
    std::cout << "V#" << count << " point: " << p3d[0] << " " << p3d[1] << " "
              << p3d[2] << std::endl;
  }
}

template< typename HalfedgeGraph, typename PointMap, typename VertexColorMap >
void
Compression_Valence_Component< HalfedgeGraph, PointMap, VertexColorMap >::
    DBG_print_mesh_vertexcolor(const HalfedgeGraph &_pMesh,
                               const VertexColorMap *_v_cm,
                               const std::string &header)
{
  typedef typename boost::property_traits< VertexColorMap >::value_type Color;

  std::cout << header << "  mesh-vertexcolor:" << std::endl;

  if(!_v_cm)
  {
    std::cout << header << "THE MESH HAS NO VERTEX COLOR MAP!" << std::endl;
    return;
  }

  size_t count = 0;
  auto vertex_iterator_pair = vertices(_pMesh);
  vertex_iterator pVertex = vertex_iterator_pair.first;
  vertex_iterator pVertex_end = vertex_iterator_pair.second;
  for(; pVertex != pVertex_end; ++pVertex)
  {
    count++;
    Color c = get(*_v_cm, *pVertex);
    std::cout << "V#" << count << " color: " << c[0] << " " << c[1] << " "
              << c[2] << std::endl;
  }
}


#ifdef _MSC_VER
#pragma warning(pop)
#endif
