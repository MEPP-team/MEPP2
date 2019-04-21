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

#include "arithmetic_codec.hpp"

#include <boost/graph/graph_traits.hpp>
#include "FEVV/Wrappings/Geometry_traits.h"
#include "FEVV/Wrappings/properties.h"

#include <queue>
#include <list>

// struct of integer coordinates

////////////////////////////////////////////////////////////////////////////////////////////////////
/// \struct	Point_Int
///
/// \brief	Point.
///
////////////////////////////////////////////////////////////////////////////////////////////////////

struct Point_Int
{

  int x; ///< The x coordinate
  int y; ///< The y coordinate
  int z; ///< The z coordinate

  ////////////////////////////////////////////////////////////////////////////////////////////////////
  ///
  /// \brief	Default constructor.
  ///
  ////////////////////////////////////////////////////////////////////////////////////////////////////

  Point_Int() { x = y = z = 0; }

  ////////////////////////////////////////////////////////////////////////////////////////////////////
  ///
  /// \brief	Addition operator.
  ///
  /// \param	Pt	The point.
  ///
  /// \return	The result of the operation.
  ////////////////////////////////////////////////////////////////////////////////////////////////////

  const Point_Int operator+(const Point_Int &Pt) const
  {
    Point_Int Res;
    Res.x = x + Pt.x;
    Res.y = y + Pt.y;
    Res.z = z + Pt.z;

    return Res;
  }

  ////////////////////////////////////////////////////////////////////////////////////////////////////
  ///
  /// \brief	Negation operator.
  ///
  /// \param	Pt	The point.
  ///
  /// \return	The result of the operation.
  ////////////////////////////////////////////////////////////////////////////////////////////////////

  const Point_Int operator-(const Point_Int &Pt) const
  {
    Point_Int Res;
    Res.x = x - Pt.x;
    Res.y = y - Pt.y;
    Res.z = z - Pt.z;

    return Res;
  }

  ////////////////////////////////////////////////////////////////////////////////////////////////////
  ///
  /// \brief	Comparison operatior ==.
  ///
  /// \param	Pt	The point.
  ///
  /// \return	true if the parameters are considered equivalent.
  ////////////////////////////////////////////////////////////////////////////////////////////////////

  bool operator==(const Point_Int &Pt) const
  {
    return (x == Pt.x && y == Pt.y && z == Pt.z);
  }

  ////////////////////////////////////////////////////////////////////////////////////////////////////
  ///
  /// \brief	Comparison operatior !=.
  ///
  ///
  /// \param	Pt	The point.
  /// \return true if the parameters are not identical.
  ////////////////////////////////////////////////////////////////////////////////////////////////////

  bool operator!=(const Point_Int &Pt) const { return !(*this == Pt); }
};


////////////////////////////////////////////////////////////////////////////////////////////////////
/// \struct	Color_Unit
///
/// \brief	Struct of integer color components.
///
////////////////////////////////////////////////////////////////////////////////////////////////////

struct Color_Unit
{
  int c0;
  int c1;
  int c2;

  Color_Unit() { c0 = c1 = c2 = 0; }

  // ELO+ begin
  // constructor to easily store data in property maps
  Color_Unit(int _c0, int _c1, int _c2)
  {
    c0 = _c0;
    c1 = _c1;
    c2 = _c2;
  }
  // ELO+ end


  ////////////////////////////////////////////////////////////////////////////////////////////////////
  ///
  /// \brief	Addition operator.
  ///
  /// \param	m_color	The color.
  ///
  /// \return	The result of the operation.
  ////////////////////////////////////////////////////////////////////////////////////////////////////

  const Color_Unit operator+(const Color_Unit &m_color) const
  {
    Color_Unit Res;
    Res.c0 = c0 + m_color.c0;
    Res.c1 = c1 + m_color.c1;
    Res.c2 = c2 + m_color.c2;

    return Res;
  }

  ////////////////////////////////////////////////////////////////////////////////////////////////////
  ///
  /// \brief	Negation operator.
  ///
  /// \param	m_color	The color.
  ///
  /// \return	The result of the operation.
  ////////////////////////////////////////////////////////////////////////////////////////////////////

  const Color_Unit operator-(const Color_Unit &m_color) const
  {
    Color_Unit Res;
    Res.c0 = c0 - m_color.c0;
    Res.c1 = c1 - m_color.c1;
    Res.c2 = c2 - m_color.c2;

    return Res;
  }

  ////////////////////////////////////////////////////////////////////////////////////////////////////
  ///
  /// \brief	Comparison operatior ==.
  ///
  /// \param	m_color	The col.
  ///
  /// \return	true if the parameters are considered equivalent.
  ////////////////////////////////////////////////////////////////////////////////////////////////////

  bool operator==(const Color_Unit &m_color) const
  {
    return (c0 == m_color.c0 && c1 == m_color.c1 && c2 == m_color.c2);
  }

  ////////////////////////////////////////////////////////////////////////////////////////////////////
  ///
  /// \brief	Comparison operatior !=.
  ///
  /// \param	m_color	The col.
  ////////////////////////////////////////////////////////////////////////////////////////////////////

  bool operator!=(const Color_Unit &m_color) const
  {
    return !(*this == m_color);
  }
};

////////////////////////////////////////////////////////////////////////////////////////////////////
/// \class	Compression_Valence_Component
///
/// \brief	Compression valence component.
///
////////////////////////////////////////////////////////////////////////////////////////////////////

template< typename HalfedgeGraph, typename PointMap, typename VertexColorMap >
class Compression_Valence_Component
{
public:
  typedef boost::graph_traits< HalfedgeGraph > GraphTraits;
  typedef typename GraphTraits::vertex_descriptor vertex_descriptor;
  typedef typename GraphTraits::vertex_iterator vertex_iterator;
  typedef typename GraphTraits::face_descriptor face_descriptor;
  typedef typename GraphTraits::face_iterator face_iterator;
  typedef typename GraphTraits::halfedge_descriptor halfedge_descriptor;
  typedef typename GraphTraits::halfedge_iterator halfedge_iterator;
  typedef typename FEVV::Geometry_traits< HalfedgeGraph >::Vector Vector;
  typedef typename FEVV::Geometry_traits< HalfedgeGraph >::Point Point3d;


  ////////////////////////////////////////////////////////////////////////////////////////////////////
  ///
  /// \brief	Default Constructor.
  ///
  ////////////////////////////////////////////////////////////////////////////////////////////////////

  Compression_Valence_Component(void)
  {
    IsOneColor = false;
    Decompress_count = 0;

    Smallest_Alpha = 5000;
    Smallest_Gamma = 5000;

    Smallest_C0 = 5000;
    Smallest_C1 = 5000;
    Smallest_C2 = 5000;

    /*ColorDiffMinC0 = 5000;
    ColorDiffMinC1 = 5000;
    ColorDiffMinC2 = 5000;*/

    IsColored = false;
    GlobalCountOperation = 0;

    // from IHM
    IsCompressed = false;
    IsDecompress = false;
    // Afficher = false;

    Possible_change_sequence = true;
    // ELO  sequence mode is enabled by default on Mepp1
    // ELO  but disabled by default on Mepp2 because mesh
    // ELO  sequence is not yet supported
    Sequence = false; // ELO- Sequence = true;
    Visu_level = 0;
    Process_level = 0;
    // Writing_level = -1;

    N_Inserted_Watermarks = 0;

    m_VC[0] = 0.;
    m_VC[1] = 0.;
    m_VC[2] = 0.;
    m_Dist = 0.0005;

    Number_non_reversible_vertices = 0;
  }

  /**
   \brief	Destructor.
   */

  ~Compression_Valence_Component() {}


  // Main Function

  /**
   \brief  Main function of compression.

   \param [in,out] pMesh            The mesh.
   \param          Input_File_Name  Filename of the input file.
   \param          File_Name        Filename of the compressed(output) file.
   \param          Qbit             The qbit.
   \param          NVertices        Maximum number of vertices after
                                    simplification/compression.
   \param          Is_normal_flipping_selected
                                    The normal flipping.
   \param          Is_use_metric_selected
                                    The use metric.
   \param          Metric_thread    The metric thread.
   \param	         Is_use_forget_metric_selected
                                    The use forget metric.
   \param          Forget_value     The forget value.
   \param	         Is_compression_selected
                                    The compression selected.
   \param          Is_adaptive_quantization_selected
                                    The adaptive quantization.
   \param          Is_bijection_selected
                                    The bijection of the geometry coding.
   \return         Information of compression results.
   */

  std::string Main_Function(HalfedgeGraph &pMesh,
                            PointMap *_pm,
                            VertexColorMap *_v_cm,
                            const std::string &Input_File_Name,
                            const std::string &File_Name,
                            const int &Qbit,
                            const int &NVertices,
                            const bool Is_normal_flipping_selected,
                            const bool Is_use_metric_selected,
                            const float &Metric_thread,
                            const bool Is_use_forget_metric_selected,
                            const int &Forget_value,
                            const bool Is_compression_selected,
                            const bool Is_adaptive_quantization_selected,
                            const bool Is_bijection_selected);

  // Initialization

  /**
   \brief	Global initialization to select the input gate.

   \param [in,out] pMesh             The mesh.
   \param          Quantization_bit  Number of bits used for geometry
                                     quantization.
   \param          File_name         Filename of the file.
   */

  void Global_Initialization(HalfedgeGraph &pMesh,
                             const int &Quantization_bit,
                             const char *File_name,
                             const PointMap *pm,
                             VertexColorMap *v_cm);

  /**
          \brief	Quantize all vertices so that the new positions are
     reguliraly spaced in the 3D space.

          \param [in,out]	pMesh	The mesh.
          */

  void Quantization(HalfedgeGraph &pMesh, const PointMap *pm);

  /**
          \brief	Color initialization.

          \param [in,out]	pMesh	The mesh.
          */

  void Color_Initialization(HalfedgeGraph &pMesh, VertexColorMap *v_cm);

  /**
          \brief	Initialization to deal with a mesh composed of multiple
     separated components.

          \param [in,out]	pMesh	The mesh.
          \param	Quantization_bit		The number of quantization for
     geometry.
          */

  void Multiple_Components_Initialization(HalfedgeGraph &pMesh,
                                          const PointMap *pm,
                                          const int &Quantization_bit);

  /**
   \brief	Initialize all flags of verticeces and facets to FREE and give
  order to vertices (manifold property check). This function is called every
  conquest.

   \param [in,out]	pMesh	The mesh.
   */
  void Init(HalfedgeGraph &mesh);

#if 0 // TODO-elo-restore-if-needed
		/**
			\brief	Color quantization.

			\param [in,out]	pMesh	The mesh.
			*/

		void Color_Quantization(HalfedgeGraph &pMesh);
#endif

  /**
   \brief  Mesh Simplification which applies iteratively decimation and
           regulation in pair.

   \param [in,out] pMesh      The mesh.
   \param          NVertices  The vertices.
   \param          Is_normal_flipping_selected
                              The normal flipping.
   \param          Is_use_metric_selected
                              The use metric.
   \param          Metric_thread
                              The metric thread.
   \param          Is_use_forget_metric_selected
                              The use forget metric.
   \param          Forget_value
                              The forget value.
   */

  void Simplification(HalfedgeGraph &pMesh,
                      const PointMap *_pm,
                      const int &NVertices,
                      const bool Is_normal_flipping_selected,
                      const bool Is_use_metric_selected,
                      const float &Metric_thread,
                      const bool Is_use_forget_metric_selected,
                      const int &Forget_value);

  /**
          \brief	Removal of a set of independent vertices.

          \param [in,out]	pMesh	The mesh.
          \param	Is_normal_flipping_selected  	The normal flipping.
          \param	Is_use_metric_selected		 	The use metric.
          \param	Metric_thread	 	The metric thread.
          \param	Is_use_forget_metric_selected	The use forget metric.
          \param	Forget_value	 	The forget value.
          \param	Component_ID	 	Identifier for the component.

          \return	Number of decimated vertices.
          */

  int Decimation_Conquest(HalfedgeGraph &pMesh,
                          const PointMap *_pm,
                          const bool Is_normal_flipping_selected,
                          const bool Is_use_metric_selected,
                          const float &Metric_thread,
                          const bool Is_use_forget_metric_selected,
                          const int &Forget_value,
                          const int &Component_ID);

  /**
          \brief	Removal of a set of independent vertices.

          \param [in,out]	pMesh	The mesh.
          \param	Is_normal_flipping_selected  	The normal flipping.
          \param	Is_use_metric_selected		 	The use metric.
          \param	Metric_thread	 	The metric thread.
          \param	Is_use_forget_metric_selected	The use forget metric.
          \param	Forget_value	 	The forget value.
          \param	Component_ID	 	Identifier for the component.

          \return	Number of decimated vertices.
          */
  int Regulation(HalfedgeGraph &_pMesh,
                 const bool Normal_flipping,
                 const bool Use_metric,
                 const float &Metric_thread,
                 const bool Use_forget_metric,
                 const int &Forget_value,
                 const int &Component_ID,
                 const PointMap *_pm);


  // Adaptive Quantization

  /**
   \brief	Adaptive quantization which not only decimates a input mesh but
  also adapts quantization precision for all intermediate meshes.

   \param [in,out] pMesh      The mesh.
   \param          NVertices  The vertices.
   \param          Is_normal_flipping_selected
                              The normal flipping.
   \param          Is_use_metric_selected
                              The use metric.
   \param          Metric_thread
                              The metric thread.
   \param          Is_use_forget_metric_selected
                              The use forget metric.
   \param          Forget_value
                              The forget value.
   \param          Qbit       The qbit.
   */

  void Adaptive_Quantization(HalfedgeGraph &pMesh,
                             PointMap *_pm,
                             VertexColorMap *_v_cm,
                             const int &NVertices,
                             const bool Is_normal_flipping_selected,
                             const bool Is_use_metric_selected,
                             const float &Metric_thread,
                             const bool Is_use_forget_metric_selected,
                             const int &Forget_value,
                             const int &Qbit);


  /**
   \brief	Refines quantization precision of mesh geometry.

   \param [in,out] pMesh         The mesh.
   \param [in,out] Decoder       The decoder.
   \param          Component_ID  Component ID.
   */

  void Augment_Geometry_Quantization_Precision(HalfedgeGraph &pMesh,
                                               Arithmetic_Codec &Decoder,
                                               const int &Component_ID,
                                               PointMap *_pm);

  /**
   \brief  Decreasing of quantization resolution based on the prediction of
           PENG. Opposite function is Augment_Geometry_Quantization_Precision.

   \param [in,out] pMesh         The mesh.
   \param          Component_ID  Component ID.
   */


  void Diminush_Geometry_Quantization_Precision(HalfedgeGraph &pMesh,
                                                const int &Component_ID,
                                                PointMap *_pm);

  /**
   \brief  Calculates the edge color difference.

   \param [in,out] pMesh               The mesh.
   \param          Component_ID        Identifier for the component.
   \param [in,out] Max_color           The maximum of color difference.
   \param [in,out] Mean_color          The mean of color difference.
   \param [in,out] Number_of_vertices  Number of vertices.
   */

  void Calculate_Edge_Color_Difference(HalfedgeGraph &pMesh,
                                       const int &Component_ID,
                                       double &Max_color,
                                       double &Mean_color,
                                       int &Number_of_vertices);

  /**
   \brief Reduces quantization precision of color coordinates.

   \param [in,out]	pMesh	The mesh.
   \param	Component_ID	 	Identifier for the component.
   */

  void Diminush_Color_Quantization_Precision(HalfedgeGraph &pMesh,
                                             const int Component_ID,
                                             VertexColorMap *_v_cm);

  /**
   \brief	Refine quantization precision of mesh color coordinates.

   \param [in,out]	pMesh  	The mesh.
   \param [in,out]	Decoder	The decoder.
   \param	Component_ID	   	Identifier for the component.
   */

  void Augment_Color_Quantization_Precision(HalfedgeGraph &_pMesh,
                                            Arithmetic_Codec &Decoder,
                                            const int &Component_ID,
                                            VertexColorMap *_v_cm);

  // Compression

  /**
   \brief  Compressions.

   \param [in,out] pMesh              The mesh.
   \param          File_Name          Filename of the output file.
   \param          Qbit               The quantization of geometry.
   \param [in,out] Connectivity_size  Size of the connectivity.
   \param [in,out] Color_size         Size of the color.
   \param [in,out] Total_size         Size of the total.
   */

  void Compression(HalfedgeGraph &pMesh,
                   const char *File_Name,
                   const int &Qbit,
                   unsigned &Connectivity_size,
                   unsigned &Color_size,
                   unsigned &Total_size,
                   const PointMap *_pm);

#if 0 // TODO-elo-restore-if-needed
		/**
		 \brief	Calculates the connectivity rate.


		 \return	The calculated connectivity rate.
		 */

		int  Calculate_Connectivity_Rate(void);
#endif

  /**
   \brief	Calculates the geometry color offset range.
          This function is needed since the coder/decoder deals only with the
  positive symbol numbers.

   */
  void Calculate_Geometry_Color_Offset_Range(void);

  /**
   \brief	When the last simplification operation has to be cancelled,
          we need to remove its elements in order not to compress these
   elements.

   \param	Component_ID	Identifier for the component.
   */

  void Remove_Last_Phase_Elements(const int &Component_ID);

  /**
   \brief	Writes a base mesh.

   \param [in,out]	pMesh			 	The mesh.
   \param [in,out]	Enc				 	The encoder.
   \param [in,out]	Connectivity_size	Size of the connectivity.
   \param [in,out]	Color_size		 	Size of the color.
   \param	Num_color_base_mesh			 	Number of color in the
   base mesh.
   */

  void Write_Base_Mesh(const HalfedgeGraph &pMesh,
                       Arithmetic_Codec &Enc,
                       unsigned &Connectivity_size,
                       unsigned &Color_size,
                       const int &Num_color_base_mesh,
                       const PointMap *_pm);


  // Decompression

  /**
   \brief	Initialize the Decompression by loading the base mesh into
  pMesh. aka decompress the first (simplest) level of the mesh.

   \param [in,out]	pMesh			 	The mesh.

   \return	Information related to the decompression (number of LoDs, size of
  the input file, etc).
   */

  void Decompress_Init(HalfedgeGraph &pMesh,
                       PointMap *_pm,
                       VertexColorMap *_v_cm,
                       const std::string &Input_File_Name);

#if 0 // TODO-elo-restore-if-needed
		/**
		 \brief	Decompression from file (No creation of mesh sequence).

		 \param [in,out]	pMesh	The mesh.
		 */

		void Decompression_From_File(HalfedgeGraph &pMesh);
#endif

  /**
   \brief	Write information of each LoD at decompression in a file.
   \param [in,out]	pMesh	The mesh.

   */
  void Write_Info(HalfedgeGraph &pMesh);

#if 0 // TODO-elo-restore-if-needed
		/**
		 \brief	Decompression from sequence (Creation of mesh sequence).

		 \param [in,out]	pMesh	The mesh.
		 \param [in,out]	New_mesh	The new copied mesh.
		 */
		void Decompression_From_Sequence(HalfedgeGraph &pMesh, HalfedgeGraph &New_mesh);


		/**
		 \brief	Obtain the previous LoD in the "from file" mode (no mesh sequence).

		 \param [in,out]	pMesh	The mesh.
		 */
		void Decompression_Coarser_From_File(HalfedgeGraph &pMesh);
#endif

  /**
   \brief  Write intermediate mesh to file.
  */
  void write_intermediate_mesh(/*const*/ HalfedgeGraph &_pMesh,
                               const VertexColorMap *_v_cm);

  /**
   \brief  Copy the current mesh and vertex-color map to
           a new mesh and vertex-color map.
           If v_cm_copy == nullptr, the vertex-color map
           is NOT copied.
  */
  static void copy_mesh(const HalfedgeGraph &_pMesh,
                        const VertexColorMap *_v_cm,
                        HalfedgeGraph &mesh_copy,
                        VertexColorMap *v_cm_copy /* nullptr allowed */);

  /**
   \brief  Store a copy of the current mesh and vertex color map to be
           able to display all levels of decompression.
  */
  void keep_intermediate_mesh(
      const HalfedgeGraph &_pMesh,
      const VertexColorMap *_v_cm,
      std::vector< HalfedgeGraph * > *intermediate_meshes,
      std::vector< VertexColorMap * >
          *intermediate_vertexColorMaps /* nullptr allowed */);


  /**
   \brief	Decompression of all LoDs from file, or until a specified level.
          The finest LoD is visualized without creating mesh sequence.

   \param [in,out]	pMesh	The mesh.
   */
  std::string Decompression_All_From_File(
      HalfedgeGraph &_pMesh,
      PointMap *_pm,
      VertexColorMap *_v_cm,
      bool do_write_info,
      std::vector< HalfedgeGraph * > *intermediate_meshes /* nullptr allowed */,
      std::vector< VertexColorMap * >
          *intermediate_vertexColorMaps /* nullptr allowed */,
      int stop_level = -1 /* decompression level to stop at */,
      bool do_write_intermediate_meshes = false);


#if 0 // TODO-elo-restore-if-needed
		/**
		 \brief	To obtain a user-desired LoD from file (No creation of mesh sequence).

		 \param [in,out]	pMesh	The mesh.
		 \param             WL      The wanted level.
		 */
		void Decompression_Specific_Level_From_File(HalfedgeGraph &pMesh, const int & WL);

		/**
		 \brief	Decompression from file for JCW (No creation of mesh sequence).

		 \param [in,out]	pMesh	The mesh.
		 */
		void JCW_Decompression_From_File(HalfedgeGraph &pMesh);


		/**
		 \brief	Decompression from file without watermark extraction and geometry correction (No creation of mesh sequence).

		 \param [in,out]	pMesh	The mesh.
		 */
		void JCW_Decompression_Without_Extraction_From_File(HalfedgeGraph &pMesh);


		/**
		 \brief	Decompression from sequence for JCW (Creation of mesh sequence).

		 \param [in,out]	pMesh	The mesh.
		 \param [in,out]	New_mesh	The new copied mesh.
		 */
		void JCW_Decompression_From_Sequence(HalfedgeGraph &pMesh, HalfedgeGraph &New_mesh);


		/**
		 \brief	Decompression from sequence without watermark extraction (Creation of mesh sequence).

		 \param [in,out]	pMesh	The mesh.
		 \param [in,out]	New_mesh	The new copied mesh.
		 */

		void JCW_Decompression_Without_Extraction_From_Sequence(HalfedgeGraph &pMesh, HalfedgeGraph &New_mesh);

		/**
		 \brief	Obtain information to the main window.

		 \return Information of each LoD at decompression.

		 */
		QString Show_Text(void);
#endif


  /**
   \brief	Decompress the each step to visualize intermediate meshes

   \param [in,out]	pMesh	The mesh.
   \param	File_Name	Filename of the file.

   \return Index of the current LoD.
   */

  int Decompress_Each_Step(HalfedgeGraph &_pMesh,
                           PointMap *_pm,
                           VertexColorMap *_v_cm);

  /**
   \brief	Decoding of the decimation conquest.

   \param [in,out]	pMesh	The mesh.
   \param [in,out]	Decoder  	The decoder.
   \param	Component_ID			Component ID.
   */

  void Un_Decimation_Conquest(HalfedgeGraph &pMesh,
                              Arithmetic_Codec &Decoder,
                              const int &Component_ID,
                              PointMap *_pm,
                              VertexColorMap *_v_cm);

  /**
   \brief	Decoding of the regulation conquest.

   \param [in,out]	pMesh	The mesh.
   \param [in,out]	Decoder  	The decoder.
   \param	Component_ID				 	Component ID.
   */

  void Un_Regulation(HalfedgeGraph &_pMesh,
                     Arithmetic_Codec &Decoder,
                     const int &Component_ID,
                     PointMap *_pm,
                     VertexColorMap *_v_cm);

#if 0 // TODO-elo-restore-if-needed
      // Other functions
      // void   Separate_Components(HalfedgeGraph &pMesh);

		/**
		 \brief	Calculates the area of pMesh.

		 \param [in,out]	pMesh	The mesh.

		 \return	The calculated area.
		 */

		double Calculate_Area(HalfedgeGraph & pMesh);
#endif

  /**
   \brief	Calculate an error cause by a vertex removal and decide in order
   not to remove an important vertex.

   \param [in,out]	pMesh	The mesh.
   \param h				Input gate.
   \param Component_ID    Component ID.
   \param Mean_color      Mean color value.
   \param Mean_area       Average of facets areas.
   \return	Decision of vertex removal.
   */
  bool Error_Projected_Surface(const HalfedgeGraph &_pMesh,
                               const PointMap *_pm,
                               const halfedge_descriptor &_h,
                               const int &_Component_ID,
                               const double &Mean_color,
                               const double &Mean_area);

  /**
   \brief	Update the area of each component.

   \param [in,out]	pMesh	The mesh.
   \param Component_ID  Component ID.
   \param Number_facets Number of facets

   \return	Decision of vertex removal.
   */
  void Recalculate_Component_Area(HalfedgeGraph &pMesh,
                                  const PointMap *_pm,
                                  const int &Component_ID,
                                  int &Number_facets);


  /**
   \brief	Change from real to int point(related to quantization bit)

   \param	pt				The point.
   \param	Component_ID	Component ID.

   \return Integer coordinates.
   */
  inline Point_Int Change_Real_Int(const Point3d &pt, const int &Component_ID);

  /**
   \brief	Change from int to real point(related to quantization bit)

   \param	pt				The point(integers).
   \param	Component_ID	Component ID.

   \return	.
   */

  inline Point3d Change_Int_Real(const Point_Int &pt, const int &Component_ID);


#if 0 // TODO-elo-restore-if-needed
		/**
		 \brief	To add flags to copied mesh


		 \param [in,out]	Original_mesh	The original mesh.
		 \param [in,out]	New_mesh	 	The new mesh.
		 */

		void Attibute_Seed_Gate_Flag(HalfedgeGraph &Original_mesh, HalfedgeGraph &New_mesh);
#endif

  /**
   \brief	Calculates the current file size. Measure bits used for
  decompression.

   \return	The calculated current file size.
   */

  unsigned Calculate_Current_File_Size(void)
  {
    // To measure exact quantity of bits used for decompression.
    unsigned Adjust_value;
    if(this->IsColored)
      Adjust_value = 25 + 9 * 4;
    else
      Adjust_value = 25;

    return this->Decoder.calculate_current_decoded_size() + Adjust_value;
  }

#if 0 // TODO-elo-restore-if-needed
		/**
		 \brief	Gets a resolution change.

		 \param [in,out]	pMesh	If non-null, the mesh.
		 \param	Prec			 	The prec.

		 \return	The resolution change.
		 */

		int GetResolutionChange(HalfedgeGraph *pMesh, float Prec);

		/**
		 \brief	To stop the decoder (Used "Previous level" while decoding).

		 */

		void Stop_Decoder(void) {this->Decoder.stop_decoder(); }

		//////////////////////////////////////////
		// Joint Compression Watermarking (JCW) //
		//////////////////////////////////////////

		/**
		 \brief	Joint Compression Watermarking (JCW)

		 \param [in,out]	pMesh			 		The mesh.
		 \param Input_File_Name						The name of the input file.
		 \param Output_File_Name					The name of the output file.
		 \param Number_bins							The number of bins.
		 \param Number_regions						The number of regions.
		 \param Embedding_strength				    The strengh of embedding (number of shifted bins).
		 \param Embedding_message                   The message of watermarking
		 \param Is_complete_reversibility_selected  The selection of complete reversibility.
		 \param Is_divide_regions_selected          The selection of division of regions.
		 \param Thres_divide_regions                The threshold to divide big regions.
		 \param Qbit                                The geometry quantization.
		 \param NVertices                           The wanted number of base mesh.
		 \param Normal_flipping                     The selection of normal flipping.
		 \param Use_metric                          The selection of geometric metric use
		 \param Metric_thread                       The threshold of metric.
		 \param Use_forget_metric                   The selection of use of "forget metric"
		 \param Forget_value                        The threshold of "Use_forget_metric"

		 \return	Information of JCW.
		 */

		QString Joint_Compression_Watermarking(HalfedgeGraph &pMesh,
											  const char * Input_File_Name,
											  const char * Output_File_Name,
											  const int & Number_bins,
											  const int & Number_regions,
											  const int & Embedding_strength,
											  const char * Embedding_message,
											  const bool Is_complete_reversibility_selected,
											  const bool Is_divide_regions_selected,
											  const int & Thres_divide_regions,
											  const int &Qbit,
											  const int  & NVertices,
											  const bool Normal_flipping,
											  const bool Use_metric,
											  const float & Metric_thread,
											  const bool Use_forget_metric,
											  const int &Forget_value);


		/**
		 \brief	JCW decimation for segmentation.

		 \param [in,out]	pMesh	The mesh.
		 \param	Normal_flipping  	The normal flipping.
		 \param	Use_metric		 	The use metric.
		 \param	Metric_thread	 	The metric thread.
		 \param	Use_forget_metric	The use forget metric.
		 \param	Forget_value	 	The forget value.
		 \param	Component_ID	 	Identifier for the component.

		 \return	.
		 */

		int JCW_Decimation_For_Segmentation(HalfedgeGraph &pMesh,const bool Normal_flipping,const bool Use_metric,const float &Metric_thread, const bool Use_forget_metric,const int &Forget_value, const int & Component_ID);

		/**
		 \brief	JCW regulation for segmentation.


		 \param [in,out]	pMesh	The mesh.
		 \param	Normal_flipping  	The normal flipping.
		 \param	Use_metric		 	The use metric.
		 \param	Metric_thread	 	The metric thread.
		 \param	Use_forget_metric	The use forget metric.
		 \param	Forget_value	 	The forget value.
		 \param	Component_ID	 	Identifier for the component.

		 \return	.
		 */

		int JCW_Regulation_For_Segmentation(HalfedgeGraph &pMesh,const bool Normal_flipping,const bool Use_metric,const float &Metric_thread, const bool Use_forget_metric,const int &Forget_value, const int & Component_ID);

        /**
         \brief	JCW un regulation for insertion.

         \param [in,out]	pMesh				   	The mesh.
         \param	Component_ID					   	Component ID.
         \param [in,out]	FP_Connectivity		   	The fp connectivity.
         \param [in,out]	SP_Moved_Position	   	The sp moved position.
         \param [in,out]	SP_Original_Position   	The sp original position.
         \param [in,out]	SP_Watermarked_Position	The sp watermarked position.
         \param [in,out]	JCW_ERROR			   	The jcw error.
         */

        void JCW_Un_Regulation_For_Insertion(HalfedgeGraph &pMesh, const int & Component_ID, list<int> & FP_Connectivity, list<Point3d> & SP_Moved_Position, list<Point3d> & SP_Original_Position, list<Point_Int> & SP_Watermarked_Position, list<vector<int> > & JCW_ERROR);

		/**
		 \brief	JCW un decimation for insertion.


		 \param [in,out]	pMesh				   	The mesh.
		 \param	Component_ID					   	Component ID
		 \param [in,out]	FP_Connectivity		   	The fp connectivity.
		 \param [in,out]	SP_Moved_Position	   	The sp moved position.
		 \param [in,out]	SP_Original_Position   	The sp original position.
		 \param [in,out]	SP_Watermarked_Position	The sp watermarked position.
		 \param [in,out]	JCW_ERROR			   	The jcw error.
		 */

		void JCW_Un_Decimation_For_Insertion(HalfedgeGraph &pMesh, const int & Component_ID, list<int> & FP_Connectivity, list<Point3d> & SP_Moved_Position, list<Point3d> & SP_Original_Position, list<Point_Int> & SP_Watermarked_Position, list<vector<int> > & JCW_ERROR);

		/**
		 \brief	JCW barycenter patch before removal.

		 \param	h		 	The input gate.
		 \param	Direction	The direction.

		 \return Predicted position.
		 */

		Point3d JCW_Barycenter_Patch_Before_Removal(const Halfedge_handle & h, const int & Direction);

		/**
		 \brief	Caclulates the barycentric position for JCW after removal.


		 \param	h		 	The input gate.
		 \param	valence  	The valence.
		 \param	Direction	The direction.

		 \return Predicted position.
		 */

		Point3d JCW_Barycenter_Patch_After_Removal(const Halfedge_handle & h, const int & valence, const int & Direction);

		/**
		 \brief	JCW barycenter patch before removal.


		 \param	h		 	The input gate.
		 \param	valence  	The valence.
		 \param	Direction	The direction.

		 \return	.
		 */

		Point3d JCW_Barycenter_Patch_Before_Removal(const Halfedge_handle & h, const int & valence, const int & Direction);

		/**
		 \brief	Undecimation conquest for JCW.

		 \param [in,out]	pMesh	The mesh.
		 \param [in,out]	Decoder  	The decrement.
		 \param	Component_ID				 	Component ID.
		 */

		void JCW_Un_Decimation_Conquest(HalfedgeGraph &pMesh,Arithmetic_Codec & Decoder, const int & Component_ID);

		/**
		 \brief	JCW unregulation.


		 \param [in,out]	pMesh	The mesh.
		 \param [in,out]	Decoder  	The decrement.
		 \param	Component_ID				 	Component ID.
		 */

		void JCW_Un_Regulation(HalfedgeGraph &pMesh, Arithmetic_Codec & Decoder, const int & Component_ID);

		/**
		 \brief	JCW unregulation for region detection, in order to obtain region number of each inserted vertices

 		 \param [in,out]	pMesh	  	The mesh.
		 \param	Component_ID		  	Identifier for the component.
		 \param [in,out]	FP_connect	The fp connect.
		 \param [in,out]	FP_Geo	  	The fp geo.
		 \param [in,out]	FP_RN	  	The fp rn.
		 */

		void JCW_Un_Regulation_For_Region_Detection(HalfedgeGraph & pMesh, const int & Component_ID, list<int> & FP_connect, list<Point3d> & FP_Geo, list<int> & FP_RN);

		/**
		 \brief	JCW undecimation for region detection.


		 \param [in,out]	pMesh	  	The mesh.
		 \param	Component_ID		  	Identifier for the component.
		 \param [in,out]	FP_connect	The fp connect.
		 \param [in,out]	FP_Geo	  	The fp geo.
		 \param [in,out]	FP_RN	  	The fp rn.
		 */

		void JCW_Un_Decimation_For_Region_Detection(HalfedgeGraph & pMesh, const int & Component_ID, list<int> & FP_connect, list<Point3d> & FP_Geo, list<int> & FP_RN);

		/**
		 \brief	Evaluates the JCW robustness.


		 \return Results of robustness evaluation.
		 */

		vector<double> JCW_Evaluate_Robustness(void);

		/**
		 \brief	Sets the bin number to NB.

		 \param	NB	The nb.
		 */

		void Set_Number_Bin(const int & NB)
		{
			this->m_NumberBin = NB;
		}

		/**
		 \brief	Gets the bin number.

		 \return	The number bin.
		 */

		int Get_Number_Bin(void)
		{
			return this->m_NumberBin;
		}

		/**
		 \brief	Sets the embedding level.


		 \param	EL	The el.
		 */

		void Set_Embedding_Level(const int &EL)
		{
			this->m_EmbeddingStrength = EL;
		}

		/**
		 \brief	Gets the embedding level.


		 \return	The embedding level.
		 */

		int Get_Embedding_Level(void)
		{
			return this->m_EmbeddingStrength;
		}

		/**
		 \brief	Sets the region number.

			 \param	NR	The nr.
		 */

		void Set_Number_Region(const int &NR)
		{
			this->m_NumberRegion = NR;
		}

		/**
		 \brief	Gets the region number.

		 \return	The number region.
		 */

		int Get_Number_Region(void)
		{
			return this->m_NumberRegion;
		}

		/**
		 \brief	 calculates the mesh center for JCW.

		 \param [in,out]	pMesh	The mesh.
		 */

		void JCW_Calculate_Mesh_Center(HalfedgeGraph &pMesh);

		/**
		 \brief	Calculateq the radius for JCW.

		 \param [in,out]	pMesh	The mesh.
		 */
		void JCW_Calculate_Radius(HalfedgeGraph &pMesh);

		/**
		 \brief	JCW expand mesh.


		 \param [in,out]	pMesh	The mesh.
		 */

		void JCW_Expand_Mesh(HalfedgeGraph &pMesh);

		/**
		 \brief	quantization for JCW.


		 \param [in,out]	pMesh	The mesh.
		 */

		void JCW_Quantization(HalfedgeGraph &pMesh);

		/**
		 \brief	Generate regions on base mesh for JCW.


		 \param [in,out]	pMesh	The mesh.
		 */

		void JCW_Generate_Regions_On_Base_Mesh(HalfedgeGraph &pMesh);

        /**
         \brief	JCW region mass center insert watermark.

         \param [in,out]	pMesh				   	The mesh.
         \param [in,out]	FP_Geometry			   	The fp geometry.
         \param [in,out]	FP_Region_Number	   	The fp region number.
         \param [in,out]	SP_Watermarked_Position	The sp watermarked position.
         \param [in,out]	SP_Moved_Position	   	The sp moved position.
         \param [in,out]	SP_Original_Position   	The sp original position.
         \param [in,out]	JCW_ERROR			   	The jcw error.
         */

        void JCW_Region_Mass_Center_Insert_Watermark(HalfedgeGraph & pMesh, list<Point3d> & FP_Geometry, list<int> & FP_Region_Number, list<Point_Int> & SP_Watermarked_Position, list<Point3d> & SP_Moved_Position, list<Point3d> & SP_Original_Position, list<vector<int> > & JCW_ERROR);


		/**
		 \brief	JCW region mass center extract watermark.

		 \param [in,out]	pMesh	The mesh.

		 */

		void JCW_Region_Mass_Center_Extract_Watermark(HalfedgeGraph & pMesh);

		/**
		 \brief	JCW code difference histogram shifting.

		 \param [in,out]	pMesh	The mesh.
		 \param	Component_ID	 	Component ID.
		 */

		void JCW_Code_Difference_Histogram_Shifting(HalfedgeGraph &pMesh,const int & Component_ID);

		/**
		 \brief	JCW decode difference histogram shifting.


		 \param [in,out]	pMesh	The mesh.
		 \param	Component_ID	 	Component ID.
		 */

		void JCW_Decode_Difference_Histogram_Shifting(HalfedgeGraph &pMesh, const int & Component_ID);

		/**
		 \brief	Initializes the spherical coordinates.

		 \param [in,out]	pMesh	The mesh.
		 */

		void Initialize_Spherical_Coordinates(HalfedgeGraph &pMesh);

		/**
		 \brief	Converts a cartesian coordinates into spherical coordinates.

		 \param	Pt				   	The point.
		 \param [in,out]	Spheric	If non-null, the spheric.
		 */

		void Convert_To_Spherical(const Point3d & Pt, double * Spheric);

		/**
		 \brief	Converts spherical coordinates into cartesian coordinates.


		 \param	Spheric				 	The spheric.
		 \param [in,out]	Cartesian	If non-null, the cartesian.
		 */

		void Convert_To_Cartesian(const double * Spheric, double * Cartesian);

		/**
		 \brief	JCW decompress one level.

		 \param [in,out]	pMesh	The mesh.
		 \param	File_Name		 	Filename of the file.
		 \param	Noise_mode		 	The noise mode.

		 \return The current level.
		 */

		int  JCW_Decompress_One_Level(HalfedgeGraph &pMesh, const char* File_Name, const int & Noise_mode);

		/**
		 \brief	Jcw decompress one level without extraction.

		 \param [in,out]	pMesh	The mesh.
		 \param	File_Name		 	Filename of the file.

		 \return	The current level.
		 */

		int  JCW_Decompress_One_Level_Without_Extraction(HalfedgeGraph &pMesh, const char* File_Name);


		/**
		 \brief	Jcw divide big regions.


		 \param [in,out]	pMesh	The mesh.
		 \param Thres_divide_regions The threshold.
		 \return Number of divided regions.
		 */

		int  JCW_Divide_Big_Regions(HalfedgeGraph &pMesh, const int & Thres_divide_regions);

		/**
		 \brief	Colorify regions for JCW.

		 \param [in,out]	pMesh	The mesh.
		 */

		void JCW_Colorify_Regions(HalfedgeGraph & pMesh);

		/**
		 \brief	Reads the watermarked information.
		*/

		void    Read_Information_To_Hide(const char * Embedding_message);

		/**
		 \brief	Writes the information to hide.

		 \return Inserted to watermarking message.
		 */

		QString Write_Information_To_Hide();

		/**
		 \brief	Clears after the compression.

		 */

		void Clear_After_Compression();
#endif

private:
  std::vector< bool >
      IsClosed; ///< The is closed. To know if the mesh is open or closed.

  bool IsColored;  ///< true if is colored
  bool IsOneColor; ///< true if is one color

  float OnlyColor[3]; ///< The coordinates of color when there is only one color
  // Number of connectivity symbol types. Without boundary = 5, with boundary =
  // 7;
  // int NummberConnectivitySymbols;

  // To encode each type of operation between decimation and increase of
  // quantization resolution
  std::vector< std::list< int > > ListOperation; ///< The list of operation
  std::vector< std::list< int > >
      Connectivity; ///< The information of connectivity to compress
  std::vector< std::list< Point_Int > >
      Geometry; ///< The geometry information to compress
  std::vector< std::list< int > >
      NumberSymbol; ///< Number of symbols of each simplification step
  std::vector< std::list< int > >
      NumberVertices; ///< Number of vertices of each simplification step

  std::list< Point_Int >
      InterGeometry; ///< The intermediate information of geometry
  std::list< int >
      InterConnectivity; ///< The intermediate information of connectivity

  std::vector< std::list< int > >
      AlphaRange; ///< The range of alpha (Frenet coordinates) of each LoD
  std::vector< std::list< int > > AlphaOffset; ///< The offset of alpha
  std::vector< std::list< int > >
      GammaRange; ///< The range of gamma (Frenet coordinates) of each LoD
  std::vector< std::list< int > > GammaOffset; ///< The offset of gamma


  // Quantization
  std::vector< unsigned > Qbit;           ///< The Quantization bits
  std::vector< float > xmin;              ///< The xmin
  std::vector< float > ymin;              ///< The ymin
  std::vector< float > zmin;              ///< The zmin
  std::vector< float > xmax;              ///< The xmax
  std::vector< float > ymax;              ///< The ymax
  std::vector< float > zmax;              ///< The zmax
  std::vector< float > Quantization_Step; ///< The quantization step

  int Smallest_Alpha; ///< The smallest alpha
  int Smallest_Gamma; ///< The smallest gamma

  std::vector< double > HighestLengthBB; ///< The highest length bb
  std::vector< double > ComponentVolume; ///< The volume of each component
  std::vector< double > ComponentArea;   ///< The area of each component
  std::vector< int >
      ComponentNumberVertices; ///< The number of vertices of each components

  // Used for adatative quantization.
  std::vector< std::list< int > >
      QuantizationCorrectVector; ///< The quantization correct vector
  std::vector< std::list< int > >
      NumberQuantizationLayer; ///< Number of quantization layers

  // for color
  std::vector< std::list< int > >
      NumberProcessedVertices; ///< Number of processed vertices
  std::vector< std::list< int > >
      ColorChildcellIndex; ///< the color childcell index
  std::vector< std::list< int > >
      ColorEncoderIndex; ///< the color encoder index

  // Colors
  std::vector< std::list< Color_Unit > >
      VertexColor; ///< Contain color error of all removed vertices
  std::list< Color_Unit >
      InterVertexColor; ///< The intermediate information vertex color

  float C0_Min; ///< The C0 minimum
  float C1_Min; ///< The C1 minimum
  float C2_Min; ///< The C2 minimum

  float Color_Quantization_Step;              ///< The color quantization step
  std::vector< int > NumberColorQuantization; ///< Number of color quantizations

  int Smallest_C0; ///< The smallest value of C0 used for preventing from
                   ///< negative sylbols
  int Smallest_C1; ///< The smallest value of C1 used for preventing from
                   ///< negative sylbols
  int Smallest_C2; ///< The smallest value of C2 used for preventing from
                   ///< negative sylbols

  int C0_Range; ///< The C0 range
  int C1_Range; ///< The C1 range
  int C2_Range; ///< The C2 range

  /*int ColorDiffMinC0;///< The color difference minimum c0 for restoration of
  colors in Mapping table int ColorDiffMinC1;///< The first color difference
  minimum c1 for restoration of colors in Mapping table int ColorDiffMinC2;///<
  The second color difference minimum c2 for restoration of colors in Mapping
  table

  int ColorDiffRangeC0;///< The color difference range c 0 for restoration of
  colors in Mapping table int ColorDiffRangeC1;///< The first color difference
  range c1 for restoration of colors in Mapping table
  int ColorDiffRangeC2;///< The second color difference range c2 for restoration
  of colors in Mapping table*/

  //// mapping table
  // vector<Color_Unit> ColorArray; ///< Color table
  // vector<Color_Unit> PredictedColorArray; ///< Predicted values of colors in
  // color table using prediction based on neighbors. vector<Color_Unit>
  // DifferenceColor; ///< Difference between original and quantized vertex
  // color. need to turn into lossless coding.
  //
  // list<int>		   ColorIndex;///< List of color index
  // vector<int>		   ReorderingColorIndex;///< Vector of reordering color
  // index list<int>		   InterColorIndex;		///< Vector of
  // inter color index
  // vector<int>		   Number_color_index; ///< contain occurence of each
  // initial color vector<int>		   IsKnownIndex;///< is known index


  // Decoder

  Arithmetic_Codec Decoder; ///< The arithmetic decoder

  Adaptive_Data_Model Color_0_Model; ///< The statistical model used for color_0
  Adaptive_Data_Model Color_1_Model; ///< The statistical model used for color_1
  Adaptive_Data_Model Color_2_Model; ///< The statistical model used for color_2

  Adaptive_Data_Model Index_Model; ///< The index model

  std::vector< int >
      NumberDecimation; ///< To stock number of Decimation operation.
  std::vector< int > NumberChangeQuantization; ///< to stock number of
                                               ///< diminution of quantization.

  int DumpSymbolDecimation; ///< The dump symbol decimation in order to
                            ///< eliminate symbols of the last conquest if it is
                            ///< useless.
  int DumpSymbolRegulation; ///< The dump symbol regulation in order to
                            ///< eliminate symbols of the last conquest if it is
                            ///< useless.

  int Decompress_count; ///< Number of decompress in order to know how many
                        ///< steps to go for the decompression step by step.
  int NumberComponents; ///< Number of components in order to know how many
                        ///< steps to go for the decompression step by step.

  int GlobalCountOperation; ///< The global count operation

  std::vector< int > ComponentOperations; ///< The operations of each component

  // JCW
  int m_NumberBin;         ///< Number of bins for JCW
  int m_EmbeddingStrength; ///< Embedding strength for JCW
  int m_NumberRegion;      ///< Number of regions for JCW
  double m_VC[3];          ///< Mesh center position for JCW
  double m_Rmin; ///< Distance of farthest vertex from mesh center for JCW
  double m_Rmax; ///< Distance of nearst vertex from mesh center for JCW
  double m_Dist; ///< Distance of each bin
  std::vector< int >
      m_Number_Vertices_Per_Regions; ///< Number of vertices in each region
  std::list< int > m_Watermarks;     ///< watermark

  int N_Inserted_Watermarks; ///< Number of inserted watermarks

  std::vector< int >
      m_N_remained_vertices; ///< Number of remained vertices in each region
  std::vector< int >
      m_N_treated_vertices; ///< Number of treated vertices in each region
  std::vector< double > m_Rad_decision;

  std::list< std::vector< int > >
      m_JCW_Move_Error; ///< Stock difference related to complete reversibility
  std::list< int >
      m_N_Errors; ///< Index of error related to complete reversibility.

  Adaptive_Data_Model DM_JCW_MOVE_ERROR;

  std::list< int > JCW_Connectivity; ///< Stock connectivity information for JCW
  std::list< Point_Int > JCW_Geometry; ///< Stock geometry information for JCW

  // double LUT_CourbureClust[3*256];
  std::vector< std::vector< float > > Region_Color; ///< Color of each region
  int Number_non_reversible_vertices; ///< Number of vertices which violate
                                      ///< complte reversibility
  int Number_Save_Over_bins; ///< Number of empty bins to shift the current bins
  bool Is_Division_Big_Regions_Enabled;   ///< true if "division_big_regions"
                                          ///< option is seleted
  bool Is_Complete_Reversibility_Enabled; ///< true if "complete reversibility"
                                          ///< option is seleted
  bool Is_Bijection_Enabled; ///< true if "bijection" option is seleted
  int Division_Threshold;    ///< Threshold of region division

  // from IHM
public:
  bool IsDecompress;             ///< true if is decompress
  bool IsCompressed;             ///< true if is compressed
  int Current_level;             ///< The current level
  int Total_layer;               ///< The total number of layer
  unsigned Initial_file_size;    ///< Size of the initial file
  unsigned Compressed_file_size; ///< Size of the compressed file
  std::string File_name;         ///< Filename of the file

  std::vector< float >
      Prog; ///< Stock information of progression of decompression in
            ///< pourcentage (100 = total decompression)
  std::vector< float >
      Ratio; ///< Stock information of ration regarding size of the input file
  std::string Message; ///< Message to be shown in the main window

public: // private:
  // bool Afficher; ///<
  bool Sequence; ///< Decompression mode (sequence mode or file mode) for IHM
  bool Possible_change_sequence; ///< To disable the mode change during
                                 ///< decompression for IHM

  int Visu_level;    ///< Level of visualized LoD for IHM
  int Process_level; ///< Level of processed LoD for IHM

  FILE *Dec_Info; ///< File to write information of decompression for IHM
  // int Writing_level; ///< Level of
  std::string
      Dec_File_Info; ///< File name to write decompression information for IHM

  // ELO+
  //--------------------------------------------------------------------
  // Property maps that replace ancient vertex/halfedge/face attributes
  //--------------------------------------------------------------------

  typename FEVV::Vertex_pmap< HalfedgeGraph, Color_Unit >
      vertex_color_int; // ELO replace pVertex->color_int(...)
  typename FEVV::Vertex_pmap< HalfedgeGraph, int >
      vertex_Seed_Edge; // ELO replace pVertex->Seed_Edge
  typename FEVV::Vertex_pmap< HalfedgeGraph, int >
      vertex_Component_Number; // ELO replace pVertex->Component_Number
  typename FEVV::Vertex_pmap< HalfedgeGraph, int >
      Vertex_Flag; // ELO replace pVert->Vertex_Flag
  typename FEVV::Vertex_pmap< HalfedgeGraph, int >
      Vertex_Number; // ELO replace pVert->Vertex_Number
  typename FEVV::Vertex_pmap< HalfedgeGraph, int >
      Vertex_Sign; // ELO replace pVert->Vertex_Sign
  typename FEVV::Vertex_pmap< HalfedgeGraph, Vector >
      vertex_normal; // ELO replace pVertex->normal()
  typename FEVV::Vertex_pmap< HalfedgeGraph, int >
      vertex_Q_Index; // ELO replace pVertex->Q_Index
  typename FEVV::Vertex_pmap< HalfedgeGraph, int >
      vertex_Region_Number; // ELO replace pVert->Region_Number
  typename FEVV::Vertex_pmap< HalfedgeGraph, int >
      vertex_Removal_Order; // ELO replace pVert->Removal_Order

  typename FEVV::Face_pmap< HalfedgeGraph, int >
      facet_tag; // ELO replace pFacet->tag(...)
  typename FEVV::Face_pmap< HalfedgeGraph, int >
      facet_Component_Number; // ELO replace pFacet->Component_Number
  typename FEVV::Face_pmap< HalfedgeGraph, int >
      Facet_Flag; // ELO replace pFacet->Facet_Flag
  typename FEVV::Face_pmap< HalfedgeGraph, Vector >
      facet_normal; // ELO replace pFacet->normal()


  //-----------------------------------------------------------
  //    Functions moved from Compression_Valence_Common.h
  //-----------------------------------------------------------

  /**
    \brief	Query if 'h' is border vertex.

    \param	h	The.

    \return	true if border vertex, false if not.
   */
  bool Is_Border_Vertex(const halfedge_descriptor &h,
                        const HalfedgeGraph &mesh);

  /**
   \brief	To find a correspondent type to retriangulate.

   \param	h	   	The.
   \param	valence	The valence.

   \return	The found type.
   */
  int Find_Type(const HalfedgeGraph &_pMesh,
                const halfedge_descriptor &h,
                const int &valence);

  /**
    \brief	Query if removal of this vertex would violate the
    manifold_property or not.

    \param	h	   	The.
    \param	type   	The type.
    \param	valence	The valence.

    \return	true if manifold property violated, false if not.
   */
  bool Is_Manifold_Property_Violated(const HalfedgeGraph &_pMesh,
                                     const halfedge_descriptor &h,
                                     const int &type,
                                     const int &valence);

  /**
    \brief	Query if removal of the front vertex can cause a normal flipping
    problem.

    \param	h	   	The.
    \param	valence	The valence.

    \return	true if normal flipping occured, false if not.
   */
  bool Is_Normal_Flipping_Occured(const HalfedgeGraph &_pMesh,
                                  const PointMap *_pm,
                                  const halfedge_descriptor &h,
                                  const unsigned &valence);

  /**
    \brief	Gives a normal vector of the triangle containing the
    halfedge_handle h.


    \param	h	The.

    \return	.
   */
  Vector Triangle_Normal(const HalfedgeGraph &_pMesh,
                         const PointMap *_pm,
                         const halfedge_descriptor &h);

  /**
    \brief	Gives a normal vector of the triangle formed by three points P Q R
    in the counterclockwise way..


    \param	P	The.
    \param	Q	The.
    \param	R	The.

    \return	.
   */
  Vector Triangle_Normal(const HalfedgeGraph &_pMesh,
                         const Point3d &P,
                         const Point3d &Q,
                         const Point3d &R);

  /**
    \brief	Query if this object is geometric violated. Here, we define a
    geometric_metric to preserve maximum the intermediate meshes.


    \param	h		 	The.
    \param	type	 	The type.
    \param	valence  	The valence.
    \param	Threshold	The threshold.

    \return	true if geometric violated, false if not.
   */
  bool Is_Geometric_Metric_Violated(const HalfedgeGraph &_pMesh,
                                    const PointMap *_pm,
                                    const halfedge_descriptor &h,
                                    const int &type,
                                    const unsigned int &valence,
                                    const float &Threshold);


  /**
    \brief	Gives an area of the triangle which contain the halfedge_handle
    h.

    \param	h	The halfedge_handle.

    \return	.
   */
  // TODO-elo  share this function with Curvature filter that uses its own one
  // TODO-elo  see
  // https://github.com/MEPP-team/MEPP2/blob/master/FEVV/Filters/Generic/Manifold/Curvature/curvature.hpp#L105
  double Area_Facet_Triangle(
      const typename boost::graph_traits< HalfedgeGraph >::halfedge_descriptor
          &h,
      const HalfedgeGraph &mesh,
      const PointMap *pm);

  /**
    \brief	to calculate the area of a triangle.

    \param	P	The.
    \param	Q	The.
    \param	R	The.

    \return	.
   */
  double Area_Facet_Triangle(const Point3d &P,
                             const Point3d &Q,
                             const Point3d &R,
                             const HalfedgeGraph &_pMesh);

  /**
     Compute the volume of ABCD tetrahedron.
     Replacement of CGAL::volume(A, B, C, D).
  */
  double Volume(const Point3d &A,
                const Point3d &B,
                const Point3d &C,
                const Point3d &D,
                const HalfedgeGraph &_pMesh);

  /**
    \brief	gets a color of front vertex of h.


    \param	h	The.

    \return	The vertex color.
   */
  Color_Unit Get_Vertex_Color(const halfedge_descriptor &h,
                              const HalfedgeGraph &_pMesh);

  /**
    \brief	Gets an average vertex color prediction using Lee's method.


    \param	h	   	The.
    \param	valence	The valence.

    \return	The average vertex color lee.
   */
  Color_Unit Get_Average_Vertex_Color_Lee(const halfedge_descriptor &h,
                                          const int &valence,
                                          const HalfedgeGraph &_pMesh);

  /**
    \brief	gives the position of the barycenter of the patch for regulation
    conquest.


    \param	h	The.

    \return	.
   */
  Point3d Barycenter_Patch_Before_Removal(const halfedge_descriptor &h,
                                          const HalfedgeGraph &_pMesh,
                                          const PointMap *_pm);

  /**
    \brief	gives the position of the barycenter of the patch for decimation
    conquest.

    \param	h	   	The.
    \param	valence	The valence.

    \return	.
   */
  Point3d Barycenter_Patch_After_Removal(const halfedge_descriptor &h,
                                         const int &valence,
                                         const HalfedgeGraph &_pMesh,
                                         const PointMap *_pm);

  /**
    \brief	Retriangulates the hole left by a removal of a vertex.


    \param [in,out]	pMesh	The mesh.
    \param	ch				 	The ch.
    \param	valence			 	The valence.
    \param	Vertex_number	 	The vertex number.
    \param	Component_ID	 	Identifier for the component.
   */
  void Retriangulation(HalfedgeGraph &pMesh,
                       const halfedge_descriptor &ch,
                       const unsigned &valence,
                       const unsigned &Vertex_number,
                       const int &Component_ID,
                       const PointMap *_pm);

  /**
    \brief	calculates a normal vector of a patch caused by a removal of a
    front vertex.


    \param	const_h	The constant h.
    \param	valence	The valence.

    \return	.
   */
  Vector Normal_Patch(const halfedge_descriptor &const_h,
                      const unsigned int &valence,
                      const HalfedgeGraph &_pMesh,
                      const PointMap *_pm);

  /**
    \brief	Calculates the base vectors of new coordinates system which is
    frenet system.

    \param	h			  	The.
    \param	normal		  	The normal.
    \param [in,out]	T2	The second Vector &.

    \return	The calculated t 1 t 2.
   */
  Vector Calculate_T1_T2(const halfedge_descriptor &h,
                         const Vector &normal,
                         Vector &T2,
                         const HalfedgeGraph &_pMesh,
                         const PointMap *_pm);

  /**
    \brief	 finds a bijection through a rotation transformation in 3D with
    only integer coordinates.

    \param	Dist  	The distance.
    \param	T1	  	The first const Vector &.
    \param	T2	  	The second const Vector &.
    \param	normal	The normal.

    \return	.
   */
  Point_Int Frenet_Rotation(const Point_Int &Dist,
                            const Vector &T1,
                            const Vector &T2,
                            const Vector &normal);


  /**
    \brief	Inverse operation of frenet rotation. This permits to refind the
    original coordinates.

    \param	Frenet	The frenet.
    \param	T1	  	The first const Vector &.
    \param	T2	  	The second const Vector &.
    \param	normal	The normal.

    \return	.
   */
  Point_Int Inverse_Frenet_Rotation(const Point_Int &Frenet,
                                    const Vector &T1,
                                    const Vector &T2,
                                    const Vector &normal);

  /**
    \brief	Gets an average color of neighboring vertices ( After removal of
    front vertex of h)

    \param	h	   	The.
    \param	valence	The valence.

    \return	The average vertex color after removal.
   */
  Color_Unit
  Get_Average_Vertex_Color_After_Removal(const halfedge_descriptor &h,
                                         const int &valence,
                                         const HalfedgeGraph &_pMesh);

  /**
    \brief	Estimate geometry quantization.

    \param [in,out]	pMesh	The mesh.
    \param	volume			 	The volume.
    \param	area			 	The area.
    \param	number_vertices  	Number of vertices.

    \return	.
   */
  int Estimate_Geometry_Quantization(double volume,
                                     double area,
                                     int number_vertices);

  /**
    \brief	ADAPTIVE_QUANTIZATION : gets symbols to correct Vector of
    under_quantization.


    \param	i	The.
    \param	j	The.
    \param	k	The.

    \return	The correct vector.
   */
  int Get_Correct_Vector(int i, int j, int k);

  /**
    \brief	 Remove edges to create a hole.


    \param [in,out]	pMesh	The mesh.
    \param	h				 	The.
    \param	type			 	The type.

    \return	true if it succeeds, false if it fails.
   */
  bool Remove_Edges(HalfedgeGraph &_pMesh,
                    const halfedge_descriptor &h,
                    const int &type);

  /**
    \brief	Gets a coefficient up quantization.


    \param	Correct_symbol	The correct symbol.
    \param	coeff		  	The coeff.
   */
  void Get_Coefficient_Up_Quantization(const int &Correct_symbol, int coeff[3]);


  //-----------------------------------------------------------
  //                     Miscellaneous
  //-----------------------------------------------------------


  /**
          Compute face and vertex normals.
   */
  void compute_normals(const HalfedgeGraph &_pMesh, const PointMap *_pm);

  /**
          Compute face and vertex normals.
   */
  void compute_normals_per_facet(const HalfedgeGraph &_pMesh,
                                 const PointMap *_pm);

  /**
          Compute face and vertex normals.
   */
  void compute_normals_per_vertex(const HalfedgeGraph &_pMesh);

  /**
          Compute the normal of a face.
   */
  void compute_facet_normal(const face_descriptor &f,
                            const HalfedgeGraph &_pMesh,
                            const PointMap *_pm);

  /**
          Compute the normal of a vertex.
   */
  void compute_vertex_normal(const vertex_descriptor &v,
                             const HalfedgeGraph &_pMesh);

  /**
          Display halfedge points coordinates for debugging purpose.
   */
  void print_halfedge(const std::string &title,
                      const halfedge_descriptor &h,
                      const HalfedgeGraph &_pMesh,
                      const PointMap *_pm);

  /**
          Display vertex point coordinates for debugging purpose.
   */
  void print_vertex(const std::string &title,
                    const vertex_descriptor &v,
                    const PointMap *_pm);

  /**
          Convert vertex to string for debugging purpose.
   */
  std::string vertex_to_string(const vertex_descriptor &v, const PointMap *_pm);

  /**
          Convert halfedge to string for debugging purpose.
   */
  std::string halfedge_to_string(const halfedge_descriptor &h,
                                 const HalfedgeGraph &_pMesh,
                                 const PointMap *_pm);

  /**
          Convert edge to string for debugging purpose.
   */
  std::string edge_to_string(const halfedge_descriptor &h,
                             const HalfedgeGraph &_pMesh,
                             const PointMap *_pm);

  /**
          Display mesh geometry for debugging purpose.
   */
  void DBG_print_mesh_geometry(const HalfedgeGraph &_pMesh,
                               const PointMap *_pm,
                               const std::string &header = std::string());

  /**
          Display mesh vertices color for debugging purpose.
   */
  void DBG_print_mesh_vertexcolor(const HalfedgeGraph &_pMesh,
                                  const VertexColorMap *_v_cm,
                                  const std::string &header = std::string());

  /**
          V1 < V2 ?
   */
  bool v_inf_to_v(const vertex_descriptor &v1,
                  const vertex_descriptor &v2,
                  const PointMap *_pm);

  /**
          Populate the mesh, given vertices and faces list.
  */
  void build_mesh(HalfedgeGraph &_pMesh,
                  PointMap *_pm,
                  VertexColorMap *_v_cm,
                  const std::vector< Point3d > &vlist,
                  const std::vector< int > &flist,
                  const std::vector< float > &clist,
                  const std::vector< int > &Color_index_list);

  /**
          Initialize vertex attributes in the same way as Mepp1.
  */
  void init_vertex_attributes(const vertex_descriptor &v);

  /**
          Truncate vertices colors to 1.
  */
  void truncate_colors(const HalfedgeGraph &_pMesh,
                      VertexColorMap *_v_cm);
};

// implementation details
#include "Compression_Valence_Component.inl"


/*! \page component_Compression_Valence Documentation
 *
 * \section auth Authors
 * H. LEE, C. DIKICI, G. LAVOUE and F. DUPONT \n
 * M2DisCo Team, LIRIS, CNRS, Universite Lyon 1/INSA-Lyon, France. \n
 * \n
 * Please contact the authors \n H. LEE (hosmail@gmail.com), C. DIKICI
 (cagataydikici@yahoo.com),
 * G. LAVOUE (glavoue@liris.cnrs.fr) and F. DUPONT (fdupont@liris.cnrs.fr) \n
 * in case you have any question, comment or a bug to report.
 * \n
 * \n
 *
 * \section paper_reference Related publications
 * Here is the list of related papers if you want to cite our methods. \n
 * It is also strongly recommended that you read these papers
 * in order to know the parameters used in our methods.
 * \n
 *
 * \subsection compression_paper Progressive compression of colored meshes
 * H. LEE, G. LAVOUE, F. DUPONT, \n
 * Rate-distortion optimization for progressive compression of 3D mesh with
 color attributes, \n
 * The Visual Computer, 2011.
 * \n
 *
 * \subsection jcw_paper Joint watermarking and progressive compression of 3D
 meshes
 * H. LEE, C. DIKICI, G. LAVOUE, F. DUPONT, \n
 * Joint reversible watermarking and progressive compression of 3D meshes, \n
 * The Visual Computer 27(6-8):781-792, 2011.
 * \n
 * \n
 *
 * \section overview Overview
 * This Compression Valence component has two main functions:
 * (1) a progressive compresion and (2) a joint progressive compression and
 reversible watermarking for 3D meshes. \n
 * The progressive compression method simplifies iteratively an input mesh to
 generate differents levels of details (LoDs). \n
 * These LoDs are then transmitted progressively in a coarse-to-fine way at the
 decompression. \n
 * In particular, our method adapts the quantization precision (both geometry
 and color) to each LoD in order to optimize the rate-distortion (R-D)
 trade-off. \n \n
 * Our joint progressive compression and reversible watermarking method offers a
 possibility to embed an watermark information
 * in order to protect the ownership of the input mesh. \n
 * Hence, at each decompression step, the inserted message can be extracted also
 progressively. \n
 * The embedded watermarks are reversible, meaning that the deformation caused
 by watermarking embedding step can be removed when extracting the watermark.
 * \n
 * \n
 *
 * \section howto_use How to use this component
 * \subsection howto_progressive_compression Progressive compression
 * First of all you have to load a mesh. \n
 * To compress the input mesh, choose "Compression" and a windows appears. \n
 * (1) File name              : to choose a name for the compressed file. Please
 be aware that the file extension should be ".p3d". \n
 * (2) Mode                   : the mode "Simplification" does not generate a
 compressed file. \n
 * (3) Compression mode       : to enable or disable the use of the adaptation
 of quantization precision. \n
 * (4) Quantization bits      : the number of bits for the geometry
 quantization. \n
 * (5) # Vertices wanted      : the number of vertices of the base (the
 coarsest) mesh. \n
 * (6) Use bijection          : to enable or disable the use of the bijection
 for the geometry encoding. This bijection reduces the coding rates but it needs
 a longer processing time. \n
 * (7) Forbid normal flipping : this option is used to forbid a normal flipping
 when removing a vertex. \n
 * (8) Use Metric             : this option is used to forbid a vertex removal
 if it induces a significant deformation. The threshold value is initially set
 to 0.25. \n The use of this metric can be "forgetted" if the number of vertices
 of the current intermediate mesh is superior to an user-defined threshold. \n
 * \n
 * For the decompression, first you have to load a .p3d file. Then, the base
 mesh is rendered. \n
 * To obtain higher LoDs, you can use : \n
 * (1) "Decompression : all LoDs", to decompress all LoDs and the finest
 intermediate is visualized, \n
 * (2) "Decompression : next LoD", to decompress one level and the next LoD is
 rendered (ALT + left mouse button), \n
 * (3) "Decompression : go to specific LoD", to reach the desired LoD. \n
 * You can also visualize the previous LoDs by selecting "Decompression :
 previous LoD" (ALT + right mouse button). \n \n
 * After loading the .p3d file, you can enable or disable the option of mesh
 sequence generation with the menu "Activate/Desactivate mesh sequence
 generation". \n
 * When this option is enabled, all LoDs are stored in the memory, so that the
 navigation of differents LoDs can be performed quickly. \n
 * You can disable this option in order to save the memory.
 * In this case, when you want to visualize the previous LoD, the decompression
 is performed again until getting the previous LoD.
 * The navigation takes a longer time. \n
 * Note that this option can be modified only after the rendering of the base
 mesh and it cannot be modified when any decompression operation is achieved. \n
 * \n
 * \subsection howto_jcw Joint compression and watermarking
 * First of all you have to load a mesh. \n
 * You can apply our joint method by selecting "JCW - Compression and
 Watermarking embedding". \n
 * (1) File name              : to choose a name for the compressed file (.p3d).
 \n
 * (2) Q bits                 : the number of bits for the geometry
 quantization. \n
 * (3) # vertices             : the number of the base mesh after an iterative
 simplification. \n
 * (4) # Bins				  : the number of the histogram bins. \n
 * (5) # Regions              : the number of regions. One watermark bit is
 embedded in each region. \n
 * (6) Embedding Strength     : the number of shifted bins when
 embedding/extracting the watermark. \n
 * (7) Embedding Message      : the message to insert. When this field is empty
 or the length of the message is shorter than necessary, the message is
 generated randomly. \n
 * (8) Divide Regions         : when this option is selected, the big region is
 divided in two in order to insert more watermark bits. \n
 * (9) Complete Reversibility : When this option is checked, the initial
 positions of all vertices are exactly restored. Some extra coding bits are
 necessary.
 * \n
 * For the decompression and the watermark extraction,
 * you have to load a .p3d file. \n
 * To obtain higher LoDs, you can select "JCW - Decompression and Watermark
 extraction : next LoD". \n
 * The extracted message is shown in the status bar of the main window. \n
 * To visualize the results without watermark extraction (non authorized users),
 you can use "JCW - Decompression without extraction : next LoD". \n
 * \n
 *
 * \section last_updated Last updated
 * 18 May, 2011
 *
 */
