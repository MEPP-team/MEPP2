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

#include <iostream>
#include <exception>

#include "FEVV/Tools/IO/FileUtilities.hpp"
#include "FEVV/Tools/IO/StringUtilities.hpp"

#if defined(__APPLE__) && defined(__clang__)
// In order to silence out warnings about VTK deprecated code invocation
#define VTK_LEGACY_SILENT
// We also want to silence out vtk complains about member function not being
// marked override:
#pragma clang diagnostic push
#pragma clang diagnostic ignored "-Winconsistent-missing-override"
#endif
#include <vtkSmartPointer.h>
#include <vtkPointData.h>
#include <vtkCellData.h>
#include <vtkCellArray.h>
#include <vtkPolyData.h>
#include <vtkUnstructuredGrid.h>
#include <vtkDoubleArray.h>
#include <vtkFieldData.h>

#include <vtkPolyDataWriter.h>         // for vtk-files
#include <vtkXMLPolyDataWriter.h>      // for vtp-files
#include <vtkUnstructuredGridWriter.h> // for vtu-files
#include <vtkXMLUnstructuredGridWriter.h>
#if defined(__APPLE__) && defined(__clang__)
#pragma clang diagnostic pop
#endif

namespace FEVV {
namespace IO {
template< typename CoordType,
          typename CoordNType,
          /*typename coordT_type,*/ typename CoordCType,
          typename IndexType >
void
load_vtk_poly_data(
    const std::vector< std::vector< CoordType > > &points_coords,
    const std::vector< std::vector< CoordNType > > &normals_coords,
    const std::vector< std::vector< CoordCType > > &vertex_color_coords,
    const std::vector< std::vector< IndexType > >
        &line_indices, /// edges (can also be used to decribe a volume element
                       /// [old compatibility with Xavier Faure work])
    const std::vector< std::vector< CoordCType > > &lines_color_coords,
    const std::vector< std::vector< IndexType > >
        &face_indices, /// polygonal faces made of vertex indices
    const std::vector< std::vector< CoordCType > > &face_color_coords,
    const std::vector< std::vector< std::vector< double > > >
        &field_attributes, // array of 2D arrays
    const std::vector< std::string > &field_names,
    vtkSmartPointer< vtkPolyData > &poly_data)
{
  /////////////////////////////////////////////////////////////////////////////
  size_t nb_points = points_coords.size(), nb_normals = normals_coords.size(),
         nb_point_colors = vertex_color_coords.size(),
         nb_lines = line_indices.size(),
         nb_line_colors = lines_color_coords.size(),
         nb_faces = face_indices.size(),
         nb_face_colors = face_color_coords.size(),
         nb_arrays_in_field = field_attributes.size();
  bool has_one_name_for_each_field = (nb_arrays_in_field == field_names.size());
  /////////////////////////////////////////////////////////////////////////////
  // some checks
  assert((nb_normals == 0u) || (nb_points == nb_normals));
  /////////////////////////////////////////////////////////////////////////////
  poly_data = vtkSmartPointer< vtkPolyData >::New();
  /////////////////////////////////////////////////////////////////////////////
  // add points
  // std::cout << " before adding points..." << std::endl ;
  vtkSmartPointer< vtkPoints > ptr_points = vtkSmartPointer< vtkPoints >::New();
  for(size_t i = 0; i < nb_points;
      ++i) // vtk writer will align point coordinates using 9 columns (so 3 3D
           // points for the first line etc.)
  {
    // ptrPoints->InsertPoint( static_cast<vtkIdType>(i), points_coords[i][0],
    // points_coords[i][1], points_coords[i][2] );
    ptr_points->InsertNextPoint(
        points_coords[i][0], points_coords[i][1], points_coords[i][2]);
  }
  // ptrPoints->Modified();
  poly_data->SetPoints(ptr_points);
  // std::cout << " after adding points..." << std::endl ;
  /////////////////////////////////////////////////////////////////////////////
  // add point normals
  // std::cout << " before adding point normals..." << std::endl ;
  if((nb_normals > 0u) &&
     (nb_points ==
      nb_normals)) // vtk writer will align normal coordinates using 9 columns
                   // (so 3 3D normals for the first line etc.)
  {
    vtkSmartPointer< vtkDoubleArray > ptr_normals =
        vtkSmartPointer< vtkDoubleArray >::New();
    ptr_normals->SetName("vertex_normals"); // set the point normal name (the
                                            // default name is "normals")
    ptr_normals->SetNumberOfComponents(3);  // 3d normals (ie x,y,z)
    ptr_normals->SetNumberOfTuples(poly_data->GetNumberOfPoints());

    for(size_t i = 0; i < nb_normals; ++i)
    {
      ptr_normals->SetTuple3(
          i, normals_coords[i][0], normals_coords[i][1], normals_coords[i][2]);
      // ptrNormals->InsertTuple3(i, normals_coords[i][0], normals_coords[i][1],
      // normals_coords[i][2] ) ; ptrNormals->InsertNextTuple3(
      // normals_coords[i][0], normals_coords[i][1], normals_coords[i][2] );
    }

    poly_data->GetPointData()->SetNormals(ptr_normals);
    // std::cout << " nb of points = " << poly_data->GetNumberOfPoints() ;
    // std::cout << " ; nb of point normals = " <<
    // poly_data->GetPointData()->GetNormals()->GetNumberOfTuples() ;
    assert(poly_data->GetNumberOfPoints() ==
           poly_data->GetPointData()->GetNormals()->GetNumberOfTuples());
  }
  // std::cout << " after adding point normals..." << std::endl ;
  /////////////////////////////////////////////////////////////////////////////
  // add point colors
  if((nb_point_colors > 0u) && (vertex_color_coords[0].size() > 2u) &&
     (nb_points == nb_point_colors))
  {
    vtkSmartPointer< vtkDoubleArray > ptr_colors =
        vtkSmartPointer< vtkDoubleArray >::New();
    ptr_colors->SetName("vertex_colors"); // set the point color name
    ptr_colors->SetNumberOfComponents(
        static_cast< int >(vertex_color_coords[0].size())); // 3d or 4d
    ptr_colors->SetNumberOfTuples(nb_point_colors);

    for(size_t i = 0; i < nb_point_colors; ++i)
    {
      if(vertex_color_coords[0].size() == 3)
        ptr_colors->SetTuple3(i,
                              vertex_color_coords[i][0],
                              vertex_color_coords[i][1],
                              vertex_color_coords[i][2]);
      else
        ptr_colors->SetTuple4(i,
                              vertex_color_coords[i][0],
                              vertex_color_coords[i][1],
                              vertex_color_coords[i][2],
                              vertex_color_coords[i][3]);
    }

    poly_data->GetPointData()->AddArray(ptr_colors);
  }
  // std::cout << " after adding point colors..." << std::endl ;
  /////////////////////////////////////////////////////////////////////////////
  // add polygons
  // std::cout << " before adding polygons..." << std::endl ;
  vtkSmartPointer< vtkCellArray > polys =
      vtkSmartPointer< vtkCellArray >::New();
  for(size_t i = 0; i < nb_faces; ++i)
  {
    vtkIdType *pt_indices = new vtkIdType[face_indices[i].size()];
    typename std::vector< IndexType >::const_iterator it_begin(
        face_indices[i].begin()),
        it_end(face_indices[i].end());
    unsigned long cpt = 0;
    while(it_begin != it_end)
    {
      pt_indices[cpt++] = static_cast< vtkIdType >(*it_begin);

      ++it_begin;
    }
    polys->InsertNextCell(static_cast< int >(face_indices[i].size()),
                          pt_indices);

    delete[] pt_indices;
  }
  poly_data->SetPolys(polys);
  // std::cout << " after adding polygons..." << std::endl ;
  /////////////////////////////////////////////////////////////////////////////
  // add polygon colors
  if((nb_face_colors > 0u) && (face_color_coords[0].size() > 2u) &&
     (nb_faces == nb_face_colors))
  {
    vtkSmartPointer< vtkDoubleArray > ptr_colors =
        vtkSmartPointer< vtkDoubleArray >::New();
    ptr_colors->SetName("face_colors"); // set the point color name
    ptr_colors->SetNumberOfComponents(
        static_cast< int >(face_color_coords[0].size())); // 3d or 4d
    ptr_colors->SetNumberOfTuples(nb_face_colors);

    for(size_t i = 0; i < nb_face_colors; ++i)
    {
      if(face_color_coords[0].size() == 3)
        ptr_colors->SetTuple3(i,
                              face_color_coords[i][0],
                              face_color_coords[i][1],
                              face_color_coords[i][2]);
      else
        ptr_colors->SetTuple4(i,
                              face_color_coords[i][0],
                              face_color_coords[i][1],
                              face_color_coords[i][2],
                              face_color_coords[i][3]);
    }
    if(nb_lines == 0)
      poly_data->GetCellData()->AddArray(ptr_colors);
    else
      poly_data->GetFieldData()->AddArray(ptr_colors);
  }
  // std::cout << " after adding face colors..." << std::endl ;
  /////////////////////////////////////////////////////////////////////////////
  // add lines
  vtkSmartPointer< vtkCellArray > lines_of_volume_cell =
      vtkSmartPointer< vtkCellArray >::New();
  for(size_t i = 0; i < nb_lines; ++i)
  {
    vtkIdType *pt_indices = new vtkIdType[line_indices[i].size()];
    typename std::vector< IndexType >::const_iterator it_begin(
        line_indices[i].begin()),
        it_end(line_indices[i].end());
    unsigned long cpt = 0;
    while(it_begin != it_end)
    {
      pt_indices[cpt++] = static_cast< vtkIdType >(*it_begin);

      ++it_begin;
    }
    lines_of_volume_cell->InsertNextCell(
        static_cast< int >(line_indices[i].size()), pt_indices);

    delete[] pt_indices;
  }
  poly_data->SetLines(lines_of_volume_cell);
  /////////////////////////////////////////////////////////////////////////////////////////////
  // add lines colors
  if((nb_line_colors > 0u) && (lines_color_coords[0].size() > 2u) &&
     (nb_lines == nb_line_colors))
  {
    vtkSmartPointer< vtkDoubleArray > ptr_colors =
        vtkSmartPointer< vtkDoubleArray >::New();
    if(line_indices[0].size() > 2)
      ptr_colors->SetName("cell_colors");
    else
      ptr_colors->SetName("edge_colors");

    ptr_colors->SetNumberOfComponents(
        static_cast< int >(lines_color_coords[0].size())); // 3d or 4d
    ptr_colors->SetNumberOfTuples(nb_line_colors);

    for(size_t i = 0; i < nb_line_colors; ++i)
    {
      if(lines_color_coords[0].size() == 3)
        ptr_colors->SetTuple3(i,
                              lines_color_coords[i][0],
                              lines_color_coords[i][1],
                              lines_color_coords[i][2]);
      else
        ptr_colors->SetTuple4(i,
                              lines_color_coords[i][0],
                              lines_color_coords[i][1],
                              lines_color_coords[i][2],
                              lines_color_coords[i][3]);
    }
    if(nb_faces == 0)
      poly_data->GetCellData()->AddArray(ptr_colors);
    else
      poly_data->GetFieldData()->AddArray(ptr_colors);
  }
  /////////////////////////////////////////////////////////////////////////////
  // Set 2D arrays of FIELD data when present
  // std::cout << " before adding data fields..." << std::endl ;
  if(nb_arrays_in_field > 0u) // vtk writer will align arrays using 9 columns
  {
    int choice = 0;
    // std::cout << nbArraysInField << " data fields to add..." << std::endl ;
    for(size_t id_array = 0; id_array < nb_arrays_in_field; id_array++)
    {
      // std::cout << " data fields number " << id_array << std::endl ;
      // std::cout << field_names[id_array] << std::endl;
      vtkSmartPointer< vtkDoubleArray > ptr_data =
          vtkSmartPointer< vtkDoubleArray >::New();
      if(has_one_name_for_each_field)
      {
        if(field_names[id_array].find("POINT_DATA_") != std::string::npos)
        {
          ptr_data->SetName((field_names[id_array].substr(11)).c_str());
          choice = 1;
        }
        else if(field_names[id_array].find("CELL_DATA_") != std::string::npos)
        {
          ptr_data->SetName((field_names[id_array].substr(10)).c_str());
          choice = 2;
        }
        else if(field_names[id_array].find("ELEMENT_DATA_") != std::string::npos)
        {
          ptr_data->SetName((field_names[id_array].substr(13)).c_str());
          choice = 2;
        }		
        else
        {
          ptr_data->SetName(field_names[id_array].c_str());
          choice = 0;
        }
      }
      else
      {
        ptr_data->SetName(
            (std::string("cell_law_") + StrUtils::convert(id_array))
                .c_str()); // name of array id_array
        choice = 0;
      }
      ptr_data->SetNumberOfComponents(
          ((field_attributes[id_array].size() > 0u)
               ? static_cast< int >(field_attributes[id_array][0].size())
               : 0));
      ptr_data->SetNumberOfTuples(field_attributes[id_array].size());

      // std::cout << " it has " << ptrData->GetNumberOfTuples() << " tuples and
      // " << ptrData->GetNumberOfComponents() << " components" << std::endl ;
      // Set 2D array values:
      for(vtkIdType i = 0; i < ptr_data->GetNumberOfTuples();
          ++i) // for each line
      {
        assert(i < static_cast< vtkIdType >(field_attributes[id_array].size()));
        for(vtkIdType j = 0; j < ptr_data->GetNumberOfComponents();
            ++j) // for each column
        {
          if((j >=
              static_cast< vtkIdType >(field_attributes[id_array][i].size())) ||
             (ptr_data->GetNumberOfComponents() !=
              static_cast< vtkIdType >(field_attributes[id_array][i].size())))
          {
            std::cerr << " Trying to access index " << j << " but only "
                      << field_attributes[id_array][i].size()
                      << " components in the data field." << std::endl;
            std::cerr << " We should find " << ptr_data->GetNumberOfComponents()
                      << " components in the provided data field. "
                      << std::endl;
            assert(false);
          }
          // if( field_attributes[id_array][i][j]<0. ||
          // field_attributes[id_array][i][j]>1.0 )
          //{
          //	std::cout << "Warning for GenSim project: attributes must be
          //within [0;1] interval for the time being. " << std::endl ;
          //}
          ptr_data->SetComponent(i, j, field_attributes[id_array][i][j]);
        }
      }
      ///////////////////////////////////////////////////////////////////////////
      switch(choice)
      {
      case 1:
        poly_data->GetPointData()->AddArray(ptr_data);
        break;
      case 2:
        poly_data->GetCellData()->AddArray(ptr_data);
        break;
      default:
        poly_data->GetFieldData()->AddArray(ptr_data);
      }
    }
  }
  // std::cout << " after adding data fields..." << std::endl ;
}
template< typename CoordType,
          typename CoordNType,
          /*typename coordT_type,*/ typename CoordCType,
          typename IndexType >
void
load_vtk_unstructured_grid(
    const std::vector< std::vector< CoordType > > &points_coords,
    const std::vector< std::vector< CoordNType > > &normals_coords,
    const std::vector< std::vector< CoordCType > > &vertex_color_coords,
    const std::vector< std::vector< IndexType > >
        &line_indices, /// edges (can also be used to decribe a volume element
                       /// [old compatibility with Xavier Faure work])
    const std::vector< std::vector< CoordCType > > &lines_color_coords,
    const std::vector< std::vector< IndexType > >
        &face_indices, /// polygonal facets made of vertex indices
    const std::vector< std::vector< CoordCType > > &face_color_coords,
    const std::vector< std::vector< std::vector< double > > >
        &field_attributes, // array of 2D arrays
    const std::vector< std::string > &field_names,
    vtkSmartPointer< vtkUnstructuredGrid > &unstructured_grid)
{
  /////////////////////////////////////////////////////////////////////////////
  size_t nb_points = points_coords.size(), nb_normals = normals_coords.size(),
         nb_point_colors = vertex_color_coords.size(),
         nb_lines = line_indices.size(),
         nb_line_colors = lines_color_coords.size(),
         nb_faces = face_indices.size(),
         nb_face_colors = face_color_coords.size(),
         nb_arrays_in_field = field_attributes.size();
  bool has_one_name_for_each_field = (nb_arrays_in_field == field_names.size());
  /////////////////////////////////////////////////////////////////////////////
  // some checks
  assert((nb_normals == 0u) || (nb_points == nb_normals));
  /////////////////////////////////////////////////////////////////////////////
  unstructured_grid = vtkSmartPointer< vtkUnstructuredGrid >::New();
  /////////////////////////////////////////////////////////////////////////////
  // add points
  // std::cout << " before adding points..." << std::endl ;
  vtkSmartPointer< vtkPoints > ptr_points = vtkSmartPointer< vtkPoints >::New();
  for(size_t i = 0; i < nb_points;
      ++i) // vtk writer will align point coordinates using 9 columns (so 3 3D
           // points for the first line etc.)
  {
    // ptrPoints->InsertPoint( static_cast<vtkIdType>(i), points_coords[i][0],
    // points_coords[i][1], points_coords[i][2] );
    ptr_points->InsertNextPoint(
        points_coords[i][0], points_coords[i][1], points_coords[i][2]);
  }
  unstructured_grid->SetPoints(ptr_points);
  // std::cout << " after adding points..." << std::endl ;
  /////////////////////////////////////////////////////////////////////////////
  // add point normals
  // std::cout << " before adding point normals..." << std::endl ;
  if((nb_normals > 0u) &&
     (nb_points ==
      nb_normals)) // vtk writer will align normal coordinates using 9 columns
                   // (so 3 3D normals for the first line etc.)
  {
    vtkSmartPointer< vtkDoubleArray > ptr_normals =
        vtkSmartPointer< vtkDoubleArray >::New();
    ptr_normals->SetName("vertex_normals"); // set the point normal name (the
                                            // default name is "normals")
    ptr_normals->SetNumberOfComponents(3);  // 3d normals (ie x,y,z)
    ptr_normals->SetNumberOfTuples(unstructured_grid->GetNumberOfPoints());

    for(size_t i = 0; i < nb_normals; ++i)
    {
      ptr_normals->SetTuple3(
          i, normals_coords[i][0], normals_coords[i][1], normals_coords[i][2]);
      // ptrNormals->InsertTuple3(i, normals_coords[i][0], normals_coords[i][1],
      // normals_coords[i][2] ) ; ptrNormals->InsertNextTuple3(
      // normals_coords[i][0], normals_coords[i][1], normals_coords[i][2] );
    }

    unstructured_grid->GetPointData()->SetNormals(ptr_normals);
    // std::cout << " nb of points = " << poly_data->GetNumberOfPoints() ;
    // std::cout << " ; nb of point normals = " <<
    // poly_data->GetPointData()->GetNormals()->GetNumberOfTuples() ;
    assert(
        unstructured_grid->GetNumberOfPoints() ==
        unstructured_grid->GetPointData()->GetNormals()->GetNumberOfTuples());
  }
  // std::cout << " after adding point normals..." << std::endl ;
  /////////////////////////////////////////////////////////////////////////////
  // add point colors
  if((nb_point_colors > 0u) && (vertex_color_coords[0].size() > 2u) &&
     (nb_points == nb_point_colors))
  {
    vtkSmartPointer< vtkDoubleArray > ptr_colors =
        vtkSmartPointer< vtkDoubleArray >::New();
    ptr_colors->SetName("vertex_colors"); // set the point color name
    ptr_colors->SetNumberOfComponents(
        static_cast< int >(vertex_color_coords[0].size())); // 3d or 4d
    ptr_colors->SetNumberOfTuples(nb_point_colors);

    for(size_t i = 0; i < nb_point_colors; ++i)
    {
      if(vertex_color_coords[0].size() == 3)
        ptr_colors->SetTuple3(i,
                              vertex_color_coords[i][0],
                              vertex_color_coords[i][1],
                              vertex_color_coords[i][2]);
      else
        ptr_colors->SetTuple4(i,
                              vertex_color_coords[i][0],
                              vertex_color_coords[i][1],
                              vertex_color_coords[i][2],
                              vertex_color_coords[i][3]);
    }

    unstructured_grid->GetPointData()->AddArray(ptr_colors);
  }
  // std::cout << " after adding point colors..." << std::endl ;
  /////////////////////////////////////////////////////////////////////////////
  // add polygons
  // std::cout << " before adding polygons..." << std::endl ;
  int *type_array = NULL;
  if(nb_faces > 0 || nb_lines > 0)
    type_array = new(std::nothrow) int[nb_faces + nb_lines];
  vtkSmartPointer< vtkCellArray > cells =
      vtkSmartPointer< vtkCellArray >::New();
  for(size_t i = 0; i < nb_faces; ++i)
  {
    vtkIdType *pt_indices = new vtkIdType[face_indices[i].size()];
    typename std::vector< IndexType >::const_iterator it_begin(
        face_indices[i].begin()),
        it_end(face_indices[i].end());
    unsigned long cpt = 0;
    while(it_begin != it_end)
    {
      pt_indices[cpt++] = static_cast< vtkIdType >(*it_begin);

      ++it_begin;
    }
    cells->InsertNextCell(static_cast< int >(face_indices[i].size()),
                          pt_indices);

    switch(face_indices[i].size())
    {
    case 0:
    case 1:
    case 2:
      throw std::runtime_error(
          "Writer::load_vtkUnstructuredGrid: a face seems to be not valid.");
    case 3:
      type_array[i] = VTK_TRIANGLE;
      break;
    case 4:
      type_array[i] = VTK_QUAD;
      break;
    default:
      type_array[i] = VTK_POLYGON;
    }
    delete[] pt_indices;
  }
  /////////////////////////////////////////////////////////////////////////////
  // add polygon colors
  if((nb_face_colors > 0u) && (face_color_coords[0].size() > 2u) &&
     (nb_faces == nb_face_colors))
  {
    vtkSmartPointer< vtkDoubleArray > ptr_colors =
        vtkSmartPointer< vtkDoubleArray >::New();
    ptr_colors->SetName("face_colors"); // set the point color name
    ptr_colors->SetNumberOfComponents(
        static_cast< int >(face_color_coords[0].size())); // 3d or 4d
    ptr_colors->SetNumberOfTuples(nb_face_colors);

    for(size_t i = 0; i < nb_face_colors; ++i)
    {
      if(face_color_coords[0].size() == 3)
        ptr_colors->SetTuple3(i,
                              face_color_coords[i][0],
                              face_color_coords[i][1],
                              face_color_coords[i][2]);
      else
        ptr_colors->SetTuple4(i,
                              face_color_coords[i][0],
                              face_color_coords[i][1],
                              face_color_coords[i][2],
                              face_color_coords[i][3]);
    }
    if(nb_lines == 0)
      unstructured_grid->GetCellData()->AddArray(ptr_colors);
    else
      unstructured_grid->GetFieldData()->AddArray(ptr_colors);
  }
  // std::cout << " after adding face colors..." << std::endl ;
  /////////////////////////////////////////////////////////////////////////////
  for(size_t i = 0; i < nb_lines; ++i)
  {
    vtkIdType *pt_indices = new vtkIdType[line_indices[i].size()];
    typename std::vector< IndexType >::const_iterator it_begin(
        line_indices[i].begin()),
        it_end(line_indices[i].end());
    unsigned long cpt = 0;
    while(it_begin != it_end)
    {
      pt_indices[cpt++] = static_cast< vtkIdType >(*it_begin);

      ++it_begin;
    }
    cells->InsertNextCell(static_cast< int >(line_indices[i].size()),
                          pt_indices);
    switch(cpt)
    {
    case 4:
      type_array[nb_faces + i] = VTK_TETRA;
      break;
    case 6:
      type_array[nb_faces + i] = VTK_HEXAHEDRON;
      break;
    default:
      throw std::runtime_error(
          "Writer::load_vtkUnstructuredGrid: a face seems to be not valid.");
    }
    delete[] pt_indices;
  }
  /////////////////////////////////////////////////////////////////////////////
  unstructured_grid->SetCells(type_array, cells);
  if(type_array != NULL)
    delete type_array;
  // std::cout << " after adding cells..." << std::endl ;
  /////////////////////////////////////////////////////////////////////////////////////////////
  // add lines colors
  if((nb_line_colors > 0u) && (lines_color_coords[0].size() > 2u) &&
     (nb_lines == nb_line_colors))
  {
    vtkSmartPointer< vtkDoubleArray > ptr_colors =
        vtkSmartPointer< vtkDoubleArray >::New();
    ptr_colors->SetName("cell_colors");

    ptr_colors->SetNumberOfComponents(
        static_cast< int >(lines_color_coords[0].size())); // 3d or 4d
    ptr_colors->SetNumberOfTuples(nb_line_colors);

    for(size_t i = 0; i < nb_line_colors; ++i)
    {
      if(lines_color_coords[0].size() == 3)
        ptr_colors->SetTuple3(i,
                              lines_color_coords[i][0],
                              lines_color_coords[i][1],
                              lines_color_coords[i][2]);
      else
        ptr_colors->SetTuple4(i,
                              lines_color_coords[i][0],
                              lines_color_coords[i][1],
                              lines_color_coords[i][2],
                              lines_color_coords[i][3]);
    }
    unstructured_grid->GetCellData()->AddArray(ptr_colors);
  }
  // std::cout << " after adding line colors..." << std::endl ;
  /////////////////////////////////////////////////////////////////////////////
  // Set 2D arrays of FIELD data when present
  // std::cout << " before adding data fields..." << std::endl ;
  if(nb_arrays_in_field > 0u) // vtk writer will align arrays using 9 columns
  {
    int choice = 0;
    // std::cout << nbArraysInField << " data fields to add..." << std::endl ;
    for(size_t id_array = 0; id_array < nb_arrays_in_field; id_array++)
    {
      // std::cout << " data fields number " << id_array << std::endl ;
      vtkSmartPointer< vtkDoubleArray > ptr_data =
          vtkSmartPointer< vtkDoubleArray >::New();
      if(has_one_name_for_each_field)
      {
        if(field_names[id_array].find("POINT_DATA_") != std::string::npos)
        {
          ptr_data->SetName((field_names[id_array].substr(11)).c_str());
          choice = 1;
        }
        else if(field_names[id_array].find("CELL_DATA_") != std::string::npos)
        {
          ptr_data->SetName((field_names[id_array].substr(10)).c_str());
          choice = 2;
        }
        else if(field_names[id_array].find("ELEMENT_DATA_") != std::string::npos)
        {
          ptr_data->SetName((field_names[id_array].substr(13)).c_str());
          choice = 2;
        }		
        else
        {
          ptr_data->SetName(field_names[id_array].c_str());
          choice = 0;
        }
      }
      else
      {
        ptr_data->SetName(
            (std::string("cell_law_") + StrUtils::convert(id_array))
                .c_str()); // name of array id_array
        choice = 0;
      }
      ptr_data->SetNumberOfComponents(
          ((field_attributes[id_array].size() > 0u)
               ? static_cast< int >(field_attributes[id_array][0].size())
               : 0));
      ptr_data->SetNumberOfTuples(field_attributes[id_array].size());

      // std::cout << " it has " << ptrData->GetNumberOfTuples() << " tuples and
      // " << ptrData->GetNumberOfComponents() << " components" << std::endl ;
      // Set 2D array values:
      for(vtkIdType i = 0; i < ptr_data->GetNumberOfTuples();
          ++i) // for each line
      {
        assert(i < static_cast< vtkIdType >(field_attributes[id_array].size()));
        for(vtkIdType j = 0; j < ptr_data->GetNumberOfComponents();
            ++j) // for each column
        {
          if((j >=
              static_cast< vtkIdType >(field_attributes[id_array][i].size())) ||
             (ptr_data->GetNumberOfComponents() !=
              static_cast< vtkIdType >(field_attributes[id_array][i].size())))
          {
            std::cout << " Trying to access index " << j << " but only "
                      << field_attributes[id_array][i].size()
                      << " components in the data field." << std::endl;
            std::cout << " We should find " << ptr_data->GetNumberOfComponents()
                      << " components in the provided data field. "
                      << std::endl;
            assert(false);
          }
          // if( field_attributes[id_array][i][j]<0. ||
          // field_attributes[id_array][i][j]>1.0 )
          //{
          //	std::cout << "Warning for GenSim project: attributes must be
          //within [0;1] interval for the time being. " << std::endl ;
          //}
          ptr_data->SetComponent(i, j, field_attributes[id_array][i][j]);
        }
      }
      ///////////////////////////////////////////////////////////////////////////
      switch(choice)
      {
      case 1:
        unstructured_grid->GetPointData()->AddArray(ptr_data);
        break;
      case 2:
        unstructured_grid->GetCellData()->AddArray(ptr_data);
        // unstructured_grid->GetCellData()->SetScalars(ptrData); // only the
        // last ptrData set as scalars will be kept!
        // for the time being we only manage data as fields...
        break;
      default:
        unstructured_grid->GetFieldData()->AddArray(ptr_data);
      }
    }
  }
}

template< typename CoordType,
          typename CoordNType,
          /*typename coordT_type,*/ typename CoordCType,
          typename IndexType >
void
write_vtk_or_vtp_or_vtu_file(
    std::string file_path,
    const std::vector< std::vector< CoordType > > &points_coords,
    const std::vector< std::vector< CoordNType > > &normals_coords,
    const std::vector< std::vector< CoordCType > > &vertex_color_coords,
    const std::vector< std::vector< IndexType > >
        &line_indices, /// edges (can also be used to decribe a volume element
                       /// [old compatibility with Xavier Faure work])
    const std::vector< std::vector< CoordCType > > &lines_color_coords,
    const std::vector< std::vector< IndexType > >
        &face_indices, /// polygonal facets made of vertex indices
    const std::vector< std::vector< CoordCType > > &face_color_coords,
    const std::vector< std::vector< std::vector< double > > >
        &field_attributes, // array of 2D arrays
    const std::vector< std::string > &field_names)
{
  bool is_a_poly_data =
      ((line_indices.size() == 0) ||
       ((line_indices.size() > 0) && (face_indices.size() == 0) &&
        !FileUtils::has_extension(
            file_path,
            ".vtu"))); // isAPolyData is true if there is no line indices [thus
                       // only vertices+faces] or if there is line indices
                       // without any faces not considered as volume [thus needs
                       // to ensure extension is not .vtu]
  vtkSmartPointer< vtkPolyData > poly_data;
  vtkSmartPointer< vtkUnstructuredGrid > unstructured_grid;
  ////////////////////////////////////////////////////////////////////////////  
  std::vector<std::string> names;
  if (field_names.size() != field_attributes.size())
  {
    names.resize(field_attributes.size());
    // When data field names need to be set, we assume that nD (n>=2)
    // fields associated with positions are displacement fields while others
    // are understood as physical laws.
    for (unsigned long i(0); i<field_attributes.size(); ++i){
      if( (points_coords.size() == field_attributes[i].size()) && 
        ((points_coords.size() != face_indices.size()) || (field_attributes[i][0].size()>1)))
           names[i] = std::string("Shifting_") + convert(i);  
        else
           names[i] = std::string("cell_law_")+ convert(i);
    }
  }
  else
    names = field_names;  
  ////////////////////////////////////////////////////////////////////////////
  if(is_a_poly_data)
  {
    load_vtk_poly_data(points_coords,
                       normals_coords,
                       vertex_color_coords,
                       line_indices,
                       lines_color_coords,
                       face_indices,
                       face_color_coords,
                       field_attributes,
                       names,
                       poly_data);
  }
  else
  {
    // there is no face and no face color for unstructured_grid meshes.
    std::vector< std::vector< IndexType > > new_facets;
    std::vector< std::vector< CoordCType > > new_facets_c;
    ////////////////////////////////////////////////////////////////////////////
    load_vtk_unstructured_grid(points_coords,
                               normals_coords,
                               vertex_color_coords,
                               line_indices,
                               lines_color_coords,
                               new_facets,
                               new_facets_c,
                               field_attributes,
                               names,
                               unstructured_grid);
  }
  ////////////////////////////////////////////////////////////////////////////
  if(FileUtils::has_extension(file_path, ".vtp"))
  {
    if(is_a_poly_data)
    {
      vtkSmartPointer< vtkXMLPolyDataWriter > writer =
          vtkSmartPointer< vtkXMLPolyDataWriter >::New();
      writer->SetFileName(file_path.c_str());
#if VTK_MAJOR_VERSION <= 5
      writer->SetInput(poly_data);
      // writer->SetInputConnection(poly_data->GetProducerPort()) ;
#else
      // VTK 6 code
      writer->SetInputData(poly_data);
#endif
      writer->Update();
      writer->Write();
    }
    else
    {
      throw std::runtime_error("Writer::write_vtk_or_vtp_or_vtu_file -> file "
                               "extension .vtp is reserved for PolyData.");
    }
  }
  else
  {
    if(FileUtils::has_extension(file_path, ".vtk"))
    {
      if(is_a_poly_data)
      {
        vtkSmartPointer< vtkPolyDataWriter > writer =
            vtkSmartPointer< vtkPolyDataWriter >::New();
        writer->SetFileName(file_path.c_str());

        if(face_indices.size() > 0)
          writer->SetHeader(
              "vtk file: 2D data (surface mesh) generated by MEPP software "
              "available at https://github.com/MEPP-team/MEPP2");
        else
        {
          if(line_indices.size() == 0)
            writer->SetHeader(
                "vtk file: 1D data (point cloud) generated by MEPP software "
                "available at https://github.com/MEPP-team/MEPP2");
          else
          {
            if(line_indices[0].size() == 2)
              writer->SetHeader(
                  "vtk file: 1D data (polyline(s)) generated by MEPP software "
                  "available at https://github.com/MEPP-team/MEPP2");
            else
              writer->SetHeader(
                  "vtk file: 3D data (volume cells) generated by MEPP software "
                  "available at https://github.com/MEPP-team/MEPP2");
          }
        }

#if VTK_MAJOR_VERSION <= 5
        writer->SetInput(poly_data);
        // writer->SetInputConnection(poly_data->GetProducerPort()) ;
#else
        // VTK 6 code
        writer->SetInputData(poly_data);
#endif
        writer->Update();
        writer->Write();
      }
      else
      {
        vtkSmartPointer< vtkUnstructuredGridWriter > writer =
            vtkSmartPointer< vtkUnstructuredGridWriter >::New();
        writer->SetFileName(file_path.c_str());
        writer->SetHeader("vtk file: 3D data (volume mesh) generated by MEPP "
                          "software available at https://github.com/MEPP-team/MEPP2");

#if VTK_MAJOR_VERSION <= 5
        writer->SetInput(unstructured_grid);
        // writer->SetInputConnection(unstructured_grid->GetProducerPort()) ;
#else
        // VTK 6 code
        writer->SetInputData(unstructured_grid);
#endif
        writer->Update();
        writer->Write();
      }
    }
    else if(FileUtils::has_extension(file_path, ".vtu"))
    {
      if(!is_a_poly_data)
      {
        vtkSmartPointer< vtkXMLUnstructuredGridWriter > writer =
            vtkSmartPointer< vtkXMLUnstructuredGridWriter >::New();
        writer->SetFileName(file_path.c_str());

#if VTK_MAJOR_VERSION <= 5
        writer->SetInput(unstructured_grid);
        // writer->SetInputConnection(unstructured_grid->GetProducerPort()) ;
#else
        // VTK 6 code
        writer->SetInputData(unstructured_grid);
#endif
        writer->Update();
        writer->Write();
      }
      else
      {
        throw std::runtime_error(
            "Writer::write_vtk_or_vtp_or_vtu_file -> file extension .vtu is "
            "reserved for UnstructuredGrid.");
      }
    }
    else
    {
      throw std::runtime_error(
          "Writer::write_vtk_or_vtp_or_vtu_file -> file extension is "
          "inappropriate (neither .vtk, .vtp nor .vtu)");
    }
  }
}
} // namespace IO
} // namespace FEVV

