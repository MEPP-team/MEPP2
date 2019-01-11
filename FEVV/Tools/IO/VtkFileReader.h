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

#include <vtkPolyDataReader.h>
#include <vtkXMLPolyDataReader.h>
#include <vtkUnstructuredGridReader.h>
#include <vtkXMLUnstructuredGridReader.h>
#include <vtkDataSetSurfaceFilter.h>
#if defined(__APPLE__) && defined(__clang__)
#pragma clang diagnostic pop
#endif

#ifdef _MSC_VER
#define IO_TOOLS_EXPORT //"DataStructures/IO_Tools/IO_Tools_export.h"
#else
#define IO_TOOLS_EXPORT
#endif

namespace FEVV {
namespace IO {

using namespace StrUtils;
using namespace FileUtils;

inline std::string
get_line_containing_dataset(std::string file_path)
{
  ifstream myfile(file_path);
  if(myfile.is_open())
  {
    std::string line;
    while(myfile.good())
    {
      getline(myfile, line);
      if(line.find("DATASET") != std::string::npos)
      {
        myfile.close();
        return line;
      }
    }
    myfile.close();
  }
  else
  {
    throw std::runtime_error(
        "Reader::getLineContainingDATASET -> unable to open file");
  }
  return std::string();
}

template< typename CoordType,
          typename CoordNType,
          typename CoordCType,
          typename IndexType >
void
read_vtk_poly_data(
    const vtkSmartPointer< vtkPolyData > &poly_data,
    std::vector< std::vector< CoordType > > &points_coords,
    std::vector< std::vector< CoordNType > >
        &normals_coords, /// should have either 0 or points_coords.size() at the end
    std::vector< std::vector< CoordCType > > &vertex_color_coords,
    std::vector< std::vector< IndexType > >
        &line_indices, /// edges (can also be used to decribe a volume element)
    std::vector< std::vector< CoordCType > > &lines_color_coords,
    std::vector< std::vector< IndexType > >
        &face_indices, /// polygonal face made of vertex indices
    std::vector< std::vector< CoordCType > > &face_color_coords,
    std::vector< std::vector< std::vector< double > > >
        &field_attributes, // array of (2D) arrays: first point data, then cell
                           // data
    std::vector< std::string > &field_names)
{
  /////////////////////////////////////////////////////////////////////////////////////////////
  // clearing vectors...
  points_coords.clear();
  normals_coords.clear();
  vertex_color_coords.clear();
  line_indices.clear();
  lines_color_coords.clear();
  face_indices.clear();
  face_color_coords.clear();
  field_attributes.clear();
  field_names.clear();
  /////////////////////////////////////////////////////////////////////////////////////////////
  /////////////////////////////////////////////////////////////////////////////////////////////
  // some checks:
  // we do not manage triangle strips independently from cells yet
  if(poly_data->GetNumberOfStrips() > 0)
  {
    std::cerr << "Reader::read_vtkPolyData -> Warning: triangle strips not "
                 "taken into account yet.\n";
  }

  // each point there is at maximum one normal (some points may have no normal):
  assert((poly_data->GetPointData()->GetNormals() == NULL) ||
         (poly_data->GetNumberOfPoints() >=
          poly_data->GetPointData()->GetNormals()->GetNumberOfTuples()));
  /////////////////////////////////////////////////////////////////////////////////////////////
  // Get points and point normals
  if(poly_data->GetNumberOfPoints() > 0)
  {
    points_coords.reserve(static_cast< size_t >(
        poly_data->GetNumberOfPoints())); // use exactly the nb of points

    vtkSmartPointer< vtkPoints > ptr_points = poly_data->GetPoints();
    vtkSmartPointer< vtkDataArray > ptr_normals =
        poly_data->GetPointData()->GetNormals();
    bool use_norm =
        (ptr_normals != NULL &&
         ptr_points->GetNumberOfPoints() ==
             ptr_normals->GetNumberOfTuples()); // use normals only if one
                                                // normal per point/vertex

    if(use_norm)
      normals_coords.reserve(static_cast< size_t >(
          poly_data->GetPointData()
              ->GetNormals()
              ->GetNumberOfTuples())); // use exactly the nb of normals

    for(vtkIdType id = 0; id < ptr_points->GetNumberOfPoints(); id++)
    {
      std::vector< double > vd(ptr_points->GetPoint(id),
                               ptr_points->GetPoint(id) +
                                   3); // ** 3 coordinates of type double
                                       // hard-coded in vtkPoints class ** ;
      std::vector< CoordType > final_vc;
      final_vc.reserve(vd.size());
      std::vector< double >::iterator it(vd.begin()), ite(vd.end());
      for(; it != ite; ++it)
        final_vc.push_back(static_cast< CoordType >(*it));

      points_coords.push_back(final_vc);

      if(use_norm)
      {
        std::vector< double > vdn(ptr_normals->GetTuple(id),
                                  ptr_normals->GetTuple(id) +
                                      ptr_normals->GetNumberOfComponents());
        std::vector< CoordNType > final_vcn;
        final_vcn.reserve(vdn.size());
        it = vdn.begin();
        ite = vdn.end();
        for(; it != ite; ++it)
          final_vcn.push_back(static_cast< CoordNType >(*it));

        normals_coords.push_back(final_vcn);
      }
    }
  }
  /////////////////////////////////////////////////////////////////////////////////////////////
  // Get facet indices
  if(poly_data->GetNumberOfPolys() > 0)
  {
    face_indices.reserve(static_cast< size_t >(
        poly_data->GetNumberOfPolys())); // use exactly the nb of faces

    vtkSmartPointer< vtkCellArray > ptr_cell_polygons = poly_data->GetPolys();
    vtkIdType npts;
    vtkIdType *pts_poly = new vtkIdType[ptr_cell_polygons->GetMaxCellSize()];

    ptr_cell_polygons->InitTraversal();
    while(ptr_cell_polygons->GetNextCell(
        npts, pts_poly)) // 0 (i.e. false) is encountered at the list' end
    {
      std::vector< IndexType > facet_indices(npts);
      for(vtkIdType i = 0; i < npts; i++)
        facet_indices[i] = static_cast< IndexType >(pts_poly[i]);

      face_indices.push_back(facet_indices);
    }

    delete[] pts_poly;
  }
  /////////////////////////////////////////////////////////////////////////////////////////////
  // Get line indices
  if(poly_data->GetNumberOfLines() > 0)
  {
    line_indices.reserve(static_cast< size_t >(
        poly_data->GetNumberOfLines())); // use exactly the nb of line_indices

    vtkSmartPointer< vtkCellArray > ptr_cell_lines = poly_data->GetLines();
    vtkIdType npts;
    vtkIdType *pts_line = new vtkIdType[ptr_cell_lines->GetMaxCellSize()];

    ptr_cell_lines->InitTraversal();
    while(ptr_cell_lines->GetNextCell(
        npts, pts_line)) // 0 (i.e. false) is encountered at the list' end
    {
      std::vector< IndexType > line_indices_tmp(npts);
      for(vtkIdType i = 0; i < npts; i++)
        line_indices_tmp[i] = static_cast< IndexType >(pts_line[i]);

      line_indices.push_back(line_indices_tmp);
    }

    delete[] pts_line;
  }
  /////////////////////////////////////////////////////////////////////////////////////////////
  // Get 2D arrays of FIELD data when present
  vtkFieldData *fd = poly_data->GetFieldData();
  if(fd != NULL)
  {
    field_attributes.resize(static_cast< size_t >(
        fd->GetNumberOfArrays())); // use exactly the nb of 2D arrays for FIELD
                                   // data
    field_names.resize(static_cast< size_t >(fd->GetNumberOfArrays()));

    for(int id_array = 0; id_array < fd->GetNumberOfArrays(); id_array++)
    {
      field_names[id_array] = fd->GetArrayName(id_array);
      // for each FIELD data (arrays)
      vtkSmartPointer< vtkDataArray > ptr_data = fd->GetArray(
          id_array); // poly_data->GetFieldData()->GetArrayName(idArray) ) ;

      field_attributes[id_array].resize(
          static_cast< size_t >(ptr_data->GetNumberOfTuples()));
      // std::cout << "ptrData->GetNumberOfTuples() = " <<
      // ptrData->GetNumberOfTuples() << std::endl ;
      for(vtkIdType id_tuple_in_array = 0;
          id_tuple_in_array < ptr_data->GetNumberOfTuples();
          id_tuple_in_array++)
      {
        // field_attributes[idArray] is a 2D array
        field_attributes[id_array][id_tuple_in_array].resize(
            static_cast< size_t >(ptr_data->GetNumberOfComponents()));
        std::copy(ptr_data->GetTuple(id_tuple_in_array),
                  ptr_data->GetTuple(id_tuple_in_array) +
                      ptr_data->GetNumberOfComponents(),

                  field_attributes[id_array][id_tuple_in_array].begin());
      }
    }
  }
  /////////////////////////////////////////////////////////////////////////////////////////////
  // Get point data when present
  vtkPointData *pd = poly_data->GetPointData();
  if(pd != NULL)
  {
    size_t old_nb_e = field_attributes.size();
    field_attributes.resize(static_cast< size_t >(
        old_nb_e + pd->GetNumberOfArrays())); // use exactly the nb of 2D arrays
                                              // for FIELD data
    field_names.resize(
        static_cast< size_t >(old_nb_e + pd->GetNumberOfArrays()));
    for(int id_array = 0; id_array < pd->GetNumberOfArrays(); id_array++)
    {
      field_names[id_array + old_nb_e] =
          std::string("POINT_DATA_") + pd->GetArrayName(id_array);
      // std::cout << "POINT_DATA: array name is " << pd->GetArrayName(idArray)
      // << std::endl;
      // for each FIELD data (arrays)
      vtkSmartPointer< vtkDataArray > ptr_data = pd->GetArray(id_array);

      field_attributes[id_array + old_nb_e].resize(
          static_cast< size_t >(ptr_data->GetNumberOfTuples()));
      // std::cout << "ptrData->GetNumberOfTuples() = " <<
      // ptrData->GetNumberOfTuples() << std::endl ;
      for(vtkIdType id_tuple_in_array = 0;
          id_tuple_in_array < ptr_data->GetNumberOfTuples();
          id_tuple_in_array++)
      {
        // field_attributes[idArray] is a 2D array
        field_attributes[id_array + old_nb_e][id_tuple_in_array].resize(
            static_cast< size_t >(ptr_data->GetNumberOfComponents()));
        std::copy(
            ptr_data->GetTuple(id_tuple_in_array),
            ptr_data->GetTuple(id_tuple_in_array) +
                ptr_data->GetNumberOfComponents(),

            field_attributes[id_array + old_nb_e][id_tuple_in_array].begin());
      }
    }
  }
  ////////////////////////////////////////////////////////////////////////////
  // Get cell data when present
  vtkCellData *cd = poly_data->GetCellData();
  if(cd != NULL)
  {
    ////////////////////////////////////////////////////////////////////////////
    // arrays (SCALARS are taken into account there)
    size_t old_nb_e = field_attributes.size();
    field_attributes.resize(static_cast< size_t >(
        old_nb_e + cd->GetNumberOfArrays())); // use exactly the nb of 2D arrays
                                              // for FIELD data
    field_names.resize(
        static_cast< size_t >(old_nb_e + cd->GetNumberOfArrays()));
    for(int id_array = 0; id_array < cd->GetNumberOfArrays(); id_array++)
    {
      field_names[old_nb_e + id_array] =
          std::string("CELL_DATA_") + cd->GetArrayName(id_array);
      // std::cout << "CELL_DATA: array name is " << cd->GetArrayName(idArray)
      // << std::endl;
      // for each FIELD data (arrays)
      vtkSmartPointer< vtkDataArray > ptr_data = cd->GetArray(id_array);

      field_attributes[old_nb_e + id_array].resize(
          static_cast< size_t >(ptr_data->GetNumberOfTuples()));
      // std::cout << "ptrData->GetNumberOfTuples() = " <<
      // ptrData->GetNumberOfTuples() << std::endl ;
      for(vtkIdType id_tuple_in_array = 0;
          id_tuple_in_array < ptr_data->GetNumberOfTuples();
          id_tuple_in_array++)
      {
        // field_attributes[oldNbE+idArray] is a 2D array
        field_attributes[old_nb_e + id_array][id_tuple_in_array].resize(
            static_cast< size_t >(ptr_data->GetNumberOfComponents()));
        std::copy(
            ptr_data->GetTuple(id_tuple_in_array),
            ptr_data->GetTuple(id_tuple_in_array) +
                ptr_data->GetNumberOfComponents(),

            field_attributes[old_nb_e + id_array][id_tuple_in_array].begin());
      }
    }
  }
  if((line_indices.size() > 0) && // line_indices are present
     (line_indices[0].size() >
      2) &&                     // but at least 3 vertex for an edge (strange)
     (face_indices.size() == 0) // and no face_indices provided...
  )
  { // line_indices in 2D can sometimes be misused as polygons by some users of
    // MEPP (but for vtk line_indices refer to edges)
    std::swap(line_indices, face_indices);
  }
  ////////////////////////////////////////////////////////////////////////////
  if(field_names.size() == field_attributes.size())
  {
    size_t cpt = 0, nb_e_glob = field_names.size();
    for(; cpt < nb_e_glob; ++cpt)
    {
      if(field_names[cpt].find("colo") != std::string::npos)
      {
        if((vertex_color_coords.size() == 0) &&
           ((field_names[cpt].find("vertex") != std::string::npos) ||
            (field_names[cpt].find("point") != std::string::npos)))
        {
          size_t nb_e = points_coords.size();
          vertex_color_coords.resize(nb_e);
          for(size_t i = 0; i < nb_e; ++i)
          {
            size_t nb_c = field_attributes[cpt][i].size();
            vertex_color_coords[i].resize(nb_c);
            for(size_t j = 0; j < nb_c; ++j)
              vertex_color_coords[i][j] =
                  static_cast< CoordCType >(field_attributes[cpt][i][j]);
            // std::copy(field_attributes[cpt][i].begin(),
            // field_attributes[cpt][i].end(), vertex_color_coords[i].begin()) ;
            // std::copy(field_attributes[cpt][i].begin(),
            // field_attributes[cpt][i].end(), std::back_inserter(
            // vertex_color_coords[i] )) ;
          }
        }
        else if((face_color_coords.size() == 0) &&
                ((field_names[cpt].find("face") != std::string::npos) ||
                 (field_names[cpt].find("polygon") != std::string::npos) ||
                 (line_indices.size() == 0)))
        {
          size_t nb_e = face_indices.size();
          face_color_coords.resize(nb_e);
          for(size_t i = 0; i < nb_e; ++i)
          {
            size_t nb_c = field_attributes[cpt][i].size();
            face_color_coords[i].resize(nb_c);
            for(size_t j = 0; j < nb_c; ++j)
              face_color_coords[i][j] =
                  static_cast< CoordCType >(field_attributes[cpt][i][j]);
            // std::copy(field_attributes[cpt][i].begin(),
            // field_attributes[cpt][i].end(), face_color_coords[i].begin()) ;
            // std::copy(field_attributes[cpt][i].begin(),
            // field_attributes[cpt][i].end(), std::back_inserter(
            // face_color_coords[i] )) ;
          }
        }
        else if((lines_color_coords.size() == 0) && (line_indices.size() > 0))
        {
          size_t nb_e = line_indices.size();
          lines_color_coords.resize(nb_e);
          for(size_t i = 0; i < nb_e; ++i)
          {
            size_t nb_c = field_attributes[cpt][i].size();
            lines_color_coords[i].resize(nb_c);
            for(size_t j = 0; j < nb_c; ++j)
              lines_color_coords[i][j] =
                  static_cast< CoordCType >(field_attributes[cpt][i][j]);
            // std::copy(field_attributes[cpt][i].begin(),
            // field_attributes[cpt][i].end(), lines_color_coords[i].begin()) ;
            // std::copy(field_attributes[cpt][i].begin(),
            // field_attributes[cpt][i].end(), std::back_inserter(
            // lines_color_coords[i] )) ;
          }
        }
        field_attributes[cpt].clear(); // to not take into account twice colors!

        size_t last_index = nb_e_glob - 1;
        if(cpt < last_index)
        {
          std::swap(field_attributes[cpt], field_attributes[last_index]);
          std::swap(field_names[cpt], field_names[last_index]);
        }
        field_attributes.resize(last_index);
        field_names.resize(last_index);
        --nb_e_glob; // update end iterator
        --cpt;
      }
      else if(field_names[cpt].find("norm") != std::string::npos)
      { // we do not take into account normals twice!
        field_attributes[cpt].clear(); // to not take into account twice colors!

        size_t last_index = nb_e_glob - 1;
        if(cpt < last_index)
        {
          std::swap(field_attributes[cpt], field_attributes[last_index]);
          std::swap(field_names[cpt], field_names[last_index]);
        }
        field_attributes.resize(last_index);
        field_names.resize(last_index);
        --nb_e_glob; // update end iterator
        --cpt;
      }
    }
  }
}

template< typename CoordType,
          typename CoordNType,
          typename CoordCType,
          typename IndexType >
void
read_vtk_unstructured_grid(
    const vtkSmartPointer< vtkUnstructuredGrid > &unstructured_grid,
    std::vector< std::vector< CoordType > > &points_coords,
    std::vector< std::vector< CoordNType > >
        &normals_coords, /// should have either 0 or points_coords.size() at the
                         /// end
    std::vector< std::vector< CoordCType > > &vertex_color_coords,
    std::vector< std::vector< IndexType > >
        &line_indices, /// edges (can also be used to decribe a volume element)
    std::vector< std::vector< CoordCType > > &lines_color_coords,
    std::vector< std::vector< IndexType > >
        &face_indices, /// polygonal facets made of vertex indices
    std::vector< std::vector< CoordCType > > &face_color_coords,
    std::vector< std::vector< std::vector< double > > >
        &field_attributes, // array of (2D) arrays: first point data, then cell
                           // data
    std::vector< std::string > &field_names)
{
  /////////////////////////////////////////////////////////////////////////////////////////////
  // clearing vectors...
  points_coords.clear();
  normals_coords.clear();
  vertex_color_coords.clear();
  line_indices.clear();
  lines_color_coords.clear();
  face_indices.clear();
  face_color_coords.clear();
  field_attributes.clear();
  /////////////////////////////////////////////////////////////////////////////////////////////

  // each point there is at maximum one normal (some points may have no normal):
  assert(
      (unstructured_grid->GetPointData()->GetNormals() == NULL) ||
      (unstructured_grid->GetNumberOfPoints() >=
       unstructured_grid->GetPointData()->GetNormals()->GetNumberOfTuples()));
  /////////////////////////////////////////////////////////////////////////////////////////////
  // Get points and point normals
  if(unstructured_grid->GetNumberOfPoints() > 0)
  {
    points_coords.reserve(static_cast< size_t >(
        unstructured_grid
            ->GetNumberOfPoints())); // use exactly the nb of points

    vtkSmartPointer< vtkPoints > ptr_points = unstructured_grid->GetPoints();
    vtkSmartPointer< vtkDataArray > ptr_normals =
        unstructured_grid->GetPointData()->GetNormals();
    bool use_norm =
        (ptr_normals != NULL &&
         ptr_points->GetNumberOfPoints() ==
             ptr_normals->GetNumberOfTuples()); // use normals only if one
                                                // normal per point/vertex

    if(use_norm)
      normals_coords.reserve(static_cast< size_t >(
          unstructured_grid->GetPointData()
              ->GetNormals()
              ->GetNumberOfTuples())); // use exactly the nb of normals

    for(vtkIdType id = 0; id < ptr_points->GetNumberOfPoints(); id++)
    {
      std::vector< double > vd(ptr_points->GetPoint(id),
                               ptr_points->GetPoint(id) +
                                   3); // ** 3 coordinates of type double
                                       // hard-coded in vtkPoints class ** ;
      std::vector< CoordType > final_vc;
      final_vc.reserve(vd.size());
      std::vector< double >::iterator it(vd.begin()), ite(vd.end());
      for(; it != ite; ++it)
        final_vc.push_back(static_cast< CoordType >(*it));

      points_coords.push_back(final_vc);

      if(use_norm)
      {
        std::vector< double > vdn(ptr_normals->GetTuple(id),
                                  ptr_normals->GetTuple(id) +
                                      ptr_normals->GetNumberOfComponents());
        std::vector< CoordNType > final_vcn;
        final_vcn.reserve(vdn.size());
        it = vdn.begin();
        ite = vdn.end();
        for(; it != ite; ++it)
          final_vcn.push_back(static_cast< CoordNType >(*it));

        normals_coords.push_back(final_vcn);
      }
    }
  }
  /////////////////////////////////////////////////////////////////////////////////////////////
  // Get facet indices
  if(unstructured_grid->GetNumberOfCells() > 0)
  {
    face_indices.reserve(static_cast< size_t >(
        unstructured_grid->GetNumberOfCells())); // use exactly the nb of faces
    line_indices.reserve(
        static_cast< size_t >(unstructured_grid->GetNumberOfCells()));

    size_t nb_face_cell = 0, nb_vol_cell = 0;

    vtkSmartPointer< vtkCellArray > ptr_cell_polygons =
        unstructured_grid->GetCells();
    vtkIdType npts;
    vtkIdType *pts_poly = new vtkIdType[ptr_cell_polygons->GetMaxCellSize()];

    ptr_cell_polygons->InitTraversal();
    vtkIdType cell_id = 0;
    while(ptr_cell_polygons->GetNextCell(
        npts, pts_poly)) // 0 (i.e. false) is encountered at the list' end
    {
      std::vector< IndexType > facet_indices(npts);
      for(vtkIdType i = 0; i < npts; i++)
      {
        facet_indices[i] = static_cast< IndexType >(pts_poly[i]);
        // std::cout << static_cast<Index>( ptsPoly[i] ) << std::endl ;
      }

      switch(unstructured_grid->GetCellType(cell_id))
      {
      case VTK_VERTEX: // vertex (nothing to do)
        throw std::runtime_error("Reader::read_vtkUnstructuredGrid -> vertex "
                                 "cells are not managed.");
        // break;
      case VTK_LINE: // line/edge (nothing to do, because edge alone are not of
                     // interest for numerical simulation)
        throw std::runtime_error(
            "Reader::read_vtkUnstructuredGrid -> edge cells are not managed.");
        // break;
      case VTK_TRIANGLE: // triangle
        face_indices.push_back(facet_indices);
        nb_face_cell++;
        break;
      case VTK_POLYGON: // polygon
        face_indices.push_back(facet_indices);
        nb_face_cell++;
        break;
      case VTK_QUAD: // quad
        face_indices.push_back(facet_indices);
        nb_face_cell++;
        break;
      case VTK_TETRA: // tetra
      {
        line_indices.push_back(facet_indices);
        nb_vol_cell++;
        std::vector< IndexType > face1, face2, face3, face4;
        face1.push_back(facet_indices[0]);
        face1.push_back(facet_indices[1]);
        face1.push_back(facet_indices[3]);
        face_indices.push_back(face1);
        nb_face_cell++; // ok

        face2.push_back(facet_indices[1]);
        face2.push_back(facet_indices[2]);
        face2.push_back(facet_indices[3]);
        face_indices.push_back(face2);
        nb_face_cell++; // ok

        face3.push_back(facet_indices[3]);
        face3.push_back(facet_indices[2]);
        face3.push_back(facet_indices[0]);
        face_indices.push_back(face3);
        nb_face_cell++; // ok

        face4.push_back(facet_indices[2]);
        face4.push_back(facet_indices[1]);
        face4.push_back(facet_indices[0]);
        face_indices.push_back(face4);
        nb_face_cell++; // ok
      }
      break;
      case VTK_HEXAHEDRON: // hexa
      {
        line_indices.push_back(facet_indices);
        nb_vol_cell++;
        std::vector< IndexType > face1, face2, face3, face4, face5, face6;
        face1.push_back(facet_indices[0]);
        face1.push_back(facet_indices[1]);
        face1.push_back(facet_indices[5]);
        face1.push_back(facet_indices[4]);
        face_indices.push_back(face1);
        nb_face_cell++; // ok

        face2.push_back(facet_indices[1]);
        face2.push_back(facet_indices[2]);
        face2.push_back(facet_indices[6]);
        face2.push_back(facet_indices[5]);
        face_indices.push_back(face2);
        nb_face_cell++; // ok

        face3.push_back(facet_indices[2]);
        face3.push_back(facet_indices[3]);
        face3.push_back(facet_indices[7]);
        face3.push_back(facet_indices[6]);
        face_indices.push_back(face3);
        nb_face_cell++; // ok

        face4.push_back(facet_indices[3]);
        face4.push_back(facet_indices[0]);
        face4.push_back(facet_indices[4]);
        face4.push_back(facet_indices[7]);
        face_indices.push_back(face4);
        nb_face_cell++; // ok

        face5.push_back(facet_indices[4]);
        face5.push_back(facet_indices[5]);
        face5.push_back(facet_indices[6]);
        face5.push_back(facet_indices[7]);
        face_indices.push_back(face5);
        nb_face_cell++; // ok

        face6.push_back(facet_indices[3]);
        face6.push_back(facet_indices[2]);
        face6.push_back(facet_indices[1]);
        face6.push_back(facet_indices[0]);
        face_indices.push_back(face6);
        nb_face_cell++; // ok
      }
      break;
      /*
      case VTK_WEDGE: //wedge/triangular prism
      line_indices.push_back( facetIndices ) ; nbVolCell++; break ;
      case VTK_PYRAMID: //pyramid
      line_indices.push_back( facetIndices ) ; nbVolCell++; break ;*/
      default:
        throw std::runtime_error(
            "Reader::read_vtkUnstructuredGrid -> cell type is not known.");
      }

      cell_id++;
    }
    face_indices.resize(nb_face_cell);
    line_indices.resize(nb_vol_cell);

    delete[] pts_poly;
  }
  /////////////////////////////////////////////////////////////////////////////////////////////
  // Get 2D arrays of FIELD data when present
  vtkFieldData *fd = unstructured_grid->GetFieldData();
  if(fd != NULL)
  {
    field_attributes.resize(static_cast< size_t >(
        fd->GetNumberOfArrays())); // use exactly the nb of 2D arrays for FIELD
                                   // data
    field_names.resize(static_cast< size_t >(fd->GetNumberOfArrays()));

    for(int id_array = 0; id_array < fd->GetNumberOfArrays(); id_array++)
    {
      field_names[id_array] = fd->GetArrayName(id_array);
      // for each FIELD data (arrays)
      vtkSmartPointer< vtkDataArray > ptr_data = fd->GetArray(
          id_array); // poly_data->GetFieldData()->GetArrayName(idArray) ) ;

      field_attributes[id_array].resize(
          static_cast< size_t >(ptr_data->GetNumberOfTuples()));
      // std::cout << "ptrData->GetNumberOfTuples() = " <<
      // ptrData->GetNumberOfTuples() << std::endl ;
      for(vtkIdType id_tuple_in_array = 0;
          id_tuple_in_array < ptr_data->GetNumberOfTuples();
          id_tuple_in_array++)
      {
        // field_attributes[idArray] is a 2D array
        field_attributes[id_array][id_tuple_in_array].resize(
            static_cast< size_t >(ptr_data->GetNumberOfComponents()));
        std::copy(ptr_data->GetTuple(id_tuple_in_array),
                  ptr_data->GetTuple(id_tuple_in_array) +
                      ptr_data->GetNumberOfComponents(),

                  field_attributes[id_array][id_tuple_in_array].begin());
      }
    }
  }
  /////////////////////////////////////////////////////////////////////////////////////////////
  // Get point data when present
  vtkPointData *pd = unstructured_grid->GetPointData();
  if(pd != NULL)
  {
    size_t old_nb_e = field_attributes.size();
    field_attributes.resize(static_cast< size_t >(
        old_nb_e + pd->GetNumberOfArrays())); // use exactly the nb of 2D arrays
                                              // for FIELD data
    field_names.resize(
        static_cast< size_t >(old_nb_e + pd->GetNumberOfArrays()));
    for(int id_array = 0; id_array < pd->GetNumberOfArrays(); id_array++)
    {
      field_names[id_array + old_nb_e] =
          std::string("POINT_DATA_") + pd->GetArrayName(id_array);
      // std::cout << "POINT_DATA: array name is " << pd->GetArrayName(idArray)
      // << std::endl;
      // for each FIELD data (arrays)
      vtkSmartPointer< vtkDataArray > ptr_data = pd->GetArray(id_array);

      field_attributes[id_array + old_nb_e].resize(
          static_cast< size_t >(ptr_data->GetNumberOfTuples()));
      // std::cout << "ptrData->GetNumberOfTuples() = " <<
      // ptrData->GetNumberOfTuples() << std::endl ;
      for(vtkIdType id_tuple_in_array = 0;
          id_tuple_in_array < ptr_data->GetNumberOfTuples();
          id_tuple_in_array++)
      {
        // field_attributes[idArray] is a 2D array
        field_attributes[id_array + old_nb_e][id_tuple_in_array].resize(
            static_cast< size_t >(ptr_data->GetNumberOfComponents()));
        std::copy(
            ptr_data->GetTuple(id_tuple_in_array),
            ptr_data->GetTuple(id_tuple_in_array) +
                ptr_data->GetNumberOfComponents(),

            field_attributes[id_array + old_nb_e][id_tuple_in_array].begin());
      }
    }
  }
  ////////////////////////////////////////////////////////////////////////////
  // Get cell data when present
  vtkCellData *cd = unstructured_grid->GetCellData();
  if(cd != NULL)
  {
    ////////////////////////////////////////////////////////////////////////////
    // arrays (SCALARS are taken into account there)
    size_t old_nb_e = field_attributes.size();
    field_attributes.resize(static_cast< size_t >(
        old_nb_e + cd->GetNumberOfArrays())); // use exactly the nb of 2D arrays
                                              // for FIELD data
    field_names.resize(
        static_cast< size_t >(old_nb_e + cd->GetNumberOfArrays()));
    for(int id_array = 0; id_array < cd->GetNumberOfArrays(); id_array++)
    {
      field_names[old_nb_e + id_array] =
          std::string("CELL_DATA_") + cd->GetArrayName(id_array);
      // std::cout << "CELL_DATA: array name is " << cd->GetArrayName(idArray)
      // << std::endl;
      // for each FIELD data (arrays)
      vtkSmartPointer< vtkDataArray > ptr_data = cd->GetArray(id_array);

      field_attributes[old_nb_e + id_array].resize(
          static_cast< size_t >(ptr_data->GetNumberOfTuples()));
      // std::cout << "ptrData->GetNumberOfTuples() = " <<
      // ptrData->GetNumberOfTuples() << std::endl ;
      for(vtkIdType id_tuple_in_array = 0;
          id_tuple_in_array < ptr_data->GetNumberOfTuples();
          id_tuple_in_array++)
      {
        // field_attributes[oldNbE+idArray] is a 2D array
        field_attributes[old_nb_e + id_array][id_tuple_in_array].resize(
            static_cast< size_t >(ptr_data->GetNumberOfComponents()));
        std::copy(
            ptr_data->GetTuple(id_tuple_in_array),
            ptr_data->GetTuple(id_tuple_in_array) +
                ptr_data->GetNumberOfComponents(),

            field_attributes[old_nb_e + id_array][id_tuple_in_array].begin());
      }
    }
  }
  ////////////////////////////////////////////////////////////////////////////
  if(field_names.size() == field_attributes.size())
  {
    size_t cpt = 0, nb_e_glob = field_names.size();
    for(; cpt < nb_e_glob; ++cpt)
    {
      if(field_names[cpt].find("colo") != std::string::npos)
      {
        if((vertex_color_coords.size() == 0) &&
           ((field_names[cpt].find("vertex") != std::string::npos) ||
            (field_names[cpt].find("point") != std::string::npos)))
        {
          size_t nb_e = points_coords.size();
          vertex_color_coords.resize(nb_e);
          for(size_t i = 0; i < nb_e; ++i)
          {
            size_t nb_c = field_attributes[cpt][i].size();
            vertex_color_coords[i].resize(nb_c);
            for(size_t j = 0; j < nb_c; ++j)
              vertex_color_coords[i][j] =
                  static_cast< CoordCType >(field_attributes[cpt][i][j]);
            // std::copy(field_attributes[cpt][i].begin(),
            // field_attributes[cpt][i].end(), vertex_color_coords[i].begin()) ;
            // std::copy(field_attributes[cpt][i].begin(),
            // field_attributes[cpt][i].end(), std::back_inserter(
            // vertex_color_coords[i] )) ;
          }
        }
        else if((face_color_coords.size() == 0) &&
                ((field_names[cpt].find("face") != std::string::npos) ||
                 (field_names[cpt].find("polygon") != std::string::npos) ||
                 (line_indices.size() == 0)))
        {
          size_t nb_e = face_indices.size();
          face_color_coords.resize(nb_e);
          for(size_t i = 0; i < nb_e; ++i)
          {
            size_t nb_c = field_attributes[cpt][i].size();
            face_color_coords[i].resize(nb_c);
            for(size_t j = 0; j < nb_c; ++j)
              face_color_coords[i][j] =
                  static_cast< CoordCType >(field_attributes[cpt][i][j]);
            // std::copy(field_attributes[cpt][i].begin(),
            // field_attributes[cpt][i].end(), face_color_coords[i].begin()) ;
            // std::copy(field_attributes[cpt][i].begin(),
            // field_attributes[cpt][i].end(), std::back_inserter(
            // face_color_coords[i] )) ;
          }
        }
        else if((lines_color_coords.size() == 0) && (line_indices.size() > 0))
        {
          size_t nb_e = line_indices.size();
          lines_color_coords.resize(nb_e);
          for(size_t i = 0; i < nb_e; ++i)
          {
            size_t nb_c = field_attributes[cpt][i].size();
            lines_color_coords[i].resize(nb_c);
            for(size_t j = 0; j < nb_c; ++j)
              lines_color_coords[i][j] =
                  static_cast< CoordCType >(field_attributes[cpt][i][j]);
            // std::copy(field_attributes[cpt][i].begin(),
            // field_attributes[cpt][i].end(), lines_color_coords[i].begin()) ;
            // std::copy(field_attributes[cpt][i].begin(),
            // field_attributes[cpt][i].end(), std::back_inserter(
            // lines_color_coords[i] )) ;
          }
        }
        field_attributes[cpt].clear(); // to not take into account twice colors!

        size_t last_index = nb_e_glob - 1;
        if(cpt < last_index)
        {
          std::swap(field_attributes[cpt], field_attributes[last_index]);
          std::swap(field_names[cpt], field_names[last_index]);
        }
        field_attributes.resize(last_index);
        field_names.resize(last_index);
        --nb_e_glob; // update end iterator
        --cpt;
      }
      else if(field_names[cpt].find("norm") != std::string::npos)
      { // we do not take into account normals twice!
        field_attributes[cpt].clear(); // to not take into account twice colors!

        size_t last_index = nb_e_glob - 1;
        if(cpt < last_index)
        {
          std::swap(field_attributes[cpt], field_attributes[last_index]);
          std::swap(field_names[cpt], field_names[last_index]);
        }
        field_attributes.resize(last_index);
        field_names.resize(last_index);
        --nb_e_glob; // update end iterator
        --cpt;
      }
    }
  }
}

template< typename CoordType,
          typename CoordNType,
          typename CoordCType,
          typename IndexType >
void
read_vtk_or_vtp_or_vtu_file(
    std::string file_path,
    std::vector< std::vector< CoordType > > &points_coords,
    std::vector< std::vector< CoordNType > >
        &normals_coords, /// should have either 0 or points_coords.size() at the
                         /// end
    std::vector< std::vector< CoordCType > > &vertex_color_coords,
    std::vector< std::vector< IndexType > >
        &line_indices, /// edges (can also be used to decribe a volume element)
    std::vector< std::vector< CoordCType > > &lines_color_coords,
    std::vector< std::vector< IndexType > >
        &face_indices, /// polygonal facets made of vertex indices
    std::vector< std::vector< CoordCType > > &face_color_coords,
    std::vector< std::vector< std::vector< double > > >
        &field_attributes, // array of 2D arrays: first point data, then cell
                           // data
    std::vector< std::string > &field_names,
    bool if_unstructured_grid_then_take_surface = false)
{
  /////////////////////////////////////////////////////////////////////////////////////////////
  // reading and getting polydata/unstructured mesh
  bool is_a_poly_data = true;
  vtkSmartPointer< vtkPolyData > poly_data;
  vtkSmartPointer< vtkUnstructuredGrid > unstructured_grid;

  std::string dataset_line = get_line_containing_dataset(file_path);

  if(FileUtils::has_extension(file_path, ".vtp"))
  {
    vtkSmartPointer< vtkXMLPolyDataReader > reader =
        vtkSmartPointer< vtkXMLPolyDataReader >::New();
    reader->SetFileName(file_path.c_str());
    reader->Update();

    poly_data = vtkSmartPointer< vtkPolyData >::New();
    poly_data->DeepCopy(reader->GetOutput());
  }
  else if(FileUtils::has_extension(file_path, ".vtu"))
  {
    vtkSmartPointer< vtkXMLUnstructuredGridReader > reader =
        vtkSmartPointer< vtkXMLUnstructuredGridReader >::New();
    reader->SetFileName(file_path.c_str());
    reader->Update();

    if(if_unstructured_grid_then_take_surface)
    {
      vtkSmartPointer< vtkDataSetSurfaceFilter > surface_filter =
          vtkSmartPointer< vtkDataSetSurfaceFilter >::New();
#if VTK_MAJOR_VERSION <= 5
      surfaceFilter->SetInput(reader->GetOutput());
#else
      surface_filter->SetInputData(reader->GetOutput());
#endif
      surface_filter->Update();

      poly_data = vtkSmartPointer< vtkPolyData >::New();
      poly_data->DeepCopy(surface_filter->GetOutput());
    }
    else
    {
      is_a_poly_data = false;
      unstructured_grid = vtkSmartPointer< vtkUnstructuredGrid >::New();
      unstructured_grid->DeepCopy(reader->GetOutput());
    }
  }
  else
  {
    if(FileUtils::has_extension(file_path, ".vtk"))
    {
      if(dataset_line.find("POLYDATA") != std::string::npos)
      {
        vtkSmartPointer< vtkPolyDataReader > reader =
            vtkSmartPointer< vtkPolyDataReader >::New();
        reader->SetFileName(file_path.c_str());
        reader->Update();

        poly_data = vtkSmartPointer< vtkPolyData >::New();
        poly_data->DeepCopy(reader->GetOutput());
      }
      else
      {
        if(dataset_line.find("UNSTRUCTURED_GRID") != std::string::npos)
        {
          vtkSmartPointer< vtkUnstructuredGridReader > reader =
              vtkSmartPointer< vtkUnstructuredGridReader >::New();
          reader->SetFileName(file_path.c_str());
          reader->Update();

          if(if_unstructured_grid_then_take_surface)
          {
            vtkSmartPointer< vtkDataSetSurfaceFilter > surface_filter =
                vtkSmartPointer< vtkDataSetSurfaceFilter >::New();
#if VTK_MAJOR_VERSION <= 5
            surfaceFilter->SetInput(reader->GetOutput());
#else
            surface_filter->SetInputData(reader->GetOutput());
#endif
            surface_filter->Update();

            poly_data = vtkSmartPointer< vtkPolyData >::New();
            poly_data->DeepCopy(surface_filter->GetOutput());
          }
          else
          {
            is_a_poly_data = false;
            unstructured_grid = vtkSmartPointer< vtkUnstructuredGrid >::New();
            unstructured_grid->DeepCopy(reader->GetOutput());
          }
        }
        else
        {
          throw std::runtime_error(
              "Reader::read_vtk_or_vtp_or_vtu_file -> file dataset type cannot "
              "be processed (neither POLYDATA nor UNSTRUCTURED_GRID)");
        }
      }
    }
    else
    {
      throw std::runtime_error(
          "Reader::read_vtk_or_vtp_or_vtu_file -> file extension is "
          "inappropriate (neither .vtk nor .vtp)");
    }
  }
  /////////////////////////////////////////////////////////////////////////////
  if(is_a_poly_data)
    return read_vtk_poly_data< CoordType,
                               CoordNType,
                               /*coordT_type,*/ CoordCType,
                               IndexType >(poly_data,
                                           points_coords,
                                           normals_coords,
                                           vertex_color_coords,
                                           line_indices,
                                           lines_color_coords,
                                           face_indices,
                                           face_color_coords,
                                           field_attributes,
                                           field_names);
  else
    return read_vtk_unstructured_grid< CoordType,
                                       CoordNType,
                                       /*coordT_type,*/ CoordCType,
                                       IndexType >(unstructured_grid,
                                                   points_coords,
                                                   normals_coords,
                                                   vertex_color_coords,
                                                   line_indices,
                                                   lines_color_coords,
                                                   face_indices,
                                                   face_color_coords,
                                                   field_attributes,
                                                   field_names);
  /////////////////////////////////////////////////////////////////////////////
}


} // namespace IO
} // namespace FEVV

