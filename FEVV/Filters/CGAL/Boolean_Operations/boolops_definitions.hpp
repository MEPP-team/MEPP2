// Copyright (c) 2012-2019 University of Lyon and CNRS (France).
// All rights reserved.
//
// This file is part of MEPP2; you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published 
// by the Free Software Foundation; either version 3 of the License, 
// or (at your option) any later version.
//
// This file is provided AS IS with NO WARRANTY OF ANY KIND, INCLUDING THE
// WARRANTY OF DESIGN, MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE.

#pragma once

/*!
 * \file Boolean_Operations_Definitions.h
 * \brief Definitions and functions for convertion between exact and non-exact numbers
 * \author Cyril Leconte
 */

#include <CGAL/Gmpq.h>
#include <CGAL/Lazy_exact_nt.h>

#include <chrono> // for time measurement


#include "boolops_enriched_polyhedron.hpp"


typedef typename EnrichedPolyhedron::Halfedge_handle Halfedge_handle;
typedef typename EnrichedPolyhedron::Point_3         Point3d;


/*!
 * \def BOOLEAN_OPERATIONS_DEBUG
 * \brief Enable debug info
 */
//#define BOOLEAN_OPERATIONS_DEBUG

/*!
 * \def BOOLEAN_OPERATIONS_DEBUG_VERBOSE
 * \brief Enable more debug info to compare with Mepp1
 */
//#define BOOLEAN_OPERATIONS_DEBUG_VERBOSE

/*!
 * \def BOOLEAN_OPERATIONS_TIME
 * \brief Enables computation time measuring
 */
#define BOOLEAN_OPERATIONS_TIME

/*!
 * \enum Bool_Op
 * \brief The three Boolean operations
 */
enum Bool_Op {UNION, INTER, MINUS};

/*!
 * \typedef num_type
 * \brief exact number type
 */
typedef CGAL::Lazy_exact_nt<CGAL::Gmpq>    num_type;

/*!
 * \typedef Exact_Kernel
 * \brief Kernel using exact number type
 */
typedef CGAL::Simple_cartesian<num_type>  Exact_Kernel;

/*!
 * \typedef Vector_exact
 * \brief 3d vector using exact number type
 */
typedef CGAL::Vector_3<Exact_Kernel>    Vector_exact;

/*!
 * \typedef Point3d_exact
 * \brief 3d point using exact number type
 */
typedef CGAL::Point_3<Exact_Kernel>    Point3d_exact;

/**
 * \fn inline Point3d_exact point_to_exact(Point3d &p)
 * \brief Convertion from a Point3d (double) to a Point3d_exact (exact)
 * \param p : The Point3d
 * \return The conversion in Point3d_exact.
 */
inline Point3d_exact point_to_exact(Point3d &p)
{
  return Point3d_exact(p.x(),p.y(),p.z());
}

/**
 * \fn inline Point3d point_to_double(Point3d_exact &pe)
 * \brief Convertion from a Point3d_exact (exact) to a Point3d (double)
 *
 * \param pe : The Point3d_exact
 * \return The conversion in Point3d (double).
 */
inline Point3d point_to_double(Point3d_exact &pe)
{
  return Point3d(to_double(pe.x()),to_double(pe.y()),to_double(pe.z()));
}

/**
 * \fn inline Vector_exact Compute_Normal_direction(Halfedge_handle &he)
 * \brief Compute a vector in the same direction as the normal vector
 *
 * \param he : A Halfedge incident to the facet
 * \return The normal direction (exact).
 */
inline Vector_exact Compute_Normal_direction(Halfedge_handle he)   // MT: suppression référence
{
  return CGAL::cross_product(
      point_to_exact(he->next()->vertex()->point()) -
          point_to_exact(he->vertex()->point()),
      point_to_exact(he->next()->next()->vertex()->point()) -
          point_to_exact(he->vertex()->point()));
}

/**
 * \fn inline double tr(double &n)
 * \brief Truncate a number to 1/1000
 *        (only if BOOLEAN_OPERATIONS_DEBUG is enable)
 * \param n : The input number in double
 * \return The truncation in double.
 */
inline double tr(double &n)
{
  return floor(n*1000)/1000;
}

/**
 *
 * \brief  Measure time since starting time and reset starting time.
 *
 * \param  Starting time
 * \return Elapsed time in seconds since starting time.
 */
inline double
get_time_and_reset(
    std::chrono::time_point< std::chrono::steady_clock > &time_start)
{
    auto time_now = std::chrono::steady_clock::now();
    std::chrono::duration< double > duration = time_now - time_start;
    time_start = time_now;

    return duration.count();
}