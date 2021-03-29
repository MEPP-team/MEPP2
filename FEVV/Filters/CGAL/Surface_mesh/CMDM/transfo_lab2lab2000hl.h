// Copyright (c) 2012-2020 University of Lyon and CNRS (France).
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

#include <FEVV/Filters/CGAL/Surface_mesh/CMDM/resources/L_data.h>
#include <FEVV/Filters/CGAL/Surface_mesh/CMDM/resources/RegularGrid_0_0_1.h>
#include <FEVV/Filters/CGAL/Surface_mesh/CMDM/resources/RegularGrid_0_0_2.h>
#include <FEVV/Filters/CGAL/Surface_mesh/CMDM/resources/RegularGridInit_0_0_1.h>
#include <FEVV/Filters/CGAL/Surface_mesh/CMDM/resources/RegularGridInit_0_0_2.h>

inline double
interpolate1_computevalue(double x0, double x1, double y0, double y1, double x)
{
  return y0 + ((x - x0) / (x1 - x0)) * (y1 - y0);
}

inline double
interpolate1_process(std::vector< double > &init_grid,
                     std::vector< double > &val_grid,
                     double x) // std::vector<double> & result)
{

  auto upper = std::upper_bound(init_grid.begin(), init_grid.end(), x);
  size_t position =
      static_cast< size_t >(std::distance(init_grid.begin(), upper));

  // notes-VVI-ELO on possible bugs:
  // - std::upper_bound must be applied on a sorted container ; is init_grid
  //   really sorted?
  // - if upper is init_grid.begin() (aka the first element of init_grid is
  //   greater than x), then position equals 0 ; the lines below accessing
  //   init_grid[position - 1] and val_grid[position - 1] will cause a
  //   container overflow

  if(position < init_grid.size())
  {
    double res = interpolate1_computevalue(init_grid[position - 1],
                                           init_grid[position],
                                           val_grid[position - 1],
                                           val_grid[position],
                                           x);
    return res;
  }
  else // if upper out of bounds returning maxval
  {

    return val_grid[position - 1];
  }
}

inline double
interpolate2_computevalue(double q_00,
                          double q_10,
                          double q_01,
                          double q_11,
                          double x0,
                          double x1,
                          double y0,
                          double y1,
                          double x,
                          double y)
{
  // https://helloacm.com
  double x1x0, y1y0, x1x, y1y, yy0, xx0;
  x1x0 = x1 - x0;
  y1y0 = y1 - y0;
  x1x = x1 - x;
  y1y = y1 - y;
  yy0 = y - y0;
  xx0 = x - x0;
  return 1.0 / (x1x0 * y1y0) *
         (q_00 * x1x * y1y + q_10 * xx0 * y1y + q_01 * x1x * yy0 +
          q_11 * xx0 * yy0);


  // return  y0 + ((x - x0) / (x1 - x0)) * (y1 - y0);
}

inline std::pair< double, double >
interpolate2_process(
    std::vector< std::vector< std::pair< double, double > > > &init_grid_AB,
    double x,
    double y) // std::vector<double> & result)
{
  // Finding range
  int x_min = (int)std::floor(x);
  int x_max = x_min + 1; // int x_max = (int)std::ceil(x);
  int y_min = (int)std::floor(y);
  int y_max = y_min + 1; // int y_max = (int)std::ceil(y);

  // If there is no need to interpolate
  if(x_min == x_max && y_min == y_max)
  {
    // return the value
    return std::make_pair(init_grid_AB[x_min + 128][y_min + 128].first,
                          init_grid_AB[x_min + 128][y_min + 128].second);
  }


  double a_interpolated =
      interpolate2_computevalue(init_grid_AB[x_min + 128][y_min + 128].first,
                                init_grid_AB[x_max + 128][y_min + 128].first,
                                init_grid_AB[x_min + 128][y_max + 128].first,
                                init_grid_AB[x_max + 128][y_max + 128].first,
                                x_min,
                                x_max,
                                y_min,
                                y_max,
                                x,
                                y);

  double b_interpolated =
      interpolate2_computevalue(init_grid_AB[x_min + 128][y_min + 128].second,
                                init_grid_AB[x_max + 128][y_min + 128].second,
                                init_grid_AB[x_min + 128][y_max + 128].second,
                                init_grid_AB[x_max + 128][y_max + 128].second,
                                x_min,
                                x_max,
                                y_min,
                                y_max,
                                x,
                                y);

  return std::make_pair(a_interpolated, b_interpolated);
}

inline void
initMatLABCH(std::vector< double > &init_grid_L,
             std::vector< double > &grid_L,
             std::vector< std::vector< std::pair< double, double > > >
                 &init_grid_AB) //, std::vector<std::vector<double>> & grid_AB)
{
  int size_L = 100001;
  int size_row = 257;
  int size_col = 257;
  init_grid_L.assign(size_L, 0);
  grid_L.assign(size_L, 0);

  int size_tabAB = 66049;
  init_grid_AB.assign(size_col,
                      std::vector< std::pair< double, double > >(
                          size_row, std::make_pair(0.0, 0.0)));

  //L_data
  //int cpt = 0; 
  for(int cpt = 0; cpt < size_L; cpt++)
    grid_L[cpt] = L_data[cpt];


  // init grid L (from matlab code)
  for(size_t i = 0; i < init_grid_L.size(); i++)
  {
    init_grid_L[i] = i * 0.001;
  }

  // inigt grid AB
  for(int cpt2 = 0; cpt2 < size_tabAB; cpt2++)
  {
    init_grid_AB[(int)(RegularGridInit_0_0_1[cpt2] + 128)][(int)(RegularGridInit_0_0_2[cpt2] + 128)].first = RegularGrid_0_0_1[cpt2];
    init_grid_AB[(int)(RegularGridInit_0_0_1[cpt2] + 128)][(int)(RegularGridInit_0_0_2[cpt2] + 128)].second = RegularGrid_0_0_2[cpt2];
  }
}


inline void
lab2000hl_conversion(
    double L,
    double A,
    double B,
    double *lab2000hl,
    std::vector< double > &init_grid_L,
    std::vector< double > &grid_L,
    std::vector< std::vector< std::pair< double, double > > > &init_grid_AB)
{
  lab2000hl[0] = interpolate1_process(init_grid_L, grid_L, L);
  std::pair< double, double > me_a_b = interpolate2_process(init_grid_AB, A, B);
  lab2000hl[1] = me_a_b.first;
  lab2000hl[2] = me_a_b.second;
}
