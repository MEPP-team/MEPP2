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

#include <Eigen/Dense>
#include <Eigen/Eigenvalues>
#include <cmath> // for std::abs()
#include <vector>


using EigenMatrix = Eigen::Matrix3d;
using EigenvectorsType = Eigen::Matrix3d;
using EigenvalueType = Eigen::Vector3d;

/*
  Compute eigenvalues and eigenvectors of a double symmetric matrix

  \param  A  (input) double symmetric matrix
  \param  U  (ouput) eigenvectors
  \param  AD (ouput) eigenvalues, unsorted
  \return  0 if ok
*/
inline
int
eigen_val_vect_compute(const EigenMatrix &a, EigenvectorsType &u, EigenvalueType &ad)
{
  Eigen::EigenSolver< EigenMatrix > solver(a, true);

  // check for error
  if( solver.info() != Eigen::ComputationInfo::Success )
    return -1;

  // ELO-note:
  //  solver.eigenvalues() and solver.eigenvectors() return
  //  matrices of std::complex values ; here we keep only
  //  the real part of this values
  ad = solver.eigenvalues().real();
  u = solver.eigenvectors().real();

  return 0;
}

/*
  Sort eigenvalues by descending order of their absolute values,
  and sort eigenvectors accordingly

  \param  AD  eigenvalues ; at function return the values are sorted
              by descending order of their absolute value
  \param  U   eigenvectors ; at function return the vectors are sorted
              in the same order as eigenvalues
*/
inline
void
eigen_val_vect_sort(EigenvalueType &ad, EigenvectorsType &u)
{
  // sort eigenvalues indices by descending order of their absolute values
  std::vector< size_t > sorted_indices;
  sorted_indices.push_back(0);
  for(int i = 1; i < 3; i++)
  {
    // find the place where to insert current eigenvalue
    auto it = sorted_indices.begin();
    for(; it != sorted_indices.end(); ++it)
    {
      if(std::abs(ad(*it)) < std::abs(ad(i)))
        break;
    }

    // insert current eigenvalue index here
    sorted_indices.insert(it, i);
  }

  // sort eigenvalues and eigenvectors
  EigenvalueType sorted_eigen_val;
  EigenvectorsType sorted_eigen_vect;

  for(int i = 0; i < 3; i++)
  {
    sorted_eigen_val(i) = ad(sorted_indices[i]);
    sorted_eigen_vect.col(i) = u.col(sorted_indices[i]);
  }

  ad = sorted_eigen_val;
  u = sorted_eigen_vect;
}
