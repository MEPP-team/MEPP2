#ifndef _EXTRACT_VPROPRES
#define _EXTRACT_VPROPRES


#include <Eigen/Dense>

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
int
val_pro(const EigenMatrix &a, EigenvectorsType &u, EigenvalueType &ad);

/*
  Sort eigenvalues by descending order of their absolute values,
  and sort eigenvectors accordingly

  \param  AD  eigenvalues ; at function return the values are sorted
              by descending order of their absolute value
  \param  U   eigenvectors ; at function return the vectors are sorted
              in the same order as eigenvalues
*/
void
eig_srt(EigenvalueType &ad, EigenvectorsType &u);


#endif
