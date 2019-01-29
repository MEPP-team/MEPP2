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

#include "genericparametricmodel.h"

#include "contrastmasking.h"
#include "contrastsensitivity.h"

//-----------------------------------------------------------------------------------

struct NWHWD16_Param
{
  SarkisonCSF csf;
  DalyMasking masking;
};

/// implements the threshold mode in [ eq. (8) ]:
/// "Just noticeable distortion profile for flat-shaded 3D mesh surfaces." IEEE
/// transactions on visualization and computer graphics 22.11 (2016).
///
/// The model takes:
/// - a Sarkison CSF and Daly Masking models as parameters
/// - contrast and frequency values as input
/// Output:
/// - a contrast threshold

class NWHWD16_Threshold
    : public GenericParametricModel< NWHWD16_Param, Eigen::Vector2d, double >
{

  //-----------------------------------------------------------------------------------

public:
  NWHWD16_Threshold()
      : GenericParametricModel< ParameterType, InputType, OutputType >()
  {
    default_params();
  }

  NWHWD16_Threshold(const ParameterType &param)
      : GenericParametricModel< ParameterType, InputType, OutputType >(param)
  {
  }

  virtual ~NWHWD16_Threshold() {}

  //-----------------------------------------------------------------------------------

public:
  virtual void compute(const InputType &in, OutputType &out) const
  {
    double c = in(0);
    double f = in(1);
    double _csf = m_param.csf(f);
    out = m_param.masking(c * _csf) / _csf;
  }

  //-----------------------------------------------------------------------------------

protected:
  /// initialise the default params as in
  /// "Just noticeable distortion profile for flat-shaded 3D mesh surfaces."
  /// IEEE transactions on visualization and computer graphics 22.11 (2016).
  virtual void default_params()
  {
    m_param.csf = SarkisonCSF(SarkisonCSF::ParameterType(-15.13, 0.0096, 0.64));
    m_param.masking =
        DalyMasking(DalyMasking::ParameterType(0.0078, 88.29, 1.0, 4.207));
  }
};

//=================================================================================//
