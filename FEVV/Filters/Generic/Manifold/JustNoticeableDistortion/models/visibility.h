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
#include "psychometricfunction.h"

//-----------------------------------------------------------------------------------

/// implements the visibility model, i.e., eq. (11) in:
/// "Just noticeable distortion profile for flat-shaded 3D mesh surfaces." IEEE
/// transactions on visualization and computer graphics 22.11 (2016).
///
/// The model takes:
/// - a WeibulPsychometricFunction as a parameter
/// - the change in contrast and threshold as input
/// and outputs:
/// - a visibility probability between 0 and 1

class VisibilityModel
    : public GenericParametricModel< WeibulPsychometricFunction,
                                     Eigen::Vector2d,
                                     double >
{

  //-----------------------------------------------------------------------------------

public:
  VisibilityModel()
      : GenericParametricModel< ParameterType, InputType, OutputType >()
  {
    default_params();
  }

  VisibilityModel(const ParameterType &param)
      : GenericParametricModel< ParameterType, InputType, OutputType >(param)
  {
  }

  virtual ~VisibilityModel() {}

  //-----------------------------------------------------------------------------------

public:
  virtual void compute(const InputType &in, OutputType &out) const
  {
    double dc = in(0); // change in contrast
    double T = in(1);  // threshold

    if(T == 0.)
      out = OutputType(dc == T);

    out = m_param(dc / T);
  }

  //-----------------------------------------------------------------------------------

protected:
  /// initialise the weibul distribution with beta = 3.5;
  virtual void default_params() { m_param = ParameterType(3.5); }
};
