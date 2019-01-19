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

#include <cmath>

#include "genericparametricmodel.h"

//-----------------------------------------------------------------------------------

/// implements a weibul psychometric function
/// 1. - exp( -pow(x, beta) )
///
/// The model needs 1 parameter (beta),
///           takes 1 input
///       and outputs a value between 0 and 1

class WeibulPsychometricFunction
    : public GenericParametricModel< double, double, double >
{

  //-----------------------------------------------------------------------------------

public:
  WeibulPsychometricFunction()
      : GenericParametricModel< ParameterType, InputType, OutputType >()
  {
    default_params();
  }

  WeibulPsychometricFunction(const ParameterType &param)
      : GenericParametricModel< ParameterType, InputType, OutputType >(param)
  {
  }

  virtual ~WeibulPsychometricFunction() {}

  //-----------------------------------------------------------------------------------

public:
  virtual void compute(const InputType &in, OutputType &out) const
  {
    out = 1. - std::exp(-std::pow(in, m_param));
  }

  //-----------------------------------------------------------------------------------

protected:
  /// initialise beta to 3.5
  virtual void default_params() { m_param = 3.5; }
};
