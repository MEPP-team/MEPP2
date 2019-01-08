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

