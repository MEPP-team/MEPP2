#pragma once

#include <cmath>

#include "genericparametricmodel.h"

//-----------------------------------------------------------------------------------

/// implements the contrast masking model of Scott Daly
/// "Visible differences predictor: an algorithm for the assessment of image
/// fidelity." Human Vision, Visual Processing, and Digital Display III. Vol.
/// 1666. and used in eq. (7) "Just noticeable distortion profile for
/// flat-shaded 3D mesh surfaces." IEEE transactions on visualization and
/// computer graphics 22.11 (2016).
///
/// The model requires 4 parameters, and takes the contrast value as input

class DalyMasking
    : public GenericParametricModel< Eigen::Vector4d, double, double >
{

  //-----------------------------------------------------------------------------------

public:
  DalyMasking()
      : GenericParametricModel< ParameterType, InputType, OutputType >()
  {
    default_params();
  }

  DalyMasking(const ParameterType &param)
      : GenericParametricModel< ParameterType, InputType, OutputType >(param)
  {
  }

  virtual ~DalyMasking() {}

  //-----------------------------------------------------------------------------------

public:
  virtual void compute(const InputType &in, OutputType &out) const
  {
    out =
        pow(1. + pow(m_param(0) * pow(m_param(1) * in, m_param(2)), m_param(3)),
            1. / m_param(3));
  }

  //-----------------------------------------------------------------------------------

protected:
  /// initialises the daly masking model as described in:
  /// "Just noticeable distortion profile for flat-shaded 3D mesh surfaces."
  /// IEEE transactions on visualization and computer graphics 22.11 (2016).
  virtual void default_params()
  {
    m_param = ParameterType(0.0078, 88.29, 1.0, 4.207);
  }
};

