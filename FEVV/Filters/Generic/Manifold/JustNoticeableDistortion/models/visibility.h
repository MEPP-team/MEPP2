#ifndef JNDModels_VISIBILITY_H
#define JNDModels_VISIBILITY_H

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

#endif // VISIBLITY_H
