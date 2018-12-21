#pragma once

#include <vector>
#include <Eigen/Core>

template< typename TParam, typename TIn = double, typename TOut = double >
class GenericParametricModel
{

  //-----------------------------------------------------------------------------------

public:
  typedef TParam ParameterType;
  typedef TIn InputType;
  typedef TOut OutputType;

  //-----------------------------------------------------------------------------------

public:
  GenericParametricModel()
  {
    //    default_params();
  }

  GenericParametricModel(const ParameterType &param) : m_param(param) {}

  virtual ~GenericParametricModel() {}

  //-----------------------------------------------------------------------------------

public:
  void setParameters(const ParameterType &param) { m_param = param; }

  const ParameterType &param() const { return m_param; } // read only
  ParameterType &param() { return m_param; }

  //-----------------------------------------------------------------------------------

public:
  virtual void compute(const InputType &in, OutputType &out) const = 0;

  //-----------------------------------------------------------------------------------

public:
  OutputType compute(const InputType &in) const
  {
    OutputType out;
    compute(in, out);

    return out;
  }

  void compute(const std::vector< InputType > &in,
               std::vector< OutputType > &out) const
  {
    out.clear();

    for(const InputType i : in)
      out.push_back(compute(i));
  }

  OutputType operator()(const InputType &in) const { return compute(in); }

  void operator()(const InputType &in, OutputType &out) const
  {
    out = compute(in);
  }

  //-----------------------------------------------------------------------------------

protected:
  /// sets default values to m_params
  virtual void default_params() = 0;

  //-----------------------------------------------------------------------------------

protected:
  ParameterType m_param;
};
