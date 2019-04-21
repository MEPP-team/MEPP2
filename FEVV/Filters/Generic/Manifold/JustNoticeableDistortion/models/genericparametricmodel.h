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
