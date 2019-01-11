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

#include <string>
#include <boost/shared_ptr.hpp>
#include <boost/filesystem.hpp>

namespace FEVV {
namespace DataStructures {

template< typename TMeshInput >
class MeshWriterInterface
{
public:
  // Basic type definitions.
  typedef TMeshInput input_type;
  typedef boost::shared_ptr< input_type > ptr_input;
  typedef MeshWriterInterface< input_type > self;

  virtual void write(const ptr_input, const std::string &) = 0;
};

} // namespace DataStructures
} // namespace FEVV

