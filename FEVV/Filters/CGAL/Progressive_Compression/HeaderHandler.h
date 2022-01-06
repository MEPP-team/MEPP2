// Copyright (c) 2012-2022 University of Lyon and CNRS (France).
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
#include<iostream>
#include "FEVV/Filters/CGAL/Manifold/Progressive_Compression/Utils/Parameters.h"

namespace FEVV {
namespace Filters {

class HeaderHandler // Manage file headers (compression and decompression)
{
public:
  HeaderHandler(); 
  ~HeaderHandler(); 

private:
	VKEPT_POSITION _vkept;
	PREDICTION_TYPE _predictor;

	void encode_header(std::string filepath) {

	}

	void decode_header() 
	{

	}
};
}
}