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