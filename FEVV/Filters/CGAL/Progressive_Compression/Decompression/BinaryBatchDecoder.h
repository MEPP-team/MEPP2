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

#include <boost/graph/graph_traits.hpp>
#include <boost/graph/properties.hpp>
#include "FEVV/Wrappings/Geometry_traits.h"

#include <list>
#include <vector>
#include <utility>

#ifdef _MSC_VER
#pragma warning(push)
#pragma warning(disable : 4146 26812 26451)
#endif
#include "BinaryBatchDecoder_draco_nowarning.h"
#if defined _MSC_VER
#pragma warning(pop)
#endif

namespace FEVV {
namespace Filters {

/// Class used to decode binary data from draco: bitsmasks, geometric and 
/// attribute info.
template<
typename HalfedgeGraph,
typename PointMap,
typename Vector = typename FEVV::Geometry_traits< HalfedgeGraph >::Vector >
class BinaryBatchDecoder
{
public:
BinaryBatchDecoder() {};
~BinaryBatchDecoder() {};
BinaryBatchDecoder(draco::DecoderBuffer &buffer, int bit_quantization)
    : _buffer(buffer)
{
  _bit_quantization = bit_quantization;
};



void DecodeBitmask(std::list<bool> &bitmask)
{
	int size_bitmask;
	draco::RAnsBitDecoder decoder;
	draco::DecodeVarint(&size_bitmask,&_buffer);
    
	decoder.StartDecoding(&_buffer);
  for(int i = 0; i < size_bitmask; i++)
  {
    bool next = decoder.DecodeNextBit();
    bitmask.push_back(next);
  }
  decoder.EndDecoding();
}

void DecodeResiduals(std::list<std::vector<Vector>> &residuals, int nb_residuals)
{
	int size_residuals;
	draco::DecodeVarint(&size_residuals,&_buffer);
	
	draco::SymbolBitDecoder symbolBitDecoder;
	
    std::vector<uint32_t> residuals_uint;	
	int nb = size_residuals * 3 * nb_residuals;
    residuals_uint.reserve(nb);
	
	symbolBitDecoder.StartDecoding(&_buffer);
	for (int i = 0; i < nb; i++)
	{
      uint32_t temp;
      symbolBitDecoder.DecodeLeastSignificantBits32(_bit_quantization + 1, &temp);
      residuals_uint.push_back(temp);
	}
	symbolBitDecoder.EndDecoding();
	
	std::vector<int32_t> residuals_int;	
    residuals_int.resize(nb);
	draco::ConvertSymbolsToSignedInts(residuals_uint.data(), nb, residuals_int.data());
	
	for (int i = 0; i < size_residuals; i++)
	{
		std::vector<Vector> paire;
		paire.reserve(nb_residuals);
		for (int j = 0; j < nb_residuals; ++j)
		{
			paire.push_back(Vector(residuals_int[3 * (nb_residuals*i + j) + 0], 
                                   residuals_int[3 * (nb_residuals*i + j) + 1], 
                                   residuals_int[3 * (nb_residuals*i + j) + 2]));
		}

		residuals.push_back(std::move(paire));
	}
}

	private:
        draco::DecoderBuffer &_buffer;
        int _bit_quantization;
};


}
}
