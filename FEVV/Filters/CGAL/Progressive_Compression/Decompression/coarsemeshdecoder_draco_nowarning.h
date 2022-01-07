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

#ifdef __GNUC__
// disable all GCC warnings in the rest of the current header
#pragma GCC system_header
#endif //__GNUC__

#include <draco/compression/decode.h>
#include <draco/compression/bit_coders/rans_bit_decoder.h>
#include <draco/core/decoder_buffer.h>
#include <draco/core/varint_decoding.h>
#include <draco/compression/bit_coders/symbol_bit_decoder.h>
#include <draco/mesh/mesh.h>
