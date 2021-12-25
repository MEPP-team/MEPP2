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
