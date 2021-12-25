#pragma once

#ifdef __GNUC__
// disable all GCC warnings in the rest of the current header
#pragma GCC system_header
#endif //__GNUC__

#undef ERROR

#include <draco/compression/decode.h>
#include <draco/io/mesh_io.h>
#include <draco/io/point_cloud_io.h>
#include <draco/core/encoder_buffer.h>
#include <draco/compression/bit_coders/rans_bit_decoder.h>
#include <draco/core/decoder_buffer.h>
#include <draco/core/varint_decoding.h>
#include <draco/compression/bit_coders/symbol_bit_decoder.h>

#include <draco/compression/entropy/rans_symbol_decoder.h>
#include <draco/mesh/mesh.h>
#include <draco/compression/entropy/symbol_decoding.h>

#include <draco/io/obj_encoder.h>
#include <draco/io/parser_utils.h>

#include <draco/attributes/point_attribute.h>
#include <draco/attributes/geometry_indices.h>
