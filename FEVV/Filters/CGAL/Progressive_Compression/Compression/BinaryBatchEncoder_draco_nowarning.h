#pragma once

#ifdef __GNUC__
// disable all GCC warnings in the rest of the current header
#pragma GCC system_header
#endif //__GNUC__


#include <draco/compression/encode.h>
#include <draco/io/mesh_io.h>
#include <draco/io/point_cloud_io.h>
#include <draco/core/encoder_buffer.h>
#include <draco/compression/bit_coders/rans_bit_encoder.h>
#include <draco/core/decoder_buffer.h>
#include <draco/core/varint_encoding.h>
#include <draco/compression/bit_coders/symbol_bit_encoder.h>
#include <draco/compression/entropy/shannon_entropy.h>

#include <draco/compression/entropy/rans_symbol_encoder.h>
#include <draco/mesh/mesh.h>
#include <draco/compression/entropy/symbol_encoding.h>
#include <draco/mesh/triangle_soup_mesh_builder.h>

#include <draco/io/obj_encoder.h>

