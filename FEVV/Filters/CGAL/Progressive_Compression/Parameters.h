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

#include <iostream>

namespace FEVV {
namespace Filters {

enum class PREDICTION_TYPE {
  POSITION = 0, /// encode the 2 edge vertex positions
  DELTA, /// encode the delta between the 2 edge vertex positions
  //BARYCENTER_PREDICTOR, // (not implemented in this release)
  //BARYCENTER, // (not implemented in this release)
  //HALFRING, // (not implemented in this release)
  BUTTERFLY, /// encode the delta between the predicted delta and the delta between the 2 edge vertex positions
};

enum class VKEPT_POSITION { MIDPOINT = 0, /// midpoint position
                            HALFEDGE, /// edge target position
                            FULLEDGE /// an "optimal" position along the collapsed edge (not implemented in this release)
                           };

enum class METRIC_TYPE {
  NO_METRIC = 0,
  QEM_3D, /// memoryless variant of QEM error (Quadric Error Metric)
  EDGE_LENGTH, /// edge length metric
  VOLUME_PRESERVING /// local absolute volume error metric
};


/**
 * \brief Parameters contains the compression parameters except the
 *        stopping criteria. 
 */
class Parameters // Local Parameters
{
private:
  PREDICTION_TYPE _prediction;
  VKEPT_POSITION _vkept_pos;
  METRIC_TYPE _metric;
  bool _first_call;       
  bool _allow_duplicates; // Global parameters
  int _quantization_bits;
  
public:
  PREDICTION_TYPE get_prediction() const { return _prediction; } 
  VKEPT_POSITION get_vkept_position() const { return _vkept_pos; }
  METRIC_TYPE get_metric() const { return _metric; }
  bool get_first_call() const { return _first_call; }
  bool get_allow_duplicates() const { return _allow_duplicates; }
  Parameters(PREDICTION_TYPE __prediction = PREDICTION_TYPE::BUTTERFLY,
             VKEPT_POSITION __vkept_pos = VKEPT_POSITION::MIDPOINT,
             METRIC_TYPE __metric = METRIC_TYPE::EDGE_LENGTH,
             bool __first_call = true,
             bool __allow_duplicates = false,
             int __quantization_bits = 16)
  {
    _prediction = __prediction;
    _vkept_pos = __vkept_pos;
    _metric = __metric;
    _first_call = __first_call;
    _allow_duplicates = __allow_duplicates;
    _quantization_bits = __quantization_bits;
  }

  ~Parameters() {}
  void set_first_call(bool _first) { _first_call = _first; }
  int get_quantization() const { return _quantization_bits; }

  std::string get_prediction_type_as_string() const
  {
    switch(_prediction)
    {
    case PREDICTION_TYPE::POSITION:
      return "position";
      break;
    case PREDICTION_TYPE::DELTA:
      return "delta";
      break;
    /*case PREDICTION_TYPE::BARYCENTER_PREDICTOR:
      return "barycentre_predictor";
      break;
    case PREDICTION_TYPE::BARYCENTER:
      return "barycentre";
      break;
    case PREDICTION_TYPE::HALFRING:
      return "simple_barycentre_half_ring";
      break;*/
    case PREDICTION_TYPE::BUTTERFLY:
      return "butterfly";
      break;
    default:
      return "unknown";
      break;
    }
  }
  std::string get_vkept_type_as_string() const
  {
    switch(_vkept_pos)
    {
    case VKEPT_POSITION::HALFEDGE:
      return "halfedge";
      break;
    case VKEPT_POSITION::FULLEDGE:
      return "fulledge";
      break;
    case VKEPT_POSITION::MIDPOINT:
      return "midpoint";
      break;
    default:
      return "midpoint";
      break;
    }
  }
  std::string get_metric_type_as_string() const
  {
    switch(_metric)
    {
    case FEVV::Filters::METRIC_TYPE::NO_METRIC:
      return "no_metric";
      break;
    case FEVV::Filters::METRIC_TYPE::QEM_3D:
      return "qem";
      break;
    case FEVV::Filters::METRIC_TYPE::EDGE_LENGTH:
      return "edge-length";
      break;
    case FEVV::Filters::METRIC_TYPE::VOLUME_PRESERVING:
      return "volume-preserving";
      break;
    default:
      return "wut"; // On doit pas tomber là dessus normalement
      break;
    }
  }

  void set_methods(METRIC_TYPE metr, VKEPT_POSITION vk)
  {
    _metric = metr;
    _vkept_pos = vk;
  }
};

} // namespace Filters
} // namespace FEVV
