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
  POSITION = 0,
  DELTA,
  BARYCENTER_PREDICTOR,
  BARYCENTER,
  HALFRING,
  BUTTERFLY,
};

enum class VKEPT_POSITION { MIDPOINT = 0, HALFEDGE, FULLEDGE };

enum class METRIC_TYPE {
  NO_METRIC = 0,
  QEM,
  EDGE_LENGTH,
  VOLUME_PRESERVING
};


/*
\brief Parameters: contains the methods used by the compression
*/

class Parameters // Local Parameters
{
private:
  PREDICTION_TYPE Prediction;
  VKEPT_POSITION Vkept_pos;
  METRIC_TYPE Metric;
  bool FirstCall;       
  bool AllowDuplicates; // Global parameters
  int QuantizationBits;
  

public:
  PREDICTION_TYPE getPrediction() const { return Prediction; } 
  VKEPT_POSITION getVKeptPosition() const { return Vkept_pos; }
  METRIC_TYPE getMetric() const { return Metric; }
  bool getFirstCall() const { return FirstCall; }
  bool getAllowDuplicates() const { return AllowDuplicates; }
  Parameters(PREDICTION_TYPE _Prediction = PREDICTION_TYPE::BUTTERFLY,
             VKEPT_POSITION _Vkept_pos = VKEPT_POSITION::MIDPOINT,
             METRIC_TYPE _Metric = METRIC_TYPE::EDGE_LENGTH,
             bool _FirstCall = true,
             bool _AllowDuplicates = false,
             int _QuantizationBits = 16)
  {
    Prediction = _Prediction;
    Vkept_pos = _Vkept_pos;
    Metric = _Metric;
    FirstCall = _FirstCall;
    AllowDuplicates = _AllowDuplicates;
    QuantizationBits = _QuantizationBits;
  }

  ~Parameters() {}
  void setFirstCall(bool _first) { FirstCall = _first; }
  int getQuantization() const { return QuantizationBits; }

  std::string PredictionTypeToString() const
  {
    switch(Prediction)
    {
    case PREDICTION_TYPE::POSITION:
      return "position";
      break;
    case PREDICTION_TYPE::DELTA:
      return "delta";
      break;
    case PREDICTION_TYPE::BARYCENTER_PREDICTOR:
      return "barycentre_predictor";
      break;
    case PREDICTION_TYPE::BARYCENTER:
      return "barycentre";
      break;
    case PREDICTION_TYPE::HALFRING:
      return "simple_barycentre_half_ring";
      break;
    case PREDICTION_TYPE::BUTTERFLY:
      return "butterfly";
      break;
    default:
      return "butterfly";
      break;
    }
  }
  std::string VKeptToString() const
  {
    switch(Vkept_pos)
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
  std::string MetricToString() const
  {
    switch(Metric)
    {
    case FEVV::Filters::METRIC_TYPE::NO_METRIC:
      return "no_metric";
      break;
    case FEVV::Filters::METRIC_TYPE::QEM:
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

  void ChangeMethods(METRIC_TYPE metr, VKEPT_POSITION vk)
  {
    Metric = metr;
    Vkept_pos = vk;
  }
};

} // namespace Filters
} // namespace FEVV
