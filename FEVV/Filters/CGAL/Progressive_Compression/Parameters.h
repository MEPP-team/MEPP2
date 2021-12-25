#pragma once
#ifndef PARAMETERS_CPP
#define PARAMETERS_CPP
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

// FEVV::Filters::PREDICTION_TYPE
// StringToPredictionType(std::string prediction)
//{
//  if(prediction.compare("position") == 0)
//    return PREDICTION_TYPE::POSITION;
//  if(prediction.compare("delta") == 0)
//    return PREDICTION_TYPE::DELTA;
//  if(prediction.compare("barycentre_predictor") == 0)
//    return PREDICTION_TYPE::BARYCENTER_PREDICTOR;
//  if(prediction.compare("barycentre") == 0)
//    return PREDICTION_TYPE::BARYCENTER;
//  if(prediction.compare("simple_barycentre_half_ring") == 0)
//    return PREDICTION_TYPE::HALFRING;
//  if(prediction.compare("butterfly") == 0)
//    return PREDICTION_TYPE::BUTTERFLY;
//  return PREDICTION_TYPE::BUTTERFLY; // cas par défaut
//}
//
// FEVV::Filters::VKEPT_POSITION
// StringToVkept(std::string vkept)
//{
//  if(vkept.compare("halfedge") == 0)
//    return VKEPT_POSITION::HALFEDGE;
//  if(vkept.compare("midpoint") == 0)
//    return VKEPT_POSITION::MIDPOINT;
//  if(vkept.compare("fulledge") == 0)
//    return VKEPT_POSITION::FULLEDGE;
//  return VKEPT_POSITION::MIDPOINT;
//}
// FEVV::Filters::METRIC_TYPE
// StringToMetricType(std::string metric)
//{
//  if(metric.compare("qem") == 0)
//    return METRIC_TYPE::QEM;
//  if(metric.compare("no_metric") == 0)
//    return METRIC_TYPE::NO_METRIC;
//  if(metric.compare("edge_length") == 0)
//    return METRIC_TYPE::EDGE_LENGTH;
//  if(metric.compare("volume_preserving") == 0)
//    return METRIC_TYPE::VOLUME_PRESERVING;
//  return METRIC_TYPE::NO_METRIC;
//}
} // namespace Filters
} // namespace FEVV
#endif
