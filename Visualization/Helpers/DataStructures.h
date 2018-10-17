#pragma once

#if defined(VisualizationDataStructures_RECURSES)
#error Recursive header files inclusion detected in VisualizationDataStructures.h
#else // defined(VisualizationDataStructures_RECURSES)
/** Prevents recursive inclusion of headers. */
#define VisualizationDataStructures_RECURSES

#if !defined VisualizationDataStructures_h
/** Prevents repeated inclusion of headers. */
#define VisualizationDataStructures_h

#include <QMetaType>
#include "Base/Color.hpp"
#include "Visualization/BaseViewerOSG.h"
#include <osg/Node>


namespace FEVV {

class BaseViewerOSG;

namespace Helpers {

enum class DataType : unsigned char { MODEL = 1, GROUP = 2, DEBUG_OBJECT = 3 };

template< typename NodeType >
struct Model
{
  NodeType *node;
  std::string name;
  DataType type;
  unsigned int position;
  BaseViewerOSG *viewer;
};

} // namespace Helpers

} // namespace FEVV

Q_DECLARE_METATYPE(FEVV::Helpers::Model< osg::Node >)

#endif // !defined VisualizationDataStructures_h

#undef VisualizationDataStructures_RECURSES
#endif // else defined(VisualizationDataStructures_RECURSES)
