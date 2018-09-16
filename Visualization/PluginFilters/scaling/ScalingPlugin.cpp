#include "ScalingPlugin.h"

#include <QtPlugin>

#if QT_VERSION < 0x050000 // (for QT4)
Q_EXPORT_PLUGIN2(ScalingPlugin, FEVV::ScalingPlugin)
#endif
