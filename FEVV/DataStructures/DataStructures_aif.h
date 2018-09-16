#ifndef FEVV_DATASTRUCTURES_AIF_H
#define FEVV_DATASTRUCTURES_AIF_H

#include "FEVV/DataStructures/AIF/AIFMesh.hpp"
#include "FEVV/DataStructures/AIF/AIFMeshHelpers.h"
#include "FEVV/DataStructures/AIF/AIFMeshReader.hpp"

namespace FEVV {

using MeshAIF = FEVV::DataStructures::AIF::AIFMesh;
using MeshAIFPtr = FEVV::DataStructures::AIF::AIFMeshReader::ptr_output;

} // namespace FEVV


#endif // FEVV_DATASTRUCTURES_AIF_H
