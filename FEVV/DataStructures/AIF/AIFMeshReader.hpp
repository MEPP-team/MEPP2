#pragma once

#include <iostream>
#include <sstream>
#include <fstream>
#include <vector>
#include <boost/shared_ptr.hpp>

#include "FEVV/DataStructures/AIF/AIFMesh.hpp"

#include "FEVV/Tools/IO/MeshReaderInterface.h"

namespace FEVV {
namespace DataStructures {
namespace AIF {


/**
 * \class	AIFMeshReader
 * \brief	This class represents an AIFMesh object reader.
 *			An AIFMeshReader reads a mesh file (.obj, .off, .ply, .vtk,
 *.vtp, .vtu) and build the corresponding AIFMesh object. Note that some
 *advanced properties (e.g. multi-texturing) may have not yet been implemented.
 * \see		AIFMesh
 */
class AIFMeshReader : public MeshReaderInterface< AIFMesh >
{

public:
  typedef AIFMeshReader self;
  typedef MeshReaderInterface< AIFMesh > superclass;
  typedef superclass::output_type output_type;
  typedef superclass::ptr_output ptr_output;
  typedef boost::shared_ptr< self > ptr_reader;

  typedef AIFMesh::CoordinateType coord_type;        // position coordinate type
  typedef AIFMesh::NormalCoordinateType coordN_type; // normal coordinate type
  typedef AIFMesh::CoordinateType coordC_type;       // color coordinate type
  typedef AIFMesh::CoordinateType coordT_type;       // texture coordinate type
  typedef unsigned long index_type;

public:
  /*!
   * 			Read a mesh file
   * \param	filePath	Path to the input mesh file.
   * \return	The built AIFMesh mesh object.
   */
  ptr_output read(const std::string &filePath);
};

} // namespace AIF
} // namespace DataStructures
} // namespace FEVV


#include "FEVV/DataStructures/AIF/AIFMeshReader.inl"
