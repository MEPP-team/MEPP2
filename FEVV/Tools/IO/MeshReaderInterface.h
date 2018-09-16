#ifndef __MeshReaderInterface_h
#define __MeshReaderInterface_h

#include <string>
#include <boost/shared_ptr.hpp>
#include <boost/filesystem.hpp>

namespace FEVV {
namespace DataStructures {

template< typename TMeshOutput >
class MeshReaderInterface
{
public:
  // Basic type definitions.
  typedef TMeshOutput output_type;
  typedef boost::shared_ptr< output_type > ptr_output;
  typedef MeshReaderInterface< output_type > self;

  virtual ptr_output read(const std::string &) = 0;
};

} // namespace DataStructures
} // namespace FEVV

#endif
