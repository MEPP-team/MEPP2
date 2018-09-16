#ifndef __MeshWriterInterface_h
#define __MeshWriterInterface_h

#include <string>
#include <boost/shared_ptr.hpp>
#include <boost/filesystem.hpp>

namespace FEVV {
namespace DataStructures {

template< typename TMeshInput >
class MeshWriterInterface
{
public:
  // Basic type definitions.
  typedef TMeshInput input_type;
  typedef boost::shared_ptr< input_type > ptr_input;
  typedef MeshWriterInterface< input_type > self;

  virtual void write(const ptr_input, const std::string &) = 0;
};

} // namespace DataStructures
} // namespace FEVV

#endif
