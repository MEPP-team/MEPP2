if( NOT BUILD_DOCUMENTATION )
  return()
endif()

add_subdirectory( Doxygen )

### Syntactic sugar for limiting CI load. A CI platform wishing/needing to only
# build the documentation, but that does not need the binaries (libraries and
# executable) can use this target ("make doc").
# Note: for the time being, the documentation being limited to the sole
# "DocDoxygen" sub-target, the "doc" target acts as an alias of the "DocDoxygen"
# target. This will change when adding new documentation sub-targets...
add_custom_target( doc ) 
add_dependencies ( doc DocDoxygen )
