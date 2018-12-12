find_package( Doxygen )

if( NOT DOXYGEN_FOUND)
  return()
endif()

if( DOXYGEN_DOT_FOUND )
  set( DOXYGEN_HAVE_DOT "YES" )
else()
  set( DOXYGEN_HAVE_DOT "NO" )
endif()

# we run several doxygen builds in parallel, so we reduce the number of procs
# used by each doxygen to run dot
if(NOT CMAKE_VERSION VERSION_LESS 2.8.5)
  include(ProcessorCount)
  processorcount(proc_count)
  math(EXPR DOT_JOBS "${proc_count} / 4")
  if(DOT_JOBS LESS 1)
    set(DOT_JOBS 1)
  endif()
else()
  set(DOT_JOBS 2)
endif()

configure_file(
  ${CMAKE_CURRENT_SOURCE_DIR}/Doxyfile.in
  ${CMAKE_CURRENT_BINARY_DIR}/Doxyfile
  @ONLY)

# Because "@includedoc Install.md" doesn't work for including markdown
# (the markdown is indeed included but isn't interpreted properly, refer to
# http://doxygen.10944.n7.nabble.com/Using-includedoc-with-markdown-td7719.html
# we need to "manually" construct the doxygen wrapper of the Install.md
file(WRITE  ${CMAKE_CURRENT_BINARY_DIR}/InstallWrapper.dox
  "/**")
file(APPEND ${CMAKE_CURRENT_BINARY_DIR}/InstallWrapper.dox
  "@page InstallWrapperPage FEVV installation instructions\r")
file(READ ${PROJECT_SOURCE_DIR}/Install.md CONTENTS)
file(APPEND ${CMAKE_CURRENT_BINARY_DIR}/InstallWrapper.dox
  "${CONTENTS}")
file(APPEND ${CMAKE_CURRENT_BINARY_DIR}/InstallWrapper.dox
  "**/")

# Construct the list of source dependencies
# FIXME: "factorize" the list of patterns used bellow with the content
# of FILE_PATTERNS within Doxygen configuration variable
set(DEPENDENCIES)
file(GLOB_RECURSE files
  *.h
  *.hxx
  *.cpp
  *.dox
  *.md)
list(APPEND DEPENDENCIES ${files})

add_custom_command(
  OUTPUT ${CMAKE_CURRENT_BINARY_DIR}/html/index.html
  COMMAND ${DOXYGEN} ${CMAKE_CURRENT_BINARY_DIR}/Doxyfile
  COMMENT "Generating Doxygen html output..."
  DEPENDS ${CMAKE_CURRENT_BINARY_DIR}/Doxyfile ${DEPENDENCIES}
)

add_custom_target(
  DocDoxygen ALL
  DEPENDS ${CMAKE_CURRENT_BINARY_DIR}/html/index.html
)

### Someday there will be an installation stage:
# install(
#  DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/html/
#  DESTINATION ${INSTALL_DATA_DIR}/Doxygen/
#  COMPONENT Documentation
# )
