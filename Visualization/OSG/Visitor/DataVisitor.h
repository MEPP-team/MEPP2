// Copyright (c) 2012-2019 University of Lyon and CNRS (France).
// All rights reserved.
//
// This file is part of MEPP2; you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published 
// by the Free Software Foundation; either version 3 of the License, 
// or (at your option) any later version.
//
// This file is provided AS IS with NO WARRANTY OF ANY KIND, INCLUDING THE
// WARRANTY OF DESIGN, MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE.
#pragma once

#include <osg/Geode>
#include <osg/NodeVisitor>

#include <iostream>

#include "Visualization/Helpers/DataStructures.h"
#include "Visualization/BaseViewerOSG.h"

namespace FEVV {

class BaseViewerOSG;

class DataVisitor : public osg::NodeVisitor
{
public:
  typedef Helpers::Model< osg::Node > Data;
  typedef std::vector< Data > Output;

public:
  DataVisitor(FEVV::BaseViewerOSG *_viewer)
      : NodeVisitor(NodeVisitor::TRAVERSE_ALL_CHILDREN), current_level(0),
        exportVector(new Output()), viewer(_viewer)
  {
    // add for time mode - for geodes with setNodeMask(0x0) - MT 14/02/18

    // doc :
    // http://podsvirov.github.io/osg/reference/openscenegraph/a00553.html#a99be1e526672165ae46af189ed7f472c
    //
    // void osg::NodeVisitor::setNodeMaskOverride  (   Node::NodeMask    mask  )
    //   - Set the NodeMaskOverride mask.
    //   - Used in validNodeMask() to determine whether to operate on a node or
    //   its subgraph, by OR'ing NodeVisitor::_nodeMaskOverride with the Node's
    //   own Node::_nodeMask.
    //     Typically used to force on nodes which may have been switched off by
    //     their own Node::_nodeMask.
    //
    // void osg::NodeVisitor::setTraversalMask   (   Node::NodeMask    mask  )
    //   - Set the TraversalMask of this NodeVisitor.
    //   - The TraversalMask is used by the NodeVisitor::validNodeMask() method
    //   to determine whether to operate on a node and its subgraph.
    //     validNodeMask() is called automatically in the Node::accept() method
    //     before any call to NodeVisitor::apply(), apply() is only ever called
    //     if validNodeMask returns true. Note, if NodeVisitor::_traversalMask
    //     is 0 then all operations will be switched off for all nodes. Whereas
    //     setting both _traversalMask and _nodeMaskOverride to 0xffffffff will
    //     allow a visitor to work on all nodes regardless of their own
    //     Node::_nodeMask state.

    setNodeMaskOverride(0xffffffff); // simply works with that one
    // setTraversalMask(0xffffffff);
  }

  virtual ~DataVisitor()
  {
    exportVector->clear();
    delete exportVector;
  }

  virtual void apply(osg::Node &_node)
  {
    // std::string name = "";
    // for( unsigned ii = current_level; ii > 0; --ii )
    // {
    //   name += '\t';
    // }
    // name += _node.className();
    // name += ": ";
    // name += _node.getName();
    // exportVector.push_back( StringNodePair(name, &_node) );

    /// We remove from the list "debug" draws like the gizmo, the unit grid...
    if((_node.getName().size() > 0) &&
       ((_node.getName() == "Gizmo") || (_node.getName() == "UnitGrid") ||
        (_node.getName() == "Hud") ||
        (_node.getName().compare(0, 7, "Dragger") == 0)))
    {
      return;
    }

    Helpers::Model< osg::Node > result;
    result.node = &_node;
    result.name = _node.getName();
    result.type = Helpers::DataType::GROUP;
    result.position = current_level;
    result.viewer = viewer;
    exportVector->push_back(result);

    // std::cout << name << std::endl;
    ++current_level;
    traverse(_node);
    --current_level;
  }

  virtual void apply(osg::Geode &_geode)
  {
    // std::string name = "";
    // for( unsigned ii = current_level; ii > 0; --ii )
    // {
    //   name += '\t';
    // }
    // name += _geode.className();
    // name += ": ";
    // name += _geode.getName();
    // exportVector.push_back( StringNodePair(name, &_geode) );

    Helpers::Model< osg::Node > result;
    result.node = &_geode;
    result.name = _geode.getName();
    result.type = Helpers::DataType::MODEL;
    result.position = current_level;
    result.viewer = viewer;
    exportVector->push_back(result);

    // std::cout << name << std::endl;
  }

  Output *exportResults() { return exportVector; }

  void reset() { exportVector->clear(); }

private:
  unsigned int current_level;
  Output *exportVector;
  BaseViewerOSG *viewer;
};

} // namespace FEVV
