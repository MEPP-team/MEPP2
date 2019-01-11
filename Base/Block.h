// Copyright (c) 2012-2019 University of Lyon and CNRS (France).
// All rights reserved.
//
// This file is part of MEPP2; you can redistribute it and/or modify
// it under the terms of the GNU Lesser General Public License as
// published by the Free Software Foundation; either version 3 of
// the License, or (at your option) any later version.
//
// This file is provided AS IS with NO WARRANTY OF ANY KIND, INCLUDING THE
// WARRANTY OF DESIGN, MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE.
#pragma once

#include <iostream>
#include <chrono>
#include <ctime>

namespace FEVV {

/**
 * class Block
 * \brief This class may be used to draw blocks in your code. It will show the
 * elapsed time between the begin and the end of the block. Many blocks can be
 * used in the same time, differentiate by a unique identifier (string).
 *
 * The code below shows a possible use of this class.
 * @code
 * #include "Base/Block.h"
 *
 * using namespace FEVV;
 * // ...
 * {
 *
 *   Block::begin("my_example", "Example of Block usage.");
 *   // ... some computation
 *   Block::end("my_example");
 * }
 * @endcode
 *
 */
class Block
{
public:
  /**
   * Default constructor.
   * @note Disable by default.
   */
  Block() = delete;

  /**
   * Copy constructor.
   * @note Disable by default.
   */
  Block(const Block &_Block) = delete;

  /**
   * Destructor.
   * @note Disable by default.
   */
  ~Block() = delete;

  /**
   * Begin a block.
   *
   * @param[in]    _uid        a unique identifier to distinguish blocks.
   * @param[in]    _message    a message to draw in the standard output.
   */
  static void begin(const std::string &_uid, const std::string &_message)
  {
    Block::InnerData data;
    data.message = _message;
    data.startClock = std::chrono::system_clock::now();
    Block::blocks.insert({_uid, data});

#if(defined(UNIX))
    std::cout << "\033[0m\033[1m\033[32m"; // Bold + green color
#endif
    std::cout << "################";
    for(unsigned int i = 0; i < _message.length(); ++i)
    {
      std::cout << "#";
    }
    std::cout << "#" << std::endl;
    std::cout << "# Begin block: " << _message << " #" << std::endl;
    std::cout << "################";
    for(unsigned int i = 0; i < _message.length(); ++i)
    {
      std::cout << "#";
    }
    std::cout << "#" << std::endl;
#if(defined(UNIX))
    std::cout << "\033[0m"; // reset
#endif
  }

  /**
   * End a block.
   *
   * @param[in]    _uid        a unique identifier to distinguish blocks.
   */
  static void end(const std::string &_uid)
  {
    if(Block::blocks.find(_uid) == Block::blocks.end())
    {
      std::cerr << "Can't found a block with this uid: " << _uid << std::endl;
      return;
    }

    Block::InnerData data = Block::blocks.at(_uid);
    std::chrono::time_point< std::chrono::system_clock > endClock;
    endClock = std::chrono::system_clock::now();
    std::chrono::duration< double > elapsed_seconds =
        endClock - data.startClock;

#if(defined(UNIX))
    std::cout << "\033[0m\033[1m\033[32m"; // Bold + green color
#endif
    std::cout << "################";
    for(unsigned int i = 0; i < data.message.length(); ++i)
    {
      std::cout << "#";
    }
    std::cout << "#" << std::endl;

    std::cout << "# End block:   " << data.message << " #" << std::endl;
    std::cout << "# Elapsed time: " << elapsed_seconds.count() << "s #"
              << std::endl;

    std::cout << "################";
    for(unsigned int i = 0; i < data.message.length(); ++i)
    {
      std::cout << "#";
    }
    std::cout << "#" << std::endl;
#if(defined(UNIX))
    std::cout << "\033[0m"; // reset
#endif

    Block::blocks.erase(_uid);
  }

protected:
  struct InnerData
  {
    InnerData() = default;
    ~InnerData() = default;
    std::string message;
    std::chrono::time_point< std::chrono::system_clock > startClock;
  };

  static std::map< std::string, Block::InnerData > blocks;
};

} // namespace FEVV

std::map< std::string, FEVV::Block::InnerData > FEVV::Block::blocks =
    std::map< std::string, FEVV::Block::InnerData >();
