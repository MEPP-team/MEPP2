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
#include <boost/assert.hpp>

namespace FEVV {

/**
 * class Assert
 * \brief This class may be used to draw Asserts in your code.
 *
 * The code below shows a possible use of this class.
 * @code
 * #include "Base/Assert.h"
 *
 * using namespace FEVV;
 * // ...
 * {
 *   int a = 7;
 *
 *   Assert::check( a <= 10, "my value is greater than 10", "main()"); // will
 * not raised an assert Assert::check( a >= 10, "my value is lower than 10",
 * "main()"); // will raised an assert
 * }
 * @endcode
 *
 */
class Assert
{
public:
  /**
   * Default constructor.
   * @note Disable by default.
   */
  Assert() = delete;

  /**
   * Copy constructor.
   * @note Disable by default.
   */
  Assert(const Assert &_Assert) = delete;

  /**
   * Destructor.
   * @note Disable by default.
   */
  ~Assert() = delete;

  /**
   * Check if we need to raise an assert.
   *
   * @param[in]    _check      if _check is false, assert is raised.
   * @param[in]    _message    a message to show in standard output if an
   *                           assert is raised.
   * @param[in]    _where      the function name where the assert is checked.
   *                           (default: "").
   *
   * @return       _check
   */
  static bool check(const bool _check,
                    const std::string &_message,
                    const std::string &_where = "")
  {
    std::string m;
    m = (_where != "") ? "[" + _where + "] " : "";
    m += _message;

    if(!_check)
    {
#if(defined(UNIX))
      std::cout << "\033[0m\033[1m\033[31m"; // Bold + red color
#endif
      std::cout << m;
#if(defined(UNIX))
      std::cout << "\033[0m"; // reset
#endif
      std::cout << std::endl;
    }

    // BOOST_ASSERT_MSG( _check, m.c_str() );

    return _check;
  }
};
} // namespace FEVV
