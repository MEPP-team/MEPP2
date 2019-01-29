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
#include "FEVV/Tools/Math/MatrixOperations.hpp"

// must be the last include to avoid side effect of disabling NDEBUG
#undef NDEBUG
#include <cassert>


//------------------------------------------------------------------------------

int
main(int argc, const char **argv)
{
  //---------------------- Tests ---------------------

  std::cout << "testing  maximum(vector)" << std::endl;
  assert((FEVV::Math::Vector::Stats::maximum(
              std::vector< int >{5, 8, 4, 3, -1}) == 8));

  std::cout << "testing  minimum(vector)" << std::endl;
  assert((FEVV::Math::Vector::Stats::minimum(
              std::vector< int >{5, 8, 4, 3, -1}) == -1));

  std::cout << "testing  mean(vector)" << std::endl;
  assert((FEVV::Math::Vector::Stats::mean(std::vector< int >{5, 8, 4, 3, -1}) ==
          3));

  std::cout << "testing  median(vector)" << std::endl;
  assert((FEVV::Math::Vector::Stats::median(
              std::vector< int >{5, 8, 4, 3, -1}) == 4));

  std::cout << "testing  dot_product(vector1, vector2)" << std::endl;
  assert((FEVV::Math::Vector::dot_product(
              std::vector< int >{5, 8, 4, 3, -1},
              std::vector< int >{10, -2, -3, 2, 1}) == 27));

  std::cout << "testing  cross_product(vector1, vector2)" << std::endl;
  assert((FEVV::Math::Vector::cross_product(std::vector< int >{4, 3, -1},
                                            std::vector< int >{10, -2, -3}) ==
          std::vector< int >{-11, 2, -38}));

  //--------------------------------------------------

  std::cout << "Test passed.\n";
  return 0; // all tests passed
}
