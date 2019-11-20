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

// -----------------------------------------------------------------------
#define MEPP_VERSION "v0.12.0 - 20/11/2019"

#define MAINWINDOW_TITLE QObject::tr("MEPP2 : 3D MEsh Processing Platform")

#if _WIN64 || __amd64__
#define PORTABLE_64_BIT
#define ARCHITECTURE QObject::tr("64 bits")
#else
#define PORTABLE_32_BIT
#define ARCHITECTURE QObject::tr("32 bits")
#endif

#ifndef NDEBUG
#define BUILD_TYPE QObject::tr("DEBUG")
#else
#define BUILD_TYPE QObject::tr("RELEASE")
#endif

#define MEPP_READER QObject::tr("USE MEPP GENERIC READER")

//#ifndef CGAL_VERSION_STR
//#define CGAL_xstr(s) #s
//#define CGAL_str(s) CGAL_xstr(s)
//#define CGAL_VERSION_STR CGAL_str(CGAL_VERSION)
//#endif

#define ORGANIZATION QObject::tr("LIRIS")
#define APPLICATION QObject::tr("MEPP2")
// -----------------------------------------------------------------------
