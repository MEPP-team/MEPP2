#pragma once

#if defined(mepp_RECURSES)
#error Recursive header files inclusion detected in mepp.h
#else // defined(mepp_RECURSES)
/** Prevents recursive inclusion of headers. */
#define mepp_RECURSES

#if !defined mepp_h
/** Prevents repeated inclusion of headers. */
#define mepp_h

// -----------------------------------------------------------------------
#define MEPP_VERSION "v0.7.2 - 15/11/2018"

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

#endif // !defined mepp_h

#undef mepp_RECURSES
#endif // else defined(mepp_RECURSES)
