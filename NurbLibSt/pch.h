// pch.h: This is a precompiled header file.
// Files listed below are compiled only once, improving build performance for future builds.
// This also affects IntelliSense performance, including code completion and many code browsing features.
// However, files listed here are ALL re-compiled if any one of them is updated between builds.
// Do not add files here that you will be updating frequently as this negates the performance advantage.

#ifndef PCH_H
#define PCH_H

#ifndef VC_EXTRALEAN
#define VC_EXTRALEAN		// Exclude rarely-used stuff from Windows headers
#endif

#ifndef WINVER				
#define WINVER _WIN32_WINNT_WS08		
#endif


#define _CRT_SECURE_CPP_OVERLOAD_STANDARD_NAMES 1

// add headers that you want to pre-compile here
#include <SDKDDKVer.h>
#include <assert.h>
#include <math.h>
#include <tchar.h>
#include "framework.h"
#include "NurbsLibSt.h"

#endif //PCH_H
