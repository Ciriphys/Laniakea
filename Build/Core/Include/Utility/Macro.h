#pragma once

#ifdef LNK_BUILD_DLL
	#ifdef LNK_WIN
		#define LNK_API __declspec(dllexport)
	#else
		#define LNK_API
	#endif
#else
	#ifdef LNK_WIN
		#define LNK_API __declspec(dllimport)
	#else
		#define LNK_API
	#endif
#endif

#ifdef LNK_WIN
	constexpr auto LNK_PROJECT_PATH = "Core\\";
#else
	constexpr auto LNK_PROJECT_PATH = "Core/";

	#ifndef __FUNCSIG__
	#define __FUNCSIG__ __FUNCTION__
	#endif
#endif

#define EXTRACT_BYTE(from, n) ((from >> 8 * (n - 1)) & 0xff)
#define SQR(x) ((x) * (x))

#ifndef M_PI 
	#define M_PI (3.14159265358979)
#endif
