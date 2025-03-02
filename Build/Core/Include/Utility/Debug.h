#pragma once

#include "lnkpch.h"

#include <Utility/Macro.h>

#ifdef LNK_DEBUG 
	#define LNK_ASSERT(expression, message) if(!(expression)) Laniakea::Assert(#expression, Laniakea::RelativeToBuildPath(__FILE__).c_str(), __FUNCTION__, __LINE__, message)
	#define LNK_DEBUG_INFO(function) DisplayFunctionInfo(#function, __FUNCTION__); function
#else 
	#define LNK_ASSERT(expression, message)
	#define LNk_DEBUG_INFO(function) function;
#endif

namespace Laniakea
{
	LNK_API void Assert(const char* expression, const char* file, const char* function, int line, const char* message);
	LNK_API void DisplayFunctionInfo(const char* functionName, const char* callerFunction);
	LNK_API std::string RelativeToBuildPath(std::string file);

	LNK_API void DisplayLog(const char* msg);
}