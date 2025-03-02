#include "lnkpch.h"

#include <Utility/Debug.h>

namespace Laniakea {

	void Assert(const char* expression, const char* file, const char* function, int line, const char* message)
	{
		printf("Assertion failed (%s) in file %s,\nin function %s, at line %d:\n%s", expression, file, function, line, message);
		abort();
	}

	void DisplayFunctionInfo(const char* functionName, const char* callerFunction)
	{
		printf("%s requested a debug information call for function \"%s\".\n", callerFunction, functionName);
	}

	std::string RelativeToBuildPath(std::string file)
	{
		std::string path;
		std::string prjPath = LNK_PROJECT_PATH;
		unsigned long long int index = file.find(prjPath) + prjPath.length();
		path = file.substr(index);

		return path;
	}

	LNK_API void DisplayLog(const char* msg)
	{
		std::cout << "[LOG]: ";
		std::cout << msg << std::endl;
	}
}