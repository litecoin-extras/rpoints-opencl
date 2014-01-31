#define __CL_ENABLE_EXCEPTIONS
#include "OpenCLFramework.h"
#include <iostream>
#include <sstream>
#include <fstream>

/**
OpenCL Basic Functions by Evil-Knievel
**/


using namespace std;
OpenCLFramework::OpenCLFramework()
{
}


OpenCLFramework::~OpenCLFramework()
{
}

std::string  OpenCLFramework::intToHex(unsigned int p){
	char str[64];
	sprintf(str, "%x", p);
	return std::string(str);
}

cl::Program& OpenCLFramework::buildProgramFromSource(cl::Context* context, std::string filename, std::string buildOptions) {
	// Read source file
	std::ifstream sourceFile(filename.c_str());
	if (sourceFile.fail())
		throw cl::Error(1, "Failed to open OpenCL source file");
	std::string sourceCode(
		std::istreambuf_iterator<char>(sourceFile),
		(std::istreambuf_iterator<char>()));

	std::cout << "(Kernel) source code size in bytes: " << sourceCode.length() << std::endl;
	std::cout << "(Kernel) starting to build for GPU device." << std::endl;
	cl::Program::Sources source(1, std::make_pair(sourceCode.c_str(), strlen(sourceCode.c_str())));
	std::cout << "(Kernel) cl source object has size: " << source.size() << std::endl;


	// Make program of the source code in the context
	cl::Program program = cl::Program(*context, source);

	VECTOR_CLASS<cl::Device> devices = context->getInfo<CL_CONTEXT_DEVICES>();

	// Build program for these specific devices
	try{
		program.build(devices, buildOptions.c_str());
	}
	catch (cl::Error error) {
		//if(error.err() == CL_BUILD_PROGRAM_FAILURE) {
		std::cout << "Build log:" << std::endl << program.getBuildInfo<CL_PROGRAM_BUILD_LOG>(devices[0]) << std::endl;
		//}
		throw error;
	}

	return program;

}