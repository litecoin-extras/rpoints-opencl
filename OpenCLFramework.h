#pragma once
/**
OpenCL Basic Functions by Evil-Knievel
**/
//#define __CL_ENABLE_EXCEPTIONS
#include <CL/cl.hpp>

typedef cl_uint bn_word;
#define BN_NBITS 256
#define BN_WSHIFT 5
#define BN_WBITS (1 << BN_WSHIFT)
#define BN_NWORDS ((BN_NBITS/8) / sizeof(bn_word))
#define BN_WORDMAX 0xffffffff

#define MODULUS_BYTES \
	0xfffffc2f, 0xfffffffe, 0xffffffff, 0xffffffff, \
	0xffffffff, 0xffffffff, 0xffffffff, 0xffffffff

typedef struct {
	bn_word d[BN_NWORDS];
} bignum;

class OpenCLFramework
{
public:
	OpenCLFramework();
	std::string intToHex(unsigned int p);
	cl::Program buildProgramFromSource(cl::Context* context, std::string filename, std::string buildOptions);

	~OpenCLFramework();
};

