#include <stdio.h>
#include <time.h>
#include <set>
#include <thread>
#include <list>
#include <mutex>
#include <iostream>


#ifdef USE_GMP
#include <gmp.h>
#else //USE_GMP
#include "mpir.h"
#endif //USE_GMP
#include "targets.h"

#include "ActiveSocket.h"

//#include "vld.h"

static const char *VERSION = "0.8-opencl";

void leftmost_bit(mpz_t &out, const mpz_t &x)
{
	//assert(mpz_cmp_ui(x, 0) > 0);
	mpz_t r;
	mpz_t r2;
	mpz_init(r);
	mpz_set_ui(r, 1);
	mpz_init(r2);

	while(mpz_cmp(r, x) <= 0)
	{
		mpz_mul_ui(r2, r, 2);
		mpz_set(r, r2);
	}

	mpz_div_ui(out, r, 2);

	mpz_clear(r2);
	mpz_clear(r);
}

// horrible raw cpp includes! I'm so lazy. deal with it.
#include "OpenCLFramework.h"
#include "result.h"
#include "curvefp.h"
#include "point.h"

static std::mutex countMutex;
static std::mutex resultMutex;
static std::mutex statsMutex;
static std::list<Result> resultList;
static unsigned int tryCount = 0;
static unsigned int resultCount = 0;
static unsigned int duplicateCount = 0;
static unsigned int rejectCount = 0;
static unsigned int queuedCount = 0;

void work()
{
	/**
	Evil-Knievel: Before starting OpenCL we must get the starting point, these methods are 1:1 from gh2k's port!	
	**/



	static const unsigned int tries = 100;

	uint64 seed64 = (uint64)time(0) << 32 | static_cast<uint32>(std::hash<std::thread::id>()(std::this_thread::get_id()));
	uint32 seed32 = static_cast<uint32>(time(0)) ^ static_cast<uint32>(std::hash<std::thread::id>()(std::this_thread::get_id()));
	srand(seed32);

	std::set<uint32> targets;
	initTargets(targets);

	mpz_t seed;
	mpz_init(seed);
	mpz_set_ui(seed, seed64);

	gmp_randstate_t randstate;
	gmp_randinit_mt(randstate);
	gmp_randseed(randstate, seed);

	mpz_t p;
	mpz_t b;
	mpz_t a;
	mpz_init(a);
	mpz_set_ui(a, 0);
	mpz_init(b);
	mpz_set_ui(b, 7);
	mpz_init(p);
	mpz_set_str(p, "FFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFEFFFFFC2F", 16);

	mpz_t gx;
	mpz_t gy;
	mpz_t r;
	mpz_init(gx);
	mpz_init(gy);
	mpz_init(r);
	mpz_set_str(gx, "79BE667EF9DCBBAC55A06295CE870B07029BFCDB2DCE28D959F2815B16F81798", 16);
	mpz_set_str(gy, "483ada7726a3c4655da4fbfc0e1108a8fd17b448a68554199c47d08ffb10d4b8", 16);
	mpz_set_str(r, "FFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFEBAAEDCE6AF48A03BBFD25E8CD0364141", 16);

	CurveFp curve(p, a, b);
	Point point(curve, gx, gy, r);
	Point g(point);

	mpz_t priv;
	mpz_init(priv);
	mpz_urandomb(priv, randstate, 256);


	mpz_t temppriv;
	mpz_init_set(temppriv, priv);

	Point myp = g * priv;

	

	/**
	Evil-Knievel: Create all OpenCL Related Entities here. We should think about moving it to OpenCLFramework.cpp
	but I am a C++ noob so I make it quick and dirty.
	**/

	// Create Function Instance
	OpenCLFramework clFunctions;

	// open context
	cl::Context context(CL_DEVICE_TYPE_GPU);
	std::vector<cl::Device> devices= context.getInfo<CL_CONTEXT_DEVICES>();
	std::cout << "(opencl) will use the following device: " << devices[0].getInfo<CL_DEVICE_NAME>()  << std::endl;
	
	//Make a queue to put jobs on the first compute device
	cl::CommandQueue queue(context, devices[0], CL_QUEUE_PROFILING_ENABLE);


	// create kernel from cl file
	cl::Program program = clFunctions.buildProgramFromSource(&context, "kernel/evilknievel.cl", "");

	// Host Memory Stuff that will be mapped to a cl::Buffer somewhere
	bignum points_in[2]; // <------ here goes the starting point, e.g. random myp point
	bignum points_found[2];
	cl_uint found_cell_index = 0;
	cl_uint found = 0;

	// create buffers for CL Device
	cl::Buffer found_buffer(context, CL_MEM_WRITE_ONLY, sizeof(cl_uint)* 1);
	cl::Buffer found_cell(context, CL_MEM_READ_WRITE, sizeof(cl_uint)* 1);
	cl::Buffer points_in_buffer(context, CL_MEM_READ_WRITE  /* | CL_MEM_COPY_HOST_PTR */, sizeof(bignum)* 2);
	cl::Buffer found_point(context, CL_MEM_WRITE_ONLY  /* | CL_MEM_COPY_HOST_PTR */, sizeof(bignum)* 2);

	// Write necessary buffers
	queue.enqueueWriteBuffer(found_cell, CL_TRUE, 0, sizeof(cl_uint)* 1, &found_cell_index);
	queue.enqueueWriteBuffer(points_in_buffer, CL_TRUE, 0, sizeof(bignum)* 2, points_in);

	// Prepare Opencl Kernel for execution
	cl::Kernel kernel = cl::Kernel(program, "search_triplets");
	kernel.setArg(0, found_buffer);
	kernel.setArg(1, points_in_buffer);
	kernel.setArg(2, found_point);
	kernel.setArg(3, found_cell);
	unsigned int i_global_work_size = 81920; // <------- these are hardcoded, maybe we need some autodetection like cgminer does !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	unsigned int i_local_work_size = 256;  // <------- these are hardcoded, maybe we need some autodetection like cgminer does !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	/* it is ugly to create this mpz_number here, but we have to do this after setting global work size*/
	mpz_t work_size_incrementor;
	mpz_init(work_size_incrementor);
	mpz_set_ui(work_size_incrementor, i_global_work_size);
	// init a fixed point to be added each iteration
	Point g_worksize_increment = g*work_size_incrementor;
	const cl::NDRange global_work_offset(0, 0);
	const cl::NDRange global_work_size(i_global_work_size, 1, 1);
	const cl::NDRange local_work_size(i_local_work_size, 1, 1);




	while( true )
	{
		
		// Launch the kernel and do 81920 parallel calculations
		int ret = queue.enqueueNDRangeKernel(kernel, global_work_offset, global_work_size, local_work_size, NULL, NULL);
		if (ret != CL_SUCCESS){
			std::cerr << "clEnqueueNDRangeKernel fail. Error " << ret << std::endl;
			if (ret == CL_OUT_OF_RESOURCES){
				std::cerr << "Kernel out of resources!!! Try lowering global/local worksize in the source code! You GPU seems to be 'low end'" << std::endl;
			}
		}

		// Now check if the kernel reported a found point back
		queue.enqueueReadBuffer(found_buffer, CL_TRUE, 0, sizeof(cl_int)* 1, &found);
		queue.enqueueReadBuffer(points_in_buffer, CL_TRUE, 0, sizeof(bignum)* 2, points_in);

		if (found > 0){
			// Yes, the GPU actually did find a value triplet.
			
			// Now let us read the result from the GPU!!!
			// WARNING: AT THE MOMENT THE GPU CAN ONLY FIND ONE TRIPLET PER KERNEL/RUN and we have only one found and one points_found variable shared with the kernel
			queue.enqueueReadBuffer(found_cell, CL_TRUE, 0, sizeof(cl_uint)* 1, &found_cell_index); // <------ found index is the thread number which actually found the triplet
			// (RATHER THAN GETTING THIS INFO FROM THE GPU, I TRY TO CALC THE X,Y A FEW LINES BELOW THIS ON MY OWN! NOT SURE IF ITS THE BEST IDEA ....... queue.enqueueReadBuffer(found_point, CL_TRUE, 0, sizeof(bignum)* 2, points_found);  // <----- and here the found x and y are given back

			found = 0; // reset found
			queue.enqueueWriteBuffer(found_buffer, CL_TRUE, 0, sizeof(cl_int)* 1, &found); // also reset found variable on GPU memory

			/**
			At this point we can reconstruct the triplet from points_found[0], points_found[1] and temppriv!!!
			in face we dont take points_found for now but calculate the found point from myp and temppriv and found_cell_index
			**/

			// IMPORTANT!! THIS HERE IS VERY INEFFICIENT, MAYBE YOU GUYS HAVE A BETTER IDEA!

			Point myp_found = myp;
			mpz_t found_cell_index_mpz;
			mpz_init(found_cell_index_mpz);
			mpz_set_ui(found_cell_index_mpz, found_cell_index);

			Point adding_found_index_offset = g*found_cell_index_mpz;
			myp_found.add(adding_found_index_offset);
			mpz_add(temppriv, temppriv, found_cell_index_mpz);

			resultMutex.lock();
			resultList.push_back(Result(myp_found.x(), myp_found.y(), temppriv));
			resultMutex.unlock();

			mpz_add(temppriv, temppriv, found_cell_index_mpz); // revert temppriv as it is going to be incremented at the end of the loop anyway
			mpz_clear(found_cell_index_mpz);

		}

		// Before doing anything else, increment temppriv by the number of iterations we have done so far!
		mpz_add(temppriv, temppriv, work_size_incrementor);
		// and also update our starting point here so we keep track how far we moved already
		myp.add(g_worksize_increment);

		countMutex.lock();
		tryCount += i_global_work_size;	// We have actually did so many tries as the global_work_size says as we are having this many parallel threads on the GPU
		countMutex.unlock();



		/**
		resultMutex.lock();
		resultList.push_back(Result(myp.x(), myp.y(), temppriv));
		resultMutex.unlock();
			
		countMutex.lock();
		tryCount += tries;
		countMutex.unlock();**/
	}

	mpz_clear(work_size_incrementor);
	mpz_clear(gx);
	mpz_clear(gy);
	mpz_clear(r);
	mpz_clear(p);
	mpz_clear(b);
	mpz_clear(a);
	mpz_clear(priv);
	mpz_clear(temppriv);
	mpz_clear(seed);
}

int submit_share(const std::string &btc, const std::string &result)
{
	std::string outbuf;
	outbuf += "result ";
	outbuf += btc;
	outbuf += " ";
	outbuf += result;
	outbuf += "\n";

	CActiveSocket socket;
	socket.Initialize();

	if (!socket.Open((const uint8 *)"stargate.bitwarrant.com", 4444))
	//if (!socket.Open((const uint8 *)"localhost", 4444))
		return -1;

	if (!socket.Send((const uint8 *)outbuf.data(), outbuf.size()))
		return -1;

	int r = 0;

	char resp[1024] = {};
	if (socket.Receive(1024) > 0)
	{
		memcpy(resp, socket.GetData(), socket.GetBytesReceived());
		
		if (strstr(resp, "submission is great"))
			r = 1;
		if (strstr(resp, "submission is already here"))
			r = 2;
	}
	socket.Close();

	return r;
}

void submitter(const std::string btc)
{
	std::list<Result> queuedResults;
	
    std::chrono::milliseconds sleeptime( 100 );

	while( true )
	{
		resultMutex.lock();
		queuedResults.splice(queuedResults.end(), resultList);
		resultMutex.unlock();
		
		std::list<Result> heldResults;

		if (queuedResults.size())
		{
			Result nextResult(queuedResults.front());
			queuedResults.pop_front();

			std::string rs(nextResult.to_s());
			int r = submit_share(btc, rs);

			statsMutex.lock();
			switch (r)
			{
			case -1:
				heldResults.push_back(nextResult);
				break;
			case 1:
				++ resultCount;
				break;
			case 0:
				++ rejectCount;
				break;
			case 2:
				++ duplicateCount;
				break;
			}
			queuedCount = static_cast<unsigned int>(queuedResults.size() + heldResults.size());
			statsMutex.unlock();

			queuedResults.splice(queuedResults.end(), heldResults);
		} else {
			std::this_thread::sleep_for( sleeptime );
		}
	}
}

int main(int argc, char **argv)
{
	if ( argc != 2 )
	{
		printf("Usage: %s bitcoinaddress\n", argv[0]);
		return -1;
	}

	std::string btc(argv[1]);

	

	std::list<Result> held_results;

	printf("C++ BTC address space analyzer by gh2k (modified to support OpenCL), v%s\n", VERSION);
	puts("Tips: 1gh2k13JzLTxDz2vmo8RH7HcjTLLg5kdc ;-)");
	puts("");

	printf("Starting OpenCL kernel. Hang on to your hats...\n");
	printf("Sending payouts to %s\n", btc.data());

	std::list<std::thread *> threads;
	// EVILKNIEVEL: For threads only one worker thread, as the parallelization will be done on the GPU itself!
	threads.push_back(new std::thread(work));
	
	threads.push_back(new std::thread(submitter, btc));
	
    std::chrono::milliseconds sleeptime( 100 );

	time_t secs = time(0);
	time_t starttime = time(0);

	while( true )
	{
		std::this_thread::sleep_for( sleeptime );
		time_t curtime = time(0);
		if (curtime > secs + 10)
		{
			secs = curtime;
			
			countMutex.lock();
			unsigned int curtries = tryCount;
			tryCount = 0;
			countMutex.unlock();

			statsMutex.lock();
			double resultsPerHour = static_cast<double>(resultCount) / (static_cast<double>(curtime - starttime) / 60.0 / 60.0);

			printf("%d keys/sec - %d accepted (%.2f / hour) - %d rej, %d dup, %d que\n", curtries / 10, resultCount, resultsPerHour, rejectCount, duplicateCount, queuedCount);
			statsMutex.unlock();
		}
	}
}
