#include <stdio.h>
#include <time.h>
#include <set>
#include <thread>
#include <list>
#include <mutex>

#ifdef USE_GMP
#include <gmp.h>
#else //USE_GMP
#include "mpir.h"
#endif //USE_GMP
#include "targets.h"

#include "ActiveSocket.h"

//#include "vld.h"

static const char *VERSION = "0.7";

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
	mpz_t additor1;
	mpz_init(additor1);
	mpz_t additor2;
	mpz_init(additor2);
	mpz_t additor3;
	mpz_init(additor3);
	
	mpz_urandomb(priv, randstate, 256);
	mpz_urandomb(additor1, randstate, 230);
	mpz_urandomb(additor2, randstate, 230);
	mpz_urandomb(additor3, randstate, 128);

	mpz_t temppriv;
	mpz_init_set(temppriv, priv);
	
	Point myp = g * priv;

	Point add1 = g * additor1;
	Point add2 = g * additor2;
	Point add3 = g * additor3;

	while( true )
	{
		unsigned int shuffle = rand() % 3;
		switch(shuffle)
		{
		case 0:
			myp.add(add1);
			mpz_add(temppriv, temppriv, additor1);
			break;
		case 2:
			myp.add(add2);
			mpz_add(temppriv, temppriv, additor2);
			break;
		default:
			myp.add(add3);
			mpz_add(temppriv, temppriv, additor3);
			break;
		}

		mpz_mod(temppriv, temppriv, r);

		for(unsigned int i = 0; i < tries; ++ i)
		{
			myp.add(g);
			mpz_add_ui(temppriv, temppriv, 1);

			uint32 match = static_cast<uint32>(mpz_get_ui(myp.x()));

			if (targets.count(match))
			{
				resultMutex.lock();
				resultList.push_back(Result(myp.x(), myp.y(), temppriv));
				resultMutex.unlock();
			}
		}
		
		countMutex.lock();
		tryCount += tries;
		countMutex.unlock();
	}

	mpz_clear(gx);
	mpz_clear(gy);
	mpz_clear(r);
	mpz_clear(p);
	mpz_clear(b);
	mpz_clear(a);
	mpz_clear(priv);
	mpz_clear(temppriv);
	mpz_clear(additor1);
	mpz_clear(additor2);
	mpz_clear(additor3);
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
	if ( argc != 3 )
	{
		printf("Usage: %s bitcoinaddress num_threads\n", argv[0]);
		return -1;
	}

	std::string btc(argv[1]);

	int threadcount = atoi(argv[2]);
	if (threadcount < 1 || threadcount > 512)
	{
		puts("Invalid thread count");
		return -1;
	}

	std::list<Result> held_results;

	printf("C++ BTC address space analyzer by gh2k, v%s\n", VERSION);
	puts("Tips: 1gh2k13JzLTxDz2vmo8RH7HcjTLLg5kdc ;-)");
	puts("");

	printf("Starting %d threads. Hang on to your hats...\n", threadcount);
	printf("Sending payouts to %s\n", btc.data());

	std::list<std::thread *> threads;
	for( int i = 0; i < threadcount; ++ i )
	{
		threads.push_back(new std::thread(work));
	}
	
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
