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

class Result
{
public:
	Result(const mpz_t &x, const mpz_t &y, const mpz_t &p)
	{
		mpz_init(m_x);
		mpz_init(m_y);
		mpz_init(m_p);
		mpz_set(m_x, x);
		mpz_set(m_y, y);
		mpz_set(m_p, p);
	}
	Result(const Result &other)
	{
		mpz_init(m_x);
		mpz_init(m_y);
		mpz_init(m_p);
		mpz_set(m_x, other.m_x);
		mpz_set(m_y, other.m_y);
		mpz_set(m_p, other.m_p);
	}
	Result()
	{
		mpz_init(m_x);
		mpz_init(m_y);
		mpz_init(m_p);
	}
	virtual ~Result()
	{
		mpz_clear(m_x);
		mpz_clear(m_y);
		mpz_clear(m_p);
	}

	Result &operator = (const Result &other)
	{
		mpz_set(m_x, other.m_x);
		mpz_set(m_y, other.m_y);
		mpz_set(m_p, other.m_p);
	}
	
	std::string to_s()
	{
		std::string r;

		char *szx = mpz_get_str(0, 16, m_x);
		r += "0x";
		r += szx;
		r += "L";
		free(szx);

		r += " ";

		char *szy = mpz_get_str(0, 16, m_y);
		r += "0x";
		r += szy;
		r += "L";
		free(szy);

		r += " ";

		char *szp = mpz_get_str(0, 16, m_p);
		r += "0x";
		r += szp;
		r += "L";
		free(szp);

		return r;
	}

private:
	mpz_t m_x;
	mpz_t m_y;
	mpz_t m_p;
};

static std::mutex countMutex;
static std::mutex resultMutex;
static std::list<Result> resultList;
static unsigned int tryCount = 0;

void inverse_mod(mpz_t &out, const mpz_t &_a, const mpz_t &m)
{
	mpz_t a;
	mpz_init(a);
	mpz_set(a, _a);

	if (mpz_cmp_ui(a, 0) < 0 || mpz_cmp(m, a) <= 0)
	{
		mpz_mod(a, _a, m);
	}

	mpz_t c;
	mpz_init(c);
	mpz_set(c, a);
	mpz_t d;
	mpz_init(d);
	mpz_set(d, m);

	mpz_t uc;
	mpz_t vc;
	mpz_t ud;
	mpz_t vd;
	mpz_init(uc);
	mpz_set_ui(uc, 1);
	mpz_init(vc);
	mpz_set_ui(vc, 0);
	mpz_init(ud);
	mpz_set_ui(ud, 0);
	mpz_init(vd);
	mpz_set_ui(vd, 1);

	mpz_t q;
	mpz_init(q);
	mpz_t tmp;
	mpz_init(tmp);
	mpz_t tmp2;
	mpz_init(tmp2);
	mpz_t tmp3;
	mpz_init(tmp3);

	while(mpz_cmp_ui(c, 0) != 0)
	{
		mpz_set(tmp, c);
		mpz_div(q, d, c);
		mpz_mod(tmp2, d, c);
		mpz_set(c,tmp2);
		mpz_set(d, tmp);

		mpz_mul(tmp, q, uc);
		mpz_set(tmp2, uc);
		mpz_sub(uc, ud, tmp);

		mpz_mul(tmp, q, vc);
		mpz_set(tmp3, vc);
		mpz_sub(vc, vd, tmp);

		mpz_set(ud, tmp2);
		mpz_set(vd, tmp3);
	}

	//assert(mpz_cmp_ui(d, 1) == 0);

	if (mpz_cmp_ui(ud, 0) > 0)
		mpz_set(out, ud);
	else
		mpz_add(out, ud, m);
	
	mpz_clear(tmp3);
	mpz_clear(tmp2);
	mpz_clear(tmp);
	mpz_clear(q);
	mpz_clear(vd);
	mpz_clear(ud);
	mpz_clear(vc);
	mpz_clear(uc);
	mpz_clear(c);
	mpz_clear(d);
	mpz_clear(a);
}

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

class CurveFp
{
public:
	CurveFp(const mpz_t &p, const mpz_t &a, const mpz_t b)
	{
		mpz_init(m_p);
		mpz_init(m_a);
		mpz_init(m_b);
		mpz_set(m_p, p);
		mpz_set(m_a, a);
		mpz_set(m_b, b);
		m_inv = false;
	}
	CurveFp(const CurveFp &other)
	{
		mpz_init(m_p);
		mpz_init(m_a);
		mpz_init(m_b);
		mpz_set(m_p, other.m_p);
		mpz_set(m_a, other.m_a);
		mpz_set(m_b, other.m_b);
		m_inv = other.m_inv;
	}
	CurveFp()
	{
		mpz_init(m_p);
		mpz_init(m_a);
		mpz_init(m_b);
		m_inv = true;
	}
	virtual ~CurveFp()
	{
		mpz_clear(m_p);
		mpz_clear(m_a);
		mpz_clear(m_b);
	}

	CurveFp &operator = (const CurveFp &other)
	{
		mpz_set(m_p, other.m_p);
		mpz_set(m_a, other.m_a);
		mpz_set(m_b, other.m_b);
		m_inv = other.m_inv;
		return *this;
	}

	bool operator == (const CurveFp &other) const
	{
		if (m_inv && other.m_inv)
			return true;

		return (m_inv == other.m_inv &&
			    mpz_cmp(m_p, other.m_p) == 0 &&
			    mpz_cmp(m_a, other.m_a) == 0 &&
				mpz_cmp(m_b, other.m_b) == 0);
	}

	bool contains_point(const mpz_t &x, const mpz_t &y) const
	{
		mpz_t ysq;
		mpz_init(ysq);
		mpz_pow_ui(ysq, y, 2);

		mpz_t xt;
		mpz_init(xt);
		mpz_pow_ui(xt, x, 3);

		mpz_t ax;
		mpz_init(ax);
		mpz_mul(ax, m_a, x);

		mpz_t xtpa;
		mpz_init(xtpa);
		mpz_add(xtpa, xt, ax);

		mpz_t xtpapb;
		mpz_init(xtpapb);
		mpz_add(xtpapb, xtpa, m_b);

		mpz_t sub;
		mpz_init(sub);
		mpz_sub(sub, ysq, xtpapb);

		mpz_t mod;
		mpz_init(mod);
		mpz_mod(mod, sub, m_p);

		bool r = mpz_cmp_ui(mod, 0) == 0;

		mpz_clear(mod);
		mpz_clear(sub);
		mpz_clear(xtpapb);
		mpz_clear(xtpa);
		mpz_clear(ax);
		mpz_clear(xt);
		mpz_clear(ysq);

		return r;
	}

	const mpz_t &p() const { return m_p; }
	const mpz_t &a() const { return m_a; }
private:
	mpz_t m_p;
	mpz_t m_a;
	mpz_t m_b;
	bool m_inv;
};

class Point
{
public:
	Point(const CurveFp &curve, const mpz_t &x, const mpz_t &y, const mpz_t &o)
		:m_curve(curve)
	{
		mpz_init(m_x);
		mpz_init(m_y);
		mpz_init(m_o);
		mpz_set(m_x, x);
		mpz_set(m_y, y);
		mpz_set(m_o, o);
		m_inv = false;

		/*assert(curve.contains_point(x, y));
		if (mpz_cmp_ui(o, 0) != 0)
		{
			//assert(Point((*this) * o).m_inv); <-- infinite recursion!
		}*/
	}
	Point(const Point &other)
		:m_curve(other.m_curve)
	{
		mpz_init(m_x);
		mpz_init(m_y);
		mpz_init(m_o);
		mpz_set(m_x, other.m_x);
		mpz_set(m_y, other.m_y);
		mpz_set(m_o, other.m_o);
		m_inv = other.m_inv;
	}
	Point()
	{
		mpz_init(m_x);
		mpz_init(m_y);
		mpz_init(m_o);
		m_inv = true;
	}
	virtual ~Point()
	{
		mpz_clear(m_x);
		mpz_clear(m_y);
		mpz_clear(m_o);
	}

	Point &operator = (const Point &other)
	{
		m_curve = other.m_curve;
		mpz_set(m_x, other.m_x);
		mpz_set(m_y, other.m_y);
		mpz_set(m_o, other.m_o);
		m_inv = other.m_inv;
		return *this;
	}
	Point operator + (const Point &other) const
	{
		if (other.m_inv)
			return *this;
		if (m_inv)
			return other;

		//assert(m_curve == other.m_curve);

		if (mpz_cmp(m_x, other.m_x) == 0)
		{
			mpz_t yadd;
			mpz_t ymod;
			mpz_init(yadd);
			mpz_init(ymod);

			mpz_add(yadd, m_y, other.m_y);
			mpz_mod(ymod, yadd, m_curve.p());

			int cmp = mpz_cmp_ui(ymod, 0);

			mpz_clear(yadd);
			mpz_clear(ymod);

			if (cmp == 0)
				return Point();
			else
				return dbl();
		}

		mpz_t ysub;
		mpz_init(ysub);
		mpz_sub(ysub, other.m_y, m_y);

		mpz_t xsub;
		mpz_init(xsub);
		mpz_sub(xsub, other.m_x, m_x);

		mpz_t imod;
		mpz_init(imod);
		inverse_mod(imod, xsub, m_curve.p());

		mpz_t mul;
		mpz_init(mul);
		mpz_mul(mul, ysub, imod);

		mpz_t l;
		mpz_init(l);
		mpz_mod(l, mul, m_curve.p());

		mpz_t lsq;
		mpz_init(lsq);
		mpz_pow_ui(lsq, l, 2);

		mpz_t lsq2;
		mpz_init(lsq2);
		mpz_sub(lsq2, lsq, m_x);

		mpz_t lsq3;
		mpz_init(lsq3);
		mpz_sub(lsq3, lsq2, other.m_x);

		mpz_t x3;
		mpz_init(x3);
		mpz_mod(x3, lsq3, m_curve.p());

		mpz_t xmx3;
		mpz_init(xmx3);
		mpz_sub(xmx3, m_x, x3);

		mpz_t ymul;
		mpz_init(ymul);
		mpz_mul(ymul, l, xmx3);

		mpz_t ys;
		mpz_init(ys);
		mpz_sub(ys, ymul, m_y);

		mpz_t y3;
		mpz_init(y3);
		mpz_mod(y3, ys, m_curve.p());

		mpz_t mpZero;
		mpz_init(mpZero);
		mpz_set_ui(mpZero, 0);

		Point r(m_curve, x3, y3, mpZero);

		mpz_clear(mpZero);
		mpz_clear(y3);
		mpz_clear(ys);
		mpz_clear(ymul);
		mpz_clear(xmx3);
		mpz_clear(x3);
		mpz_clear(lsq3);
		mpz_clear(lsq2);
		mpz_clear(lsq);
		mpz_clear(l);
		mpz_clear(mul);
		mpz_clear(imod);
		mpz_clear(xsub);
		mpz_clear(ysub);
		return r;
	}
	Point operator * (const mpz_t &mul)
	{
		if (m_inv)
			return Point();

		mpz_t e;
		mpz_init(e);
		mpz_set(e, mul);

		if (mpz_cmp_ui(m_o, 0) != 0)
			mpz_mod(e, mul, m_o);

		if (mpz_cmp_ui(e, 0) == 0)
			return Point();

		mpz_t e3;
		mpz_init(e3);
		mpz_mul_ui(e3, e, 3);

		mpz_t yneg;
		mpz_init(yneg);
		mpz_neg(yneg, m_y);

		Point negative_self(m_curve, m_x, yneg, m_o);

		mpz_t lmb;
		mpz_init(lmb);
		leftmost_bit(lmb, e3);

		mpz_t i;
		mpz_init(i);
		mpz_div_ui(i, lmb, 2);

		Point r(*this);

		mpz_t i2;
		mpz_t e3i;
		mpz_t ei;
		mpz_init(i2);
		mpz_init(e3i);
		mpz_init(ei);

		while(mpz_cmp_ui(i, 1) > 0)
		{
			r = r.dbl();

			mpz_and(e3i, e3, i);
			mpz_and(ei, e, i);

			if (mpz_cmp_ui(e3i, 0) != 0 && mpz_cmp_ui(ei, 0) == 0)
				r = r + *this;
			if (mpz_cmp_ui(e3i, 0) == 0 && mpz_cmp_ui(ei, 0) != 0)
				r = r + negative_self;
			mpz_div_ui(i2, i, 2);
			mpz_set(i, i2);
		}

		mpz_clear(ei);
		mpz_clear(e3i);
		mpz_clear(i2);
		mpz_clear(i);
		mpz_clear(lmb);
		mpz_clear(yneg);
		mpz_clear(e3);
		mpz_clear(e);

		return r;
	}

	const mpz_t &x() const { return m_x; }
	const mpz_t &y() const { return m_y; }

	Point dbl() const
	{
		if (m_inv)
			return Point();

		mpz_t xsq;
		mpz_init(xsq);
		mpz_pow_ui(xsq, m_x, 2);

		mpz_t xsq3;
		mpz_init(xsq3);
		mpz_mul_ui(xsq3, xsq, 3);

		mpz_t xsq3a;
		mpz_init(xsq3a);
		mpz_add(xsq3a, xsq3, m_curve.a());

		mpz_t ydbl;
		mpz_init(ydbl);
		mpz_mul_ui(ydbl, m_y, 2);

		mpz_t imod;
		mpz_init(imod);
		inverse_mod(imod, ydbl, m_curve.p());

		mpz_t mul;
		mpz_init(mul);
		mpz_mul(mul, xsq3a, imod);

		mpz_t l;
		mpz_init(l);
		mpz_mod(l, mul, m_curve.p());

		mpz_t lsq;
		mpz_init(lsq);
		mpz_pow_ui(lsq, l, 2);

		mpz_t x2;
		mpz_init(x2);
		mpz_mul_ui(x2, m_x, 2);

		mpz_t lsqmx2;
		mpz_init(lsqmx2);
		mpz_sub(lsqmx2, lsq, x2);

		mpz_t x3;
		mpz_init(x3);
		mpz_mod(x3, lsqmx2, m_curve.p());

		mpz_t xmx3;
		mpz_init(xmx3);
		mpz_sub(xmx3, m_x, x3);

		mpz_t ymul;
		mpz_init(ymul);
		mpz_mul(ymul, l, xmx3);

		mpz_t ys;
		mpz_init(ys);
		mpz_sub(ys, ymul, m_y);

		mpz_t y3;
		mpz_init(y3);
		mpz_mod(y3, ys, m_curve.p());

		mpz_t mpZero;
		mpz_init(mpZero);
		mpz_set_ui(mpZero, 0);

		Point r(m_curve, x3, y3, mpZero);
		
		mpz_clear(mpZero);
		mpz_clear(y3);
		mpz_clear(ys);
		mpz_clear(ymul);
		mpz_clear(xmx3);
		mpz_clear(x3);
		mpz_clear(lsqmx2);
		mpz_clear(x2);
		mpz_clear(lsq);
		mpz_clear(l);
		mpz_clear(mul);
		mpz_clear(imod);
		mpz_clear(ydbl);
		mpz_clear(xsq3a);
		mpz_clear(xsq3);
		mpz_clear(xsq);

		return r;
	}
private:
	mpz_t m_x;
	mpz_t m_y;
	mpz_t m_o;
	CurveFp m_curve;
	bool m_inv;
};

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
		//mpz_set_ui(priv, 99337);
		unsigned int shuffle = rand() % 3;
		switch(shuffle)
		{
		case 0:
			myp = myp + add1;
			mpz_add(temppriv, temppriv, additor1);
			break;
		case 2:
			myp = myp + add2;
			mpz_add(temppriv, temppriv, additor2);
			break;
		default:
			myp = myp + add3;
			mpz_add(temppriv, temppriv, additor3);
			break;
		}

		for(unsigned int i = 0; i < tries; ++ i)
		{
			myp = myp + g;
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
		return -1;

	if (!socket.Send((const uint8 *)outbuf.data(), outbuf.size()))
		return -1;

	char resp[1024] = {};
	socket.Receive(1024);
	memcpy(resp, socket.GetData(), socket.GetBytesReceived());
	int r = 0;
	if (strstr(resp, "submission is great"))
		r = 1;
	socket.Close();

	return r;
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

	puts("C++ BTC address space analyzer by gh2k");
	puts("Tips: 1gh2k13JzLTxDz2vmo8RH7HcjTLLg5kdc ;-)");
	puts("");

	printf("Starting %d threads. Hang on to your hats...\n", threadcount);
	printf("Sending payouts to %s\n", btc.data());

	std::list<std::thread *> threads;
	for( int i = 0; i < 12; ++ i )
	{
		threads.push_back(new std::thread(work));
	}
	
    std::chrono::milliseconds sleeptime( 100 );

	time_t secs = time(0);
	time_t starttime = time(0);

	size_t resultCount = 0;

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

			resultMutex.lock();
			std::list<Result> newResults = resultList;
			resultList.clear();
			resultMutex.unlock();

			while (newResults.size())
			{
				std::string rs(newResults.front().to_s());
				switch (submit_share(btc, rs))
				{
				case -1:
					puts("Share send failed. Couldn't connect.");
					break;
				case 1:
					puts("Share accepted!");
					++ resultCount;
					break;
				case 0:
					puts("Share REJECTED.");
					break;
				}

				newResults.pop_front();
			}

			double resultsPerHour = static_cast<double>(resultCount) / (static_cast<double>(curtime - starttime) / 60.0 / 60.0);

			printf("%d keys/sec - %d found (%.2f / hour)\n", curtries / 10, resultCount, resultsPerHour);
		}
	}
}
