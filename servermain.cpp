#include <stdio.h>
#include <time.h>
#include <set>
#include <thread>
#include <list>
#include <vector>
#include <mutex>
#include <istream>
#include <iterator>
#include <algorithm>
#include <sstream>

#ifdef USE_GMP
#include <gmp.h>
#else //USE_GMP
#include "mpir.h"
#endif //USE_GMP
#include "targets.h"

#include "ActiveSocket.h"
#include "PassiveSocket.h"

#include "sqlite3.h"

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

#include "curvefp.h"
#include "point.h"

#define MAX_PACKET 1024

static std::string earning;
static std::mutex earningMutex;

void record_share(sqlite3 *dbhandle, const std::string &btc, const std::string &a, const std::string &b, const std::string &c)
{
	std::string e;
	earningMutex.lock();
	e = earning;
	earningMutex.unlock();


}

int get_share_count(sqlite3 *dbhandle, const std::string &x)
{
	sqlite3_stmt *stmt = 0;
	const char *ignoreMe = 0;
	sqlite3_prepare_v2(dbhandle, "SELECT COUNT(*) from Shares WHERE X=?", 1024, &stmt, &ignoreMe);
	sqlite3_bind_text(stmt, 1, x.data(), -1, NULL);
	sqlite3_step(stmt);
	int r = sqlite3_column_int(stmt, 0);
	sqlite3_finalize(stmt);
	return r;
}

void serve(CActiveSocket *client)
{
	bool persist = false;
	
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

	std::set<uint32> targets;
	initTargets(targets);

	sqlite3 *dbhandle = 0;
	sqlite3_open_v2("shares.db", &dbhandle, SQLITE_OPEN_READWRITE | SQLITE_OPEN_FULLMUTEX, NULL);

	while (client->Receive(MAX_PACKET))
	{
		char data[1025] = {};
		memcpy(data, client->GetData(), client->GetBytesReceived());
		
		std::string recvd(data);
		if (recvd.compare("persist") == 0)
		{
			persist = true;
			client->Send((const uint8*)"OK\n", 3);
		} else if (recvd.find("result") == 0) {
			std::vector<std::string> tokens;
			std::istringstream iss(recvd);

			std::copy(std::istream_iterator<std::string>(iss),
					 std::istream_iterator<std::string>(),
					 std::back_inserter<std::vector<std::string> >(tokens));

			bool validRequest = tokens.size() == 5;

			if ( validRequest )
			{
				std::string btc(tokens[1]);
				std::string sa(tokens[2].substr(2, tokens[2].length() - 3));
				std::string sb(tokens[3].substr(2, tokens[2].length() - 3));
				std::string sc(tokens[4].substr(2, tokens[2].length() - 3));

				mpz_t a;
				mpz_t b;
				mpz_t c;
				mpz_init_set_str(a, sa.data(), 16);
				mpz_init_set_str(b, sb.data(), 16);
				mpz_init_set_str(c, sc.data(), 16);

				Point myp(g * c);
				
				uint32 match = static_cast<uint32>(mpz_get_ui(myp.x()));

				if (mpz_cmp(myp.x(), a) == 0 && mpz_cmp(myp.y(), b) == 0 && targets.count(match) > 0)
				{
					int howmany = get_share_count(dbhandle, tokens[2]);

					if ( howmany == 0 )
					{
						record_share(dbhandle, btc, sa, sb, sc);
						
						const char *successResponse = "{\"msg\": \"hint\",\"reason\": \"your submission is great.\"}\n";
						client->Send((const uint8*)successResponse, strlen(successResponse));
					} else {
						const char *duplicateResponse = "{\"msg\": \"hint\",\"reason\": \"your submission is already here.\"}\n";
						client->Send((const uint8*)duplicateResponse, strlen(duplicateResponse));
					}
				} else {
					const char *invalidResponse = "{\"msg\": \"hint\",\"reason\": \"your submission is invalid.\"}\n";
					client->Send((const uint8*)invalidResponse, strlen(invalidResponse));
				}
			}
		}

		if (!persist)
			break;
	}
	client->Close();
	sqlite3_close(dbhandle);

	delete client;

	mpz_clear(a);
	mpz_clear(b);
	mpz_clear(p);
	mpz_clear(gy);
	mpz_clear(gx);
	mpz_clear(r);
}

int main(int argc, char **argv)
{
	CPassiveSocket serverSocket;
	CActiveSocket *clientSocket = NULL;

	serverSocket.Initialize();

	serverSocket.Listen((const uint8*)"0.0.0.0", 4444);

	while( true )
	{
		if ((clientSocket = serverSocket.Accept()) != NULL)
			new std::thread(serve, clientSocket);
	}

	return 0;
}