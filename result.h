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