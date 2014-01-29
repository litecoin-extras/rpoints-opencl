
class CurveFp
{
public:
	CurveFp(const mpz_t &p, const mpz_t &a, const mpz_t b)
	{
		mpz_init_set(m_p, p);
		mpz_init_set(m_a, a);
		mpz_init_set(m_b, b);
		m_inv = false;
	}
	CurveFp(const CurveFp &other)
	{
		mpz_init_set(m_p, other.m_p);
		mpz_init_set(m_a, other.m_a);
		mpz_init_set(m_b, other.m_b);
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
