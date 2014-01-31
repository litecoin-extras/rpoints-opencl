class Point
{
public:
	Point(const CurveFp &curve, const mpz_t &x, const mpz_t &y, const mpz_t &o)
		:m_curve(curve)
	{
		mpz_init_set(m_x, x);
		mpz_init_set(m_y, y);
		mpz_init_set(m_o, o);
		m_inv = false;

		/*assert(curve.contains_point(x, y));
		if (mpz_cmp_ui(o, 0) != 0)
		{
			//assert(Point((*this) * o).m_inv); <-- infinite recursion!
		}*/
		init_sp();
	}
	Point(const Point &other)
		:m_curve(other.m_curve)
	{
		mpz_init_set(m_x, other.m_x);
		mpz_init_set(m_y, other.m_y);
		mpz_init_set(m_o, other.m_o);
		m_inv = other.m_inv;
		init_sp();
	}
	Point()
	{
		mpz_init(m_x);
		mpz_init(m_y);
		mpz_init(m_o);
		m_inv = true;
		init_sp();
	}
	virtual ~Point()
	{
		mpz_clear(m_x);
		mpz_clear(m_y);
		mpz_clear(m_o);
		done_sp();
	}

	// Scratchpad variables
	mpz_t imod;
	mpz_t l;
	mpz_t lsq;
	mpz_t lsq2;
	mpz_t lsq3;
	mpz_t mul;
	mpz_t x3;
	mpz_t xmx3;
	mpz_t xsub;
	mpz_t ymul;
	mpz_t ys;
	mpz_t ysub;

	void init_sp(){
		mpz_init(imod);
		mpz_init(l);
		mpz_init(lsq);
		mpz_init(lsq2);
		mpz_init(lsq3);
		mpz_init(mul);
		mpz_init(x3);
		mpz_init(xmx3);
		mpz_init(xsub);
		mpz_init(ymul);
		mpz_init(ys);
		mpz_init(ysub);
	}

	void done_sp(){
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
	void add(const Point &other)
	{
		if (other.m_inv)
			return;
		if (m_inv)
		{
			*this = other;
			return;
		}

		if (mpz_cmp(m_x, other.m_x) == 0)
		{
			mpz_t yadd;
			mpz_t ymod;
			mpz_init(yadd);
			mpz_init(ymod);

			mpz_add(yadd, m_y, other.m_y);
			mpz_mod(ymod, yadd, m_curve.m_p);

			int cmp = mpz_cmp_ui(ymod, 0);

			mpz_clear(yadd);
			mpz_clear(ymod);

			if (cmp == 0)
			{
				*this = Point();
				return;
			}
			else
			{
				*this = dbl();
				return;
			}
		}

		mpz_sub(ysub, other.m_y, m_y);
		mpz_sub(xsub, other.m_x, m_x);
		mpz_invert(imod, xsub, m_curve.m_p);
		mpz_mul(mul, ysub, imod);
		mpz_mod(l, mul, m_curve.m_p);
		mpz_pow_ui(lsq, l, 2);
		mpz_sub(lsq2, lsq, m_x);
		mpz_sub(lsq3, lsq2, other.m_x);
		mpz_mod(x3, lsq3, m_curve.m_p);
		mpz_sub(xmx3, m_x, x3);
		mpz_mul(ymul, l, xmx3);
		mpz_sub(ys, ymul, m_y);
		mpz_mod(m_y, ys, m_curve.m_p);
		mpz_set(m_x, x3);
		mpz_set_ui(m_o, 0);
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
				r.add(*this);
			if (mpz_cmp_ui(e3i, 0) == 0 && mpz_cmp_ui(ei, 0) != 0)
				r.add(negative_self);
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
		mpz_invert(imod, ydbl, m_curve.m_p);

		mpz_t mul;
		mpz_init(mul);
		mpz_mul(mul, xsq3a, imod);

		mpz_t l;
		mpz_init(l);
		mpz_mod(l, mul, m_curve.m_p);

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
		mpz_mod(x3, lsqmx2, m_curve.m_p);

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
		mpz_mod(y3, ys, m_curve.m_p);

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
