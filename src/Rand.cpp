#include <stdio.h>
#include <cmath>
#include <algorithm>
#include "Rand.h"
#include "MachineDefines.h"
#include "Constants.h"
#include "Error.h"
#ifdef _OPENMP
#include <omp.h>
#endif

/* RANDLIB static variables */
int32_t* Xcg1, *Xcg2;
int **SamplingQueue = nullptr;

///////////// ********* ///////////// ********* ///////////// ********* ///////////// ********* ///////////// ********* ///////////// *********
/////////////////////// NEIL rand_lib code (with some Gemma rand lib also)
///////////// ********* ///////////// ********* ///////////// ********* ///////////// ********* ///////////// ********* ///////////// *********

double ranf(void)
{
	int32_t k, s1, s2, z;
	unsigned int curntg;

#ifdef _OPENMP
	curntg = CACHE_LINE_SIZE * omp_get_thread_num();
#else
	curntg = 0;
#endif
	s1 = Xcg1[curntg];
	s2 = Xcg2[curntg];
	k = s1 / 53668;
	s1 = Xa1 * (s1 - k * 53668) - k * 12211;
	if (s1 < 0) s1 += Xm1;
	k = s2 / 52774;
	s2 = Xa2 * (s2 - k * 52774) - k * 3791;
	if (s2 < 0) s2 += Xm2;
	Xcg1[curntg] = s1;
	Xcg2[curntg] = s2;
	z = s1 - s2;
	if (z < 1) z += (Xm1 - 1);
	return ((double)z) / Xm1;
}
double ranf_mt(int tn)
{
	int32_t k, s1, s2, z;
	int curntg;

	curntg = CACHE_LINE_SIZE * tn;
	s1 = Xcg1[curntg];
	s2 = Xcg2[curntg];
	k = s1 / 53668;
	s1 = Xa1 * (s1 - k * 53668) - k * 12211;
	if (s1 < 0) s1 += Xm1;
	k = s2 / 52774;
	s2 = Xa2 * (s2 - k * 52774) - k * 3791;
	if (s2 < 0) s2 += Xm2;
	Xcg1[curntg] = s1;
	Xcg2[curntg] = s2;
	z = s1 - s2;
	if (z < 1) z += (Xm1 - 1);
	return ((double)z) / Xm1;
}

void setall(int32_t *pseed1, int32_t *pseed2)
/*
**********************************************************************
	 void setall(int32_t iseed1,int32_t iseed2)
			   SET ALL random number generators
	 Sets the initial seed of generator 1 to ISEED1 and ISEED2. The
	 initial seeds of the other generators are set accordingly, and
	 all generators states are set to these seeds.
	 This is a transcription from Pascal to C++ of routine
	 Set_Initial_Seed from the paper
	 L'Ecuyer, P. and Cote, S. "Implementing a Random Number Package
	 with Splitting Facilities." ACM Transactions on Mathematical
	 Software, 17:98-111 (1991)
	 https://citeseerx.ist.psu.edu/viewdoc/summary?doi=10.1.1.149.9439
							  Arguments
	 iseed1 -> First of two integer seeds
	 iseed2 -> Second of two integer seeds
**********************************************************************
*/
{
	int g;

	int32_t iseed1 = *pseed1;
	int32_t iseed2 = *pseed2;

	for (g = 0; g < MAX_NUM_THREADS; g++) {
		*(Xcg1 + g * CACHE_LINE_SIZE) = iseed1 = mltmod(Xa1vw, iseed1, Xm1);
		*(Xcg2 + g * CACHE_LINE_SIZE) = iseed2 = mltmod(Xa2vw, iseed2, Xm2);
	}

	*pseed1 = iseed1;
	*pseed2 = iseed2;
}
int32_t mltmod(int32_t a, int32_t s, int32_t m)
/*
**********************************************************************
	 int32_t mltmod(int32_t a, int32_t s, int32_t m)
					Returns (a * s) MOD m
	 This is a transcription from Pascal to C++ of routine
	 MULtMod_Decompos from the paper
	 L'Ecuyer, P. and Cote, S. "Implementing a Random Number Package
	 with Splitting Facilities." ACM Transactions on Mathematical
	 Software, 17:98-111 (1991)
	https://citeseerx.ist.psu.edu/viewdoc/summary?doi=10.1.1.149.9439
**********************************************************************
*/
{
	const int32_t h = 32768;
	int32_t a0, a1, k, p, q, qh, rh;
	/*
		 H = 2**((b-2)/2) where b = 32 because we are using a 32 bit
		  machine. On a different machine recompute H
	*/
	if (a <= 0 || a >= m || s <= 0 || s >= m) {
		fputs(" a, m, s out of order in mltmod - ABORT!\n", stderr);
		fprintf(stderr, " a = %12d s = %12d m = %12d\n", a, s, m);
		fputs(" mltmod requires: 0 < a < m; 0 < s < m\n", stderr);
		exit(1);
	}

	if (a < h) {
		a0 = a;
		p = 0;
	} else {
		a1 = a / h;
		a0 = a - h * a1;
		qh = m / h;
		rh = m - h * qh;
		if (a1 >= h) { // a2 == 1
			a1 -= h;
			k = s / qh;
			p = h * (s - k * qh) - k * rh;
			while (p < 0) {
				p += m;
			}
		} else {
			p = 0;
		}
		// p == (a2 * s * h) MOD m
		if (a1 != 0) {
			q = m / a1;
			k = s / q;
			p -= (k * (m - a1 * q));
			if (p > 0) p -= m;
			p += (a1 * (s - k * q));
			while (p < 0) {
				p += m;
			}
		}
		// p == ((a2 * h + a1) * s) MOD m
		k = p / qh;
		p = h * (p - k * qh) - k * rh;
		while (p < 0) {
			p += m;
		}
	}
	// p == ((a2 * h + a1) * h * s) MOD m
	if (a0 != 0) {
		q = m / a0;
		k = s / q;
		p -= (k * (m - a0 * q));
		if (p > 0) p -= m;
		p += (a0 * (s - k * q));
		while (p < 0) {
			p += m;
		}
	}
	return p;
}

int32_t ignbin(int32_t n, double pp)
{
	/*
**********************************************************************
	 int32_t ignbin(int32_t n,double pp)
					GENerate BINomial random deviate
							  Function
	 Generates a single random deviate from a binomial
	 distribution whose number of trials is N and whose
	 probability of an event in each trial is P.
							  Arguments
	 n  --> The number of trials in the binomial distribution
			from which a random deviate is to be generated.
		JJV (N >= 0)
	 pp --> The probability of an event in each trial of the
			binomial distribution from which a random deviate
			is to be generated.
		JJV (0.0 <= PP <= 1.0)
	 ignbin <-- A random deviate yielding the number of events
				from N independent trials, each of which has
				a probability of event P.
							  Method
	 This is algorithm BTPE from:
		 Kachitvichyanukul, V. and Schmeiser, B. W.
		 Binomial Random Variate Generation.
		 Communications of the ACM, 31, 2
		 (February, 1988) 216.
**********************************************************************
	 SUBROUTINE BTPEC(N,PP,ISEED,JX)
	 BINOMIAL RANDOM VARIATE GENERATOR
	 MEAN .LT. 30 -- INVERSE CDF
	   MEAN .GE. 30 -- ALGORITHM BTPE:  ACCEPTANCE-REJECTION VIA
	   FOUR REGION COMPOSITION.  THE FOUR REGIONS ARE A TRIANGLE
	   (SYMMETRIC IN THE CENTER), A PAIR OF PARALLELOGRAMS (ABOVE
	   THE TRIANGLE), AND EXPONENTIAL LEFT AND RIGHT TAILS.
	 BTPE REFERS TO BINOMIAL-TRIANGLE-PARALLELOGRAM-EXPONENTIAL.
	 BTPEC REFERS TO BTPE AND "COMBINED."  THUS BTPE IS THE
	   RESEARCH AND BTPEC IS THE IMPLEMENTATION OF A COMPLETE
	   USABLE ALGORITHM.
	 REFERENCE:  VORATAS KACHITVICHYANUKUL AND BRUCE SCHMEISER,
	   "BINOMIAL RANDOM VARIATE GENERATION,"
	   COMMUNICATIONS OF THE ACM, FORTHCOMING
	 WRITTEN:  SEPTEMBER 1980.
	   LAST REVISED:  MAY 1985, JULY 1987
	 REQUIRED SUBPROGRAM:  RAND() -- A UNIFORM (0,1) RANDOM NUMBER
						   GENERATOR
	 ARGUMENTS
	   N : NUMBER OF BERNOULLI TRIALS            (INPUT)
	   PP : PROBABILITY OF SUCCESS IN EACH TRIAL (INPUT)
	   ISEED:  RANDOM NUMBER SEED                (INPUT AND OUTPUT)
	   JX:  RANDOMLY GENERATED OBSERVATION       (OUTPUT)
	 VARIABLES
	   PSAVE: VALUE OF PP FROM THE LAST CALL TO BTPEC
	   NSAVE: VALUE OF N FROM THE LAST CALL TO BTPEC
	   XNP:  VALUE OF THE MEAN FROM THE LAST CALL TO BTPEC
	   P: PROBABILITY USED IN THE GENERATION PHASE OF BTPEC
	   FFM: TEMPORARY VARIABLE EQUAL TO XNP + P
	   M:  INTEGER VALUE OF THE CURRENT MODE
	   FM:  FLOATING POINT VALUE OF THE CURRENT MODE
	   XNPQ: TEMPORARY VARIABLE USED IN SETUP AND SQUEEZING STEPS
	   P1:  AREA OF THE TRIANGLE
	   C:  HEIGHT OF THE PARALLELOGRAMS
	   XM:  CENTER OF THE TRIANGLE
	   XL:  LEFT END OF THE TRIANGLE
	   XR:  RIGHT END OF THE TRIANGLE
	   AL:  TEMPORARY VARIABLE
	   XLL:  RATE FOR THE LEFT EXPONENTIAL TAIL
	   XLR:  RATE FOR THE RIGHT EXPONENTIAL TAIL
	   P2:  AREA OF THE PARALLELOGRAMS
	   P3:  AREA OF THE LEFT EXPONENTIAL TAIL
	   P4:  AREA OF THE RIGHT EXPONENTIAL TAIL
	   U:  A U(0,P4) RANDOM VARIATE USED FIRST TO SELECT ONE OF THE
		   FOUR REGIONS AND THEN CONDITIONALLY TO GENERATE A VALUE
		   FROM THE REGION
	   V:  A U(0,1) RANDOM NUMBER USED TO GENERATE THE RANDOM VALUE
		   (REGION 1) OR TRANSFORMED INTO THE VARIATE TO ACCEPT OR
		   REJECT THE CANDIDATE VALUE
	   IX:  INTEGER CANDIDATE VALUE
	   X:  PRELIMINARY CONTINUOUS CANDIDATE VALUE IN REGION 2 LOGIC
		   AND A FLOATING POINT IX IN THE ACCEPT/REJECT LOGIC
	   K:  ABSOLUTE VALUE OF (IX-M)
	   F:  THE HEIGHT OF THE SCALED DENSITY FUNCTION USED IN THE
		   ACCEPT/REJECT DECISION WHEN BOTH M AND IX ARE SMALL
		   ALSO USED IN THE INVERSE TRANSFORMATION
	   R: THE RATIO P/Q
	   G: CONSTANT USED IN CALCULATION OF PROBABILITY
	   MP:  MODE PLUS ONE, THE LOWER INDEX FOR EXPLICIT CALCULATION
			OF F WHEN IX IS GREATER THAN M
	   IX1:  CANDIDATE VALUE PLUS ONE, THE LOWER INDEX FOR EXPLICIT
			 CALCULATION OF F WHEN IX IS LESS THAN M
	   I:  INDEX FOR EXPLICIT CALCULATION OF F FOR BTPE
	   AMAXP: MAXIMUM ERROR OF THE LOGARITHM OF NORMAL BOUND
	   YNORM: LOGARITHM OF NORMAL BOUND
	   ALV:  NATURAL LOGARITHM OF THE ACCEPT/REJECT VARIATE V
	   X1,F1,Z,W,Z2,X2,F2, AND W2 ARE TEMPORARY VARIABLES TO BE
	   USED IN THE FINAL ACCEPT/REJECT TEST
	   QN: PROBABILITY OF NO SUCCESS IN N TRIALS
	 REMARK
	   IX AND JX COULD LOGICALLY BE THE SAME VARIABLE, WHICH WOULD
	   SAVE A MEMORY POSITION AND A LINE OF CODE.  HOWEVER, SOME
	   COMPILERS (E.G.,CDC MNF) OPTIMIZE BETTER WHEN THE ARGUMENTS
	   ARE NOT INVOLVED.
	 ISEED NEEDS TO BE DOUBLE PRECISION IF THE IMSL ROUTINE
	 GGUBFS IS USED TO GENERATE UNIFORM RANDOM NUMBER, OTHERWISE
	 TYPE OF ISEED SHOULD BE DICTATED BY THE UNIFORM GENERATOR
**********************************************************************
*****DETERMINE APPROPRIATE ALGORITHM AND WHETHER SETUP IS NECESSARY
*/
/* JJV changed initial values to ridiculous values */
	double psave = -1.0E37;
	int32_t nsave = -214748365;
	int32_t ignbin, i, ix, ix1, k, m, mp, T1;
	double al, alv, amaxp, c, f, f1, f2, ffm, fm, g, p, p1, p2, p3, p4, q, qn, r, u, v, w, w2, x, x1,
		x2, xl, xll, xlr, xm, xnp, xnpq, xr, ynorm, z, z2;

	/*
	*****SETUP, PERFORM ONLY WHEN PARAMETERS CHANGE
	JJV added checks to ensure 0.0 <= PP <= 1.0
	*/
	if (pp < 0.0F) ERR_CRITICAL("PP < 0.0 in IGNBIN");
	if (pp > 1.0F) ERR_CRITICAL("PP > 1.0 in IGNBIN");
	psave = pp;
	p = std::min(psave, 1.0 - psave);
	q = 1.0 - p;

	/*
	JJV added check to ensure N >= 0
	*/
	if (n < 0L) ERR_CRITICAL("N < 0 in IGNBIN");
	xnp = n * p;
	nsave = n;
	if (xnp < 30.0) goto S140;
	ffm = xnp + p;
	m = (int32_t)ffm;
	fm = m;
	xnpq = xnp * q;
	p1 = (int32_t)(2.195 * sqrt(xnpq) - 4.6 * q) + 0.5;
	xm = fm + 0.5;
	xl = xm - p1;
	xr = xm + p1;
	c = 0.134 + 20.5 / (15.3 + fm);
	al = (ffm - xl) / (ffm - xl * p);
	xll = al * (1.0 + 0.5 * al);
	al = (xr - ffm) / (xr * q);
	xlr = al * (1.0 + 0.5 * al);
	p2 = p1 * (1.0 + c + c);
	p3 = p2 + c / xll;
	p4 = p3 + c / xlr;
S30:
	/*
	*****GENERATE VARIATE
	*/
	u = ranf() * p4;
	v = ranf();
	/*
		 TRIANGULAR REGION
	*/
	if (u > p1) goto S40;
	ix = (int32_t)(xm - p1 * v + u);
	goto S170;
S40:
	/*
		 PARALLELOGRAM REGION
	*/
	if (u > p2) goto S50;
	x = xl + (u - p1) / c;
	v = v * c + 1.0 - std::abs(xm - x) / p1;
	if (v > 1.0 || v <= 0.0) goto S30;
	ix = (int32_t)x;
	goto S70;
S50:
	/*
		 LEFT TAIL
	*/
	if (u > p3) goto S60;
	ix = (int32_t)(xl + log(v) / xll);
	if (ix < 0) goto S30;
	v *= ((u - p2) * xll);
	goto S70;
S60:
	/*
		 RIGHT TAIL
	*/
	ix = (int32_t)(xr - log(v) / xlr);
	if (ix > n) goto S30;
	v *= ((u - p3) * xlr);
S70:
	/*
	*****DETERMINE APPROPRIATE WAY TO PERFORM ACCEPT/REJECT TEST
	*/
	k = std::abs(ix - m);
	if (k > 20 && k < xnpq / 2 - 1) goto S130;
	/*
		 EXPLICIT EVALUATION
	*/
	f = 1.0;
	r = p / q;
	g = (n + 1) * r;
	T1 = m - ix;
	if (T1 < 0) goto S80;
	else if (T1 == 0) goto S120;
	else  goto S100;
S80:
	mp = m + 1;
	for (i = mp; i <= ix; i++) f *= (g / i - r);
	goto S120;
S100:
	ix1 = ix + 1;
	for (i = ix1; i <= m; i++) f /= (g / i - r);
S120:
	if (v <= f) goto S170;
	goto S30;
S130:
	/*
		 SQUEEZING USING UPPER AND LOWER BOUNDS ON ALOG(F(X))
	*/
	amaxp = k / xnpq * ((k * (k / 3.0 + 0.625) + 0.1666666666666) / xnpq + 0.5);
	ynorm = -(k * k / (2.0 * xnpq));
	alv = log(v);
	if (alv < ynorm - amaxp) goto S170;
	if (alv > ynorm + amaxp) goto S30;
	/*
		 STIRLING'S FORMULA TO MACHINE ACCURACY FOR
		 THE FINAL ACCEPTANCE/REJECTION TEST
	*/
	x1 = ix + 1.0;
	f1 = fm + 1.0;
	z = n + 1.0 - fm;
	w = n - ix + 1.0;
	z2 = z * z;
	x2 = x1 * x1;
	f2 = f1 * f1;
	w2 = w * w;
	if (alv <= xm * log(f1 / x1) + (n - m + 0.5) * log(z / w) + (ix - m) * log(w * p / (x1 * q)) + (13860.0 -
		(462.0 - (132.0 - (99.0 - 140.0 / f2) / f2) / f2) / f2) / f1 / 166320.0 + (13860.0 - (462.0 -
		(132.0 - (99.0 - 140.0 / z2) / z2) / z2) / z2) / z / 166320.0 + (13860.0 - (462.0 - (132.0 -
			(99.0 - 140.0 / x2) / x2) / x2) / x2) / x1 / 166320.0 + (13860.0 - (462.0 - (132.0 - (99.0
				- 140.0 / w2) / w2) / w2) / w2) / w / 166320.0) goto S170;
	goto S30;
S140:
	/*
		 INVERSE CDF LOGIC FOR MEAN LESS THAN 30
	*/
	qn = pow(q, (double)n);
	r = p / q;
	g = r * (n + 1);
S150:
	ix = 0;
	f = qn;
	u = ranf();
S160:
	if (u < f) goto S170;
	if (ix > 110) goto S150;
	u -= f;
	ix += 1;
	f *= (g / ix - r);
	goto S160;
S170:
	if (psave > 0.5) ix = n - ix;
	ignbin = ix;
	return ignbin;
}
int32_t ignpoi(double mu)
/*
**********************************************************************
	 int32_t ignpoi(double mu)
					GENerate POIsson random deviate
							  Function
	 Generates a single random deviate from a Poisson
	 distribution with mean MU.
							  Arguments
	 mu --> The mean of the Poisson distribution from which
			a random deviate is to be generated.
		(mu >= 0.0)
	 ignpoi <-- The random deviate.
							  Method
	 Renames KPOIS from TOMS as slightly modified by BWB to use RANF
	 instead of SUNIF.
	 For details see:
			   Ahrens, J.H. and Dieter, U.
			   Computer Generation of Poisson Deviates
			   From Modified Normal Distributions.
			   ACM Trans. Math. Software, 8, 2
			   (June 1982),163-179
**********************************************************************
**********************************************************************


	 P O I S S O N  DISTRIBUTION


**********************************************************************
**********************************************************************

	 FOR DETAILS SEE:

			   AHRENS, J.H. AND DIETER, U.
			   COMPUTER GENERATION OF POISSON DEVIATES
			   FROM MODIFIED NORMAL DISTRIBUTIONS.
			   ACM TRANS. MATH. SOFTWARE, 8,2 (JUNE 1982), 163 - 179.

	 (SLIGHTLY MODIFIED VERSION OF THE PROGRAM IN THE ABOVE ARTICLE)

**********************************************************************
	  INTEGER FUNCTION IGNPOI(IR,MU)
	 INPUT:  IR=CURRENT STATE OF BASIC RANDOM NUMBER GENERATOR
			 MU=MEAN MU OF THE POISSON DISTRIBUTION
	 OUTPUT: IGNPOI=SAMPLE FROM THE POISSON-(MU)-DISTRIBUTION
	 MUPREV=PREVIOUS MU, MUOLD=MU AT LAST EXECUTION OF STEP P OR B.
	 TABLES: COEFFICIENTS A0-A7 FOR STEP F. FACTORIALS FACT
	 COEFFICIENTS A(K) - FOR PX = FK*V*V*SUM(A(K)*V**K)-DEL
	 SEPARATION OF CASES A AND B
*/
{
	extern double fsign(double num, double sign);
	static double a0 = -0.5;
	static double a1 = 0.3333333;
	static double a2 = -0.2500068;
	static double a3 = 0.2000118;
	static double a4 = -0.1661269;
	static double a5 = 0.1421878;
	static double a6 = -0.1384794;
	static double a7 = 0.125006;
	/* JJV changed the initial values of MUPREV and MUOLD */
	static double fact[10] = {
		1.0,1.0,2.0,6.0,24.0,120.0,720.0,5040.0,40320.0,362880.0
	};
	/* JJV added ll to the list, for Case A */
	int32_t ignpoi, j, k, kflag, l, ll, m;
	double b1, b2, c, c0, c1, c2, c3, d, del, difmuk, e, fk, fx, fy, g, omega, p, p0, px, py, q, s, t, u, v, x, xx, pp[35];

	if (mu < 10.0) goto S120;
	/*
		 C A S E  A. (RECALCULATION OF S,D,LL IF MU HAS CHANGED)
		 JJV changed l in Case A to ll
	*/
	s = sqrt(mu);
	d = 6.0 * mu * mu;
	/*
				 THE POISSON PROBABILITIES PK EXCEED THE DISCRETE NORMAL
				 PROBABILITIES FK WHENEVER K >= M(MU). LL=IFIX(MU-1.1484)
				 IS AN UPPER BOUND TO M(MU) FOR ALL MU >= 10 .
	*/
	ll = (int32_t)(mu - 1.1484);

	/*
		 STEP N. NORMAL SAMPLE - SNORM(IR) FOR STANDARD NORMAL DEVIATE
	*/
	g = mu + s * snorm();
	if (g < 0.0) goto S20;
	ignpoi = (int32_t)(g);
	/*
		 STEP I. IMMEDIATE ACCEPTANCE IF IGNPOI IS LARGE ENOUGH
	*/
	if (ignpoi >= ll) return ignpoi;
	/*
		 STEP S. SQUEEZE ACCEPTANCE - SUNIF(IR) FOR (0,1)-SAMPLE U
	*/
	fk = (double)ignpoi;
	difmuk = mu - fk;
	u = ranf();
	if (d * u >= difmuk * difmuk * difmuk) return ignpoi;
S20:
	/*
		 STEP P. PREPARATIONS FOR STEPS Q AND H.
				 (RECALCULATIONS OF PARAMETERS IF NECESSARY)
				 .3989423=(2*PI)**(-.5)  .416667E-1=1./24.  .1428571=1./7.
				 THE QUANTITIES B1, B2, C3, C2, C1, C0 ARE FOR THE HERMITE
				 APPROXIMATIONS TO THE DISCRETE NORMAL PROBABILITIES FK.
				 C=.1069/MU GUARANTEES MAJORIZATION BY THE 'HAT'-FUNCTION.
	*/
	omega = 0.3989423 / s;
	b1 = 4.166667E-2 / mu;
	b2 = 0.3 * b1 * b1;
	c3 = 0.1428571 * b1 * b2;
	c2 = b2 - 15.0 * c3;
	c1 = b1 - 6.0 * b2 + 45.0 * c3;
	c0 = 1.0 - b1 + 3.0 * b2 - 15.0 * c3;
	c = 0.1069 / mu;

	if (g < 0.0) goto S50;
	/*
				 'SUBROUTINE' F IS CALLED (KFLAG=0 FOR CORRECT RETURN)
	*/
	kflag = 0;
	goto S70;
S40:
	/*
		 STEP Q. QUOTIENT ACCEPTANCE (RARE CASE)
	*/
	if (fy - u * fy <= py * exp(px - fx)) return ignpoi;
S50:
	/*
		 STEP E. EXPONENTIAL SAMPLE - SEXPO(IR) FOR STANDARD EXPONENTIAL
				 DEVIATE E AND SAMPLE T FROM THE LAPLACE 'HAT'
				 (IF T <= -.6744 THEN PK < FK FOR ALL MU >= 10.)
	*/
	e = sexpo();
	u = ranf();
	u += (u - 1.0);
	t = 1.8 + fsign(e, u);
	if (t <= -0.6744) goto S50;
	ignpoi = (int32_t)(mu + s * t);
	fk = (double)ignpoi;
	difmuk = mu - fk;
	/*
				 'SUBROUTINE' F IS CALLED (KFLAG=1 FOR CORRECT RETURN)
	*/
	kflag = 1;
	goto S70;
S60:
	/*
		 STEP H. HAT ACCEPTANCE (E IS REPEATED ON REJECTION)
	*/
	if (c * fabs(u) > py * exp(px + e) - fy * exp(fx + e)) goto S50;
	return ignpoi;
S70:
	/*
		 STEP F. 'SUBROUTINE' F. CALCULATION OF PX,PY,FX,FY.
				 CASE IGNPOI .LT. 10 USES FACTORIALS FROM TABLE FACT
	*/
	if (ignpoi >= 10) goto S80;
	px = -mu;
	py = pow(mu, (double)ignpoi) / *(fact + ignpoi);
	goto S110;
S80:
	/*
				 CASE IGNPOI .GE. 10 USES POLYNOMIAL APPROXIMATION
				 A0-A7 FOR ACCURACY WHEN ADVISABLE
				 .8333333E-1=1./12.  .3989423=(2*PI)**(-.5)
	*/
	del = 8.333333E-2 / fk;
	del -= (4.8 * del * del * del);
	v = difmuk / fk;
	if (fabs(v) <= 0.25) goto S90;
	px = fk * log(1.0 + v) - difmuk - del;
	goto S100;
S90:
	px = fk * v * v * (((((((a7 * v + a6) * v + a5) * v + a4) * v + a3) * v + a2) * v + a1) * v + a0) - del;
S100:
	py = 0.3989423 / sqrt(fk);
S110:
	x = (0.5 - difmuk) / s;
	xx = x * x;
	fx = -0.5 * xx;
	fy = omega * (((c3 * xx + c2) * xx + c1) * xx + c0);
	if (kflag <= 0) goto S40;
	goto S60;
S120:
	/*
		 C A S E  B. (START NEW TABLE AND CALCULATE P0 IF NECESSARY)
		 JJV changed MUPREV assignment to initial value
	*/
	m = std::max(INT32_C(1), (int32_t)(mu));
	l = 0;
	p = exp(-mu);
	q = p0 = p;
S130:
	/*
		 STEP U. UNIFORM SAMPLE FOR INVERSION METHOD
	*/
	u = ranf();
	ignpoi = 0;
	if (u <= p0) return ignpoi;
	/*
		 STEP T. TABLE COMPARISON UNTIL THE END PP(L) OF THE
				 PP-TABLE OF CUMULATIVE POISSON PROBABILITIES
				 (0.458=PP(9) FOR MU=10)
	*/
	if (l == 0) goto S150;
	j = 1;
	if (u > 0.458) j = std::min(l, m);
	for (k = j; k <= l; k++) {
		if (u <= *(pp + k - 1)) goto S180;
	}
	if (l == 35) goto S130;
S150:
	/*
		 STEP C. CREATION OF NEW POISSON PROBABILITIES P
				 AND THEIR CUMULATIVES Q=PP(K)
	*/
	l += 1;
	for (k = l; k <= 35; k++) {
		p = p * mu / (double)k;
		q += p;
		*(pp + k - 1) = q;
		if (u <= q) goto S170;
	}
	l = 35;
	goto S130;
S170:
	l = k;
S180:
	ignpoi = k;
	return ignpoi;
}
int32_t ignpoi_mt(double mu, int tn)
/*
**********************************************************************
int32_t ignpoi_mt(double mu)
GENerate POIsson random deviate
Function
Generates a single random deviate from a Poisson
distribution with mean MU.
Arguments
mu --> The mean of the Poisson distribution from which
a random deviate is to be generated.
(mu >= 0.0)
ignpoi_mt <-- The random deviate.
Method
Renames KPOIS from TOMS as slightly modified by BWB to use RANF
instead of SUNIF.
For details see:
Ahrens, J.H. and Dieter, U.
Computer Generation of Poisson Deviates
From Modified Normal Distributions.
ACM Trans. Math. Software, 8, 2
(June 1982),163-179
**********************************************************************
**********************************************************************


P O I S S O N  DISTRIBUTION


**********************************************************************
**********************************************************************

FOR DETAILS SEE:

AHRENS, J.H. AND DIETER, U.
COMPUTER GENERATION OF POISSON DEVIATES
FROM MODIFIED NORMAL DISTRIBUTIONS.
ACM TRANS. MATH. SOFTWARE, 8,2 (JUNE 1982), 163 - 179.

(SLIGHTLY MODIFIED VERSION OF THE PROGRAM IN THE ABOVE ARTICLE)

**********************************************************************
INTEGER FUNCTION IGNPOI(IR,MU)
INPUT:  IR=CURRENT STATE OF BASIC RANDOM NUMBER GENERATOR
MU=MEAN MU OF THE POISSON DISTRIBUTION
OUTPUT: IGNPOI=SAMPLE FROM THE POISSON-(MU)-DISTRIBUTION
MUPREV=PREVIOUS MU, MUOLD=MU AT LAST EXECUTION OF STEP P OR B.
TABLES: COEFFICIENTS A0-A7 FOR STEP F. FACTORIALS FACT
COEFFICIENTS A(K) - FOR PX = FK*V*V*SUM(A(K)*V**K)-DEL
SEPARATION OF CASES A AND B
*/
{
	extern double fsign(double num, double sign);
	double a0 = -0.5;
	double a1 = 0.3333333;
	double a2 = -0.2500068;
	double a3 = 0.2000118;
	double a4 = -0.1661269;
	double a5 = 0.1421878;
	double a6 = -0.1384794;
	double a7 = 0.125006;
	/* JJV changed the initial values of MUPREV and MUOLD */
	double fact[10] = {
		1.0,1.0,2.0,6.0,24.0,120.0,720.0,5040.0,40320.0,362880.0
	};
	/* JJV added ll to the list, for Case A */
	int32_t ignpoi_mt, j, k, kflag, l, ll, m;
	double b1, b2, c, c0, c1, c2, c3, d, del, difmuk, e, fk, fx, fy, g, omega, p, p0, px, py, q, s, t, u, v, x, xx, pp[35];

	if (mu < 10.0) goto S120;
	/*
	C A S E  A. (RECALCULATION OF S,D,LL IF MU HAS CHANGED)
	JJV changed l in Case A to ll
	*/
	s = sqrt(mu);
	d = 6.0 * mu * mu;
	/*
	THE POISSON PROBABILITIES PK EXCEED THE DISCRETE NORMAL
	PROBABILITIES FK WHENEVER K >= M(MU). LL=IFIX(MU-1.1484)
	IS AN UPPER BOUND TO M(MU) FOR ALL MU >= 10 .
	*/
	ll = (int32_t)(mu - 1.1484);

	/*
	STEP N. NORMAL SAMPLE - SNORM(IR) FOR STANDARD NORMAL DEVIATE
	*/
	g = mu + s * snorm_mt(tn);
	if (g < 0.0) goto S20;
	ignpoi_mt = (int32_t)(g);
	/*
	STEP I. IMMEDIATE ACCEPTANCE IF IGNPOI IS LARGE ENOUGH
	*/
	if (ignpoi_mt >= ll) return ignpoi_mt;
	/*
	STEP S. SQUEEZE ACCEPTANCE - SUNIF(IR) FOR (0,1)-SAMPLE U
	*/
	fk = (double)ignpoi_mt;
	difmuk = mu - fk;
	u = ranf_mt(tn);
	if (d * u >= difmuk * difmuk * difmuk) return ignpoi_mt;
S20:
	/*
	STEP P. PREPARATIONS FOR STEPS Q AND H.
	(RECALCULATIONS OF PARAMETERS IF NECESSARY)
	.3989423=(2*PI)**(-.5)  .416667E-1=1./24.  .1428571=1./7.
	THE QUANTITIES B1, B2, C3, C2, C1, C0 ARE FOR THE HERMITE
	APPROXIMATIONS TO THE DISCRETE NORMAL PROBABILITIES FK.
	C=.1069/MU GUARANTEES MAJORIZATION BY THE 'HAT'-FUNCTION.
	*/
	omega = 0.3989423 / s;
	b1 = 4.166667E-2 / mu;
	b2 = 0.3 * b1 * b1;
	c3 = 0.1428571 * b1 * b2;
	c2 = b2 - 15.0 * c3;
	c1 = b1 - 6.0 * b2 + 45.0 * c3;
	c0 = 1.0 - b1 + 3.0 * b2 - 15.0 * c3;
	c = 0.1069 / mu;

	if (g < 0.0) goto S50;
	/*
	'SUBROUTINE' F IS CALLED (KFLAG=0 FOR CORRECT RETURN)
	*/
	kflag = 0;
	goto S70;
S40:
	/*
	STEP Q. QUOTIENT ACCEPTANCE (RARE CASE)
	*/
	if (fy - u * fy <= py * exp(px - fx)) return ignpoi_mt;
S50:
	/*
	STEP E. EXPONENTIAL SAMPLE - SEXPO(IR) FOR STANDARD EXPONENTIAL
	DEVIATE E AND SAMPLE T FROM THE LAPLACE 'HAT'
	(IF T <= -.6744 THEN PK < FK FOR ALL MU >= 10.)
	*/
	e = sexpo_mt(tn);
	u = ranf_mt(tn);
	u += (u - 1.0);
	t = 1.8 + fsign(e, u);
	if (t <= -0.6744) goto S50;
	ignpoi_mt = (int32_t)(mu + s * t);
	fk = (double)ignpoi_mt;
	difmuk = mu - fk;
	/*
	'SUBROUTINE' F IS CALLED (KFLAG=1 FOR CORRECT RETURN)
	*/
	kflag = 1;
	goto S70;
S60:
	/*
	STEP H. HAT ACCEPTANCE (E IS REPEATED ON REJECTION)
	*/
	if (c * fabs(u) > py * exp(px + e) - fy * exp(fx + e)) goto S50;
	return ignpoi_mt;
S70:
	/*
	STEP F. 'SUBROUTINE' F. CALCULATION OF PX,PY,FX,FY.
	CASE IGNPOI .LT. 10 USES FACTORIALS FROM TABLE FACT
	*/
	if (ignpoi_mt >= 10) goto S80;
	px = -mu;
	py = pow(mu, (double)ignpoi_mt) / *(fact + ignpoi_mt);
	goto S110;
S80:
	/*
	CASE IGNPOI .GE. 10 USES POLYNOMIAL APPROXIMATION
	A0-A7 FOR ACCURACY WHEN ADVISABLE
	.8333333E-1=1./12.  .3989423=(2*PI)**(-.5)
	*/
	del = 8.333333E-2 / fk;
	del -= (4.8 * del * del * del);
	v = difmuk / fk;
	if (fabs(v) <= 0.25) goto S90;
	px = fk * log(1.0 + v) - difmuk - del;
	goto S100;
S90:
	px = fk * v * v * (((((((a7 * v + a6) * v + a5) * v + a4) * v + a3) * v + a2) * v + a1) * v + a0) - del;
S100:
	py = 0.3989423 / sqrt(fk);
S110:
	x = (0.5 - difmuk) / s;
	xx = x * x;
	fx = -0.5 * xx;
	fy = omega * (((c3 * xx + c2) * xx + c1) * xx + c0);
	if (kflag <= 0) goto S40;
	goto S60;
S120:
	/*
	C A S E  B. (START NEW TABLE AND CALCULATE P0 IF NECESSARY)
	JJV changed MUPREV assignment to initial value
	*/
	m = std::max(INT32_C(1), (int32_t)(mu));
	l = 0;
	p = exp(-mu);
	q = p0 = p;
S130:
	/*
	STEP U. UNIFORM SAMPLE FOR INVERSION METHOD
	*/
	u = ranf_mt(tn);
	ignpoi_mt = 0;
	if (u <= p0) return ignpoi_mt;
	/*
	STEP T. TABLE COMPARISON UNTIL THE END PP(L) OF THE
	PP-TABLE OF CUMULATIVE POISSON PROBABILITIES
	(0.458=PP(9) FOR MU=10)
	*/
	if (l == 0) goto S150;
	j = 1;
	if (u > 0.458) j = std::min(l, m);
	for (k = j; k <= l; k++) {
		if (u <= *(pp + k - 1)) goto S180;
	}
	if (l == 35) goto S130;
S150:
	/*
	STEP C. CREATION OF NEW POISSON PROBABILITIES P
	AND THEIR CUMULATIVES Q=PP(K)
	*/
	l += 1;
	for (k = l; k <= 35; k++) {
		p = p * mu / (double)k;
		q += p;
		*(pp + k - 1) = q;
		if (u <= q) goto S170;
	}
	l = 35;
	goto S130;
S170:
	l = k;
S180:
	ignpoi_mt = k;
	return ignpoi_mt;
}
int32_t ignbin_mt(int32_t n, double pp, int tn)
{
	/*
**********************************************************************
int32_t ignbin_mt(int32_t n,double pp)
GENerate BINomial random deviate
Function
Generates a single random deviate from a binomial
distribution whose number of trials is N and whose
probability of an event in each trial is P.
Arguments
n  --> The number of trials in the binomial distribution
from which a random deviate is to be generated.
JJV (N >= 0)
pp --> The probability of an event in each trial of the
binomial distribution from which a random deviate
is to be generated.
JJV (0.0 <= PP <= 1.0)
ignbin <-- A random deviate yielding the number of events
from N independent trials, each of which has
a probability of event P.
Method
This is algorithm BTPE from:
Kachitvichyanukul, V. and Schmeiser, B. W.
Binomial Random Variate Generation.
Communications of the ACM, 31, 2
(February, 1988) 216.
**********************************************************************
SUBROUTINE BTPEC(N,PP,ISEED,JX)
BINOMIAL RANDOM VARIATE GENERATOR
MEAN .LT. 30 -- INVERSE CDF
MEAN .GE. 30 -- ALGORITHM BTPE:  ACCEPTANCE-REJECTION VIA
FOUR REGION COMPOSITION.  THE FOUR REGIONS ARE A TRIANGLE
(SYMMETRIC IN THE CENTER), A PAIR OF PARALLELOGRAMS (ABOVE
THE TRIANGLE), AND EXPONENTIAL LEFT AND RIGHT TAILS.
BTPE REFERS TO BINOMIAL-TRIANGLE-PARALLELOGRAM-EXPONENTIAL.
BTPEC REFERS TO BTPE AND "COMBINED."  THUS BTPE IS THE
RESEARCH AND BTPEC IS THE IMPLEMENTATION OF A COMPLETE
USABLE ALGORITHM.
REFERENCE:  VORATAS KACHITVICHYANUKUL AND BRUCE SCHMEISER,
"BINOMIAL RANDOM VARIATE GENERATION,"
COMMUNICATIONS OF THE ACM, FORTHCOMING
WRITTEN:  SEPTEMBER 1980.
LAST REVISED:  MAY 1985, JULY 1987
REQUIRED SUBPROGRAM:  RAND() -- A UNIFORM (0,1) RANDOM NUMBER
GENERATOR
ARGUMENTS
N : NUMBER OF BERNOULLI TRIALS            (INPUT)
PP : PROBABILITY OF SUCCESS IN EACH TRIAL (INPUT)
ISEED:  RANDOM NUMBER SEED                (INPUT AND OUTPUT)
JX:  RANDOMLY GENERATED OBSERVATION       (OUTPUT)
VARIABLES
PSAVE: VALUE OF PP FROM THE LAST CALL TO BTPEC
NSAVE: VALUE OF N FROM THE LAST CALL TO BTPEC
XNP:  VALUE OF THE MEAN FROM THE LAST CALL TO BTPEC
P: PROBABILITY USED IN THE GENERATION PHASE OF BTPEC
FFM: TEMPORARY VARIABLE EQUAL TO XNP + P
M:  INTEGER VALUE OF THE CURRENT MODE
FM:  FLOATING POINT VALUE OF THE CURRENT MODE
XNPQ: TEMPORARY VARIABLE USED IN SETUP AND SQUEEZING STEPS
P1:  AREA OF THE TRIANGLE
C:  HEIGHT OF THE PARALLELOGRAMS
XM:  CENTER OF THE TRIANGLE
XL:  LEFT END OF THE TRIANGLE
XR:  RIGHT END OF THE TRIANGLE
AL:  TEMPORARY VARIABLE
XLL:  RATE FOR THE LEFT EXPONENTIAL TAIL
XLR:  RATE FOR THE RIGHT EXPONENTIAL TAIL
P2:  AREA OF THE PARALLELOGRAMS
P3:  AREA OF THE LEFT EXPONENTIAL TAIL
P4:  AREA OF THE RIGHT EXPONENTIAL TAIL
U:  A U(0,P4) RANDOM VARIATE USED FIRST TO SELECT ONE OF THE
FOUR REGIONS AND THEN CONDITIONALLY TO GENERATE A VALUE
FROM THE REGION
V:  A U(0,1) RANDOM NUMBER USED TO GENERATE THE RANDOM VALUE
(REGION 1) OR TRANSFORMED INTO THE VARIATE TO ACCEPT OR
REJECT THE CANDIDATE VALUE
IX:  INTEGER CANDIDATE VALUE
X:  PRELIMINARY CONTINUOUS CANDIDATE VALUE IN REGION 2 LOGIC
AND A FLOATING POINT IX IN THE ACCEPT/REJECT LOGIC
K:  ABSOLUTE VALUE OF (IX-M)
F:  THE HEIGHT OF THE SCALED DENSITY FUNCTION USED IN THE
ACCEPT/REJECT DECISION WHEN BOTH M AND IX ARE SMALL
ALSO USED IN THE INVERSE TRANSFORMATION
R: THE RATIO P/Q
G: CONSTANT USED IN CALCULATION OF PROBABILITY
MP:  MODE PLUS ONE, THE LOWER INDEX FOR EXPLICIT CALCULATION
OF F WHEN IX IS GREATER THAN M
IX1:  CANDIDATE VALUE PLUS ONE, THE LOWER INDEX FOR EXPLICIT
CALCULATION OF F WHEN IX IS LESS THAN M
I:  INDEX FOR EXPLICIT CALCULATION OF F FOR BTPE
AMAXP: MAXIMUM ERROR OF THE LOGARITHM OF NORMAL BOUND
YNORM: LOGARITHM OF NORMAL BOUND
ALV:  NATURAL LOGARITHM OF THE ACCEPT/REJECT VARIATE V
X1,F1,Z,W,Z2,X2,F2, AND W2 ARE TEMPORARY VARIABLES TO BE
USED IN THE FINAL ACCEPT/REJECT TEST
QN: PROBABILITY OF NO SUCCESS IN N TRIALS
REMARK
IX AND JX COULD LOGICALLY BE THE SAME VARIABLE, WHICH WOULD
SAVE A MEMORY POSITION AND A LINE OF CODE.  HOWEVER, SOME
COMPILERS (E.G.,CDC MNF) OPTIMIZE BETTER WHEN THE ARGUMENTS
ARE NOT INVOLVED.
ISEED NEEDS TO BE DOUBLE PRECISION IF THE IMSL ROUTINE
GGUBFS IS USED TO GENERATE UNIFORM RANDOM NUMBER, OTHERWISE
TYPE OF ISEED SHOULD BE DICTATED BY THE UNIFORM GENERATOR
**********************************************************************
*****DETERMINE APPROPRIATE ALGORITHM AND WHETHER SETUP IS NECESSARY
*/
/* JJV changed initial values to ridiculous values */
	double psave = -1.0E37;
	int32_t nsave = -214748365;
	int32_t ignbin_mt, i, ix, ix1, k, m, mp, T1;
	double al, alv, amaxp, c, f, f1, f2, ffm, fm, g, p, p1, p2, p3, p4, q, qn, r, u, v, w, w2, x, x1,
		x2, xl, xll, xlr, xm, xnp, xnpq, xr, ynorm, z, z2;

	/*
	*****SETUP, PERFORM ONLY WHEN PARAMETERS CHANGE
	JJV added checks to ensure 0.0 <= PP <= 1.0
	*/
	if (pp < 0.0) ERR_CRITICAL("PP < 0.0 in IGNBIN");
	if (pp > 1.0) ERR_CRITICAL("PP > 1.0 in IGNBIN");
	psave = pp;
	p = std::min(psave, 1.0 - psave);
	q = 1.0 - p;

	/*
	JJV added check to ensure N >= 0
	*/
	if (n < 0) ERR_CRITICAL("N < 0 in IGNBIN");
	xnp = n * p;
	nsave = n;
	if (xnp < 30.0) goto S140;
	ffm = xnp + p;
	m = (int32_t)ffm;
	fm = m;
	xnpq = xnp * q;
	p1 = (int32_t)(2.195 * sqrt(xnpq) - 4.6 * q) + 0.5;
	xm = fm + 0.5;
	xl = xm - p1;
	xr = xm + p1;
	c = 0.134 + 20.5 / (15.3 + fm);
	al = (ffm - xl) / (ffm - xl * p);
	xll = al * (1.0 + 0.5 * al);
	al = (xr - ffm) / (xr * q);
	xlr = al * (1.0 + 0.5 * al);
	p2 = p1 * (1.0 + c + c);
	p3 = p2 + c / xll;
	p4 = p3 + c / xlr;
S30:
	/*
	*****GENERATE VARIATE
	*/
	u = ranf_mt(tn) * p4;
	v = ranf_mt(tn);
	/*
	TRIANGULAR REGION
	*/
	if (u > p1) goto S40;
	ix = (int32_t)(xm - p1 * v + u);
	goto S170;
S40:
	/*
	PARALLELOGRAM REGION
	*/
	if (u > p2) goto S50;
	x = xl + (u - p1) / c;
	v = v * c + 1.0 - std::abs(xm - x) / p1;
	if (v > 1.0 || v <= 0.0) goto S30;
	ix = (int32_t)x;
	goto S70;
S50:
	/*
	LEFT TAIL
	*/
	if (u > p3) goto S60;
	ix = (int32_t)(xl + log(v) / xll);
	if (ix < 0) goto S30;
	v *= ((u - p2) * xll);
	goto S70;
S60:
	/*
	RIGHT TAIL
	*/
	ix = (int32_t)(xr - log(v) / xlr);
	if (ix > n) goto S30;
	v *= ((u - p3) * xlr);
S70:
	/*
	*****DETERMINE APPROPRIATE WAY TO PERFORM ACCEPT/REJECT TEST
	*/
	k = std::abs(ix - m);
	if (k > 20 && k < xnpq / 2 - 1) goto S130;
	/*
	EXPLICIT EVALUATION
	*/
	f = 1.0;
	r = p / q;
	g = (n + 1) * r;
	T1 = m - ix;
	if (T1 < 0) goto S80;
	else if (T1 == 0) goto S120;
	else  goto S100;
S80:
	mp = m + 1;
	for (i = mp; i <= ix; i++) f *= (g / i - r);
	goto S120;
S100:
	ix1 = ix + 1;
	for (i = ix1; i <= m; i++) f /= (g / i - r);
S120:
	if (v <= f) goto S170;
	goto S30;
S130:
	/*
	SQUEEZING USING UPPER AND LOWER BOUNDS ON ALOG(F(X))
	*/
	amaxp = k / xnpq * ((k * (k / 3.0 + 0.625) + 0.1666666666666) / xnpq + 0.5);
	ynorm = -(k * k / (2.0 * xnpq));
	alv = log(v);
	if (alv < ynorm - amaxp) goto S170;
	if (alv > ynorm + amaxp) goto S30;
	/*
	STIRLING'S FORMULA TO MACHINE ACCURACY FOR
	THE FINAL ACCEPTANCE/REJECTION TEST
	*/
	x1 = ix + 1.0;
	f1 = fm + 1.0;
	z = n + 1.0 - fm;
	w = n - ix + 1.0;
	z2 = z * z;
	x2 = x1 * x1;
	f2 = f1 * f1;
	w2 = w * w;
	if (alv <= xm * log(f1 / x1) + (n - m + 0.5) * log(z / w) + (ix - m) * log(w * p / (x1 * q)) + (13860.0 -
		(462.0 - (132.0 - (99.0 - 140.0 / f2) / f2) / f2) / f2) / f1 / 166320.0 + (13860.0 - (462.0 -
		(132.0 - (99.0 - 140.0 / z2) / z2) / z2) / z2) / z / 166320.0 + (13860.0 - (462.0 - (132.0 -
			(99.0 - 140.0 / x2) / x2) / x2) / x2) / x1 / 166320.0 + (13860.0 - (462.0 - (132.0 - (99.0
				- 140.0 / w2) / w2) / w2) / w2) / w / 166320.0) goto S170;
	goto S30;
S140:
	/*
	INVERSE CDF LOGIC FOR MEAN LESS THAN 30
	*/
	qn = pow(q, (double)n);
	r = p / q;
	g = r * (n + 1);
S150:
	ix = 0;
	f = qn;
	u = ranf_mt(tn);
S160:
	if (u < f) goto S170;
	if (ix > 110) goto S150;
	u -= f;
	ix += 1;
	f *= (g / ix - r);
	goto S160;
S170:
	if (psave > 0.5) ix = n - ix;
	ignbin_mt = ix;
	return ignbin_mt;
}
double sexpo(void)
/*
**********************************************************************


	 (STANDARD-)  E X P O N E N T I A L   DISTRIBUTION


**********************************************************************
**********************************************************************

	 FOR DETAILS SEE:

			   AHRENS, J.H. AND DIETER, U.
			   COMPUTER METHODS FOR SAMPLING FROM THE
			   EXPONENTIAL AND NORMAL DISTRIBUTIONS.
			   COMM. ACM, 15,10 (OCT. 1972), 873 - 882.

	 ALL STATEMENT NUMBERS CORRESPOND TO THE STEPS OF ALGORITHM
	 'SA' IN THE ABOVE PAPER (SLIGHTLY MODIFIED IMPLEMENTATION)

	 Modified by Barry W. Brown, Feb 3, 1988 to use RANF instead of
	 SUNIF.  The argument IR thus goes away.

**********************************************************************
	 Q(N) = SUM(ALOG(2.0)**K/K!)    K=1,..,N ,      THE HIGHEST N
	 (HERE 8) IS DETERMINED BY Q(N)=1.0 WITHIN STANDARD PRECISION
*/
{
	static double q[8] = {
		0.6931472,0.9333737,0.9888778,0.9984959,0.9998293,0.9999833,0.9999986,
		.9999999
	};
	int32_t i;
	double sexpo, a, u, ustar, umin;

	a = 0.0;
	u = ranf();
	goto S30;
S20:
	a += q[0];
S30:
	u += u;
	/*
	 * JJV changed the following to reflect the true algorithm and prevent
	 * JJV unpredictable behavior if U is initially 0.5.
	 *  if(u <= 1.0) goto S20;
	 */
	if (u < 1.0) goto S20;
	u -= 1.0;
	if (u > q[0]) goto S60;
	sexpo = a + u;
	return sexpo;
S60:
	i = 1;
	ustar = ranf();
	umin = ustar;
S70:
	ustar = ranf();
	if (ustar < umin) umin = ustar;
	i += 1;
	if (u > q[i - 1]) goto S70;
	return  a + umin * q[0];
}
double sexpo_mt(int tn)
/*
**********************************************************************


(STANDARD-)  E X P O N E N T I A L   DISTRIBUTION


**********************************************************************
**********************************************************************

FOR DETAILS SEE:

AHRENS, J.H. AND DIETER, U.
COMPUTER METHODS FOR SAMPLING FROM THE
EXPONENTIAL AND NORMAL DISTRIBUTIONS.
COMM. ACM, 15,10 (OCT. 1972), 873 - 882.

ALL STATEMENT NUMBERS CORRESPOND TO THE STEPS OF ALGORITHM
'SA' IN THE ABOVE PAPER (SLIGHTLY MODIFIED IMPLEMENTATION)

Modified by Barry W. Brown, Feb 3, 1988 to use RANF instead of
SUNIF.  The argument IR thus goes away.

**********************************************************************
Q(N) = SUM(ALOG(2.0)**K/K!)    K=1,..,N ,      THE HIGHEST N
(HERE 8) IS DETERMINED BY Q(N)=1.0 WITHIN STANDARD PRECISION
*/
{
	double q[8] = {
		0.6931472,0.9333737,0.9888778,0.9984959,0.9998293,0.9999833,0.9999986,
			.9999999
	};
	int32_t i;
	double sexpo_mt, a, u, ustar, umin;

	a = 0.0;
	u = ranf_mt(tn);
	goto S30;
S20:
	a += q[0];
S30:
	u += u;
	/*
	* JJV changed the following to reflect the true algorithm and prevent
	* JJV unpredictable behavior if U is initially 0.5.
	*  if(u <= 1.0) goto S20;
	*/
	if (u < 1.0) goto S20;
	u -= 1.0;
	if (u > q[0]) goto S60;
	sexpo_mt = a + u;
	return sexpo_mt;
S60:
	i = 1;
	ustar = ranf_mt(tn);
	umin = ustar;
S70:
	ustar = ranf_mt(tn);
	if (ustar < umin) umin = ustar;
	i += 1;
	if (u > q[i - 1]) goto S70;
	return  a + umin * q[0];
}
double snorm(void)
/*
**********************************************************************


	 (STANDARD-)  N O R M A L  DISTRIBUTION


**********************************************************************
**********************************************************************

	 FOR DETAILS SEE:

			   AHRENS, J.H. AND DIETER, U.
			   EXTENSIONS OF FORSYTHE'S METHOD FOR RANDOM
			   SAMPLING FROM THE NORMAL DISTRIBUTION.
			   MATH. COMPUT., 27,124 (OCT. 1973), 927 - 937.

	 ALL STATEMENT NUMBERS CORRESPOND TO THE STEPS OF ALGORITHM 'FL'
	 (M=5) IN THE ABOVE PAPER     (SLIGHTLY MODIFIED IMPLEMENTATION)

	 Modified by Barry W. Brown, Feb 3, 1988 to use RANF instead of
	 SUNIF.  The argument IR thus goes away.

**********************************************************************
	 THE DEFINITIONS OF THE CONSTANTS A(K), D(K), T(K) AND
	 H(K) ARE ACCORDING TO THE ABOVEMENTIONED ARTICLE
*/
{
	static double a[32] = {
		0.0,3.917609E-2,7.841241E-2,0.11777,0.1573107,0.1970991,0.2372021,0.2776904,
		0.3186394,0.36013,0.4022501,0.4450965,0.4887764,0.5334097,0.5791322,
		0.626099,0.6744898,0.7245144,0.7764218,0.8305109,0.8871466,0.9467818,
		1.00999,1.077516,1.150349,1.229859,1.318011,1.417797,1.534121,1.67594,
		1.862732,2.153875
	};
	static double d[31] = {
		0.0,0.0,0.0,0.0,0.0,0.2636843,0.2425085,0.2255674,0.2116342,0.1999243,
		0.1899108,0.1812252,0.1736014,0.1668419,0.1607967,0.1553497,0.1504094,
		0.1459026,0.14177,0.1379632,0.1344418,0.1311722,0.128126,0.1252791,
		0.1226109,0.1201036,0.1177417,0.1155119,0.1134023,0.1114027,0.1095039
	};
	static double t[31] = {
		7.673828E-4,2.30687E-3,3.860618E-3,5.438454E-3,7.0507E-3,8.708396E-3,
		1.042357E-2,1.220953E-2,1.408125E-2,1.605579E-2,1.81529E-2,2.039573E-2,
		2.281177E-2,2.543407E-2,2.830296E-2,3.146822E-2,3.499233E-2,3.895483E-2,
		4.345878E-2,4.864035E-2,5.468334E-2,6.184222E-2,7.047983E-2,8.113195E-2,
		9.462444E-2,0.1123001,0.136498,0.1716886,0.2276241,0.330498,0.5847031
	};
	static double h[31] = {
		3.920617E-2,3.932705E-2,3.951E-2,3.975703E-2,4.007093E-2,4.045533E-2,
		4.091481E-2,4.145507E-2,4.208311E-2,4.280748E-2,4.363863E-2,4.458932E-2,
		4.567523E-2,4.691571E-2,4.833487E-2,4.996298E-2,5.183859E-2,5.401138E-2,
		5.654656E-2,5.95313E-2,6.308489E-2,6.737503E-2,7.264544E-2,7.926471E-2,
		8.781922E-2,9.930398E-2,0.11556,0.1404344,0.1836142,0.2790016,0.7010474
	};
	int32_t i; //made this non-static: ggilani 27/11/14
	double snorm, u, s, ustar, aa, w, y, tt; //made this non-static: ggilani 27/11/14
	u = ranf();
	s = 0.0;
	if (u > 0.5) s = 1.0;
	u += (u - s);
	u = 32.0 * u;
	i = (int32_t)(u);
	if (i == 32) i = 31;
	if (i == 0) goto S100;
	/*
									START CENTER
	*/
	ustar = u - (double)i;
	aa = *(a + i - 1);
S40:
	if (ustar <= *(t + i - 1)) goto S60;
	w = (ustar - *(t + i - 1)) * *(h + i - 1);
S50:
	/*
									EXIT   (BOTH CASES)
	*/
	y = aa + w;
	snorm = y;
	if (s == 1.0) snorm = -y;
	return snorm;
S60:
	/*
									CENTER CONTINUED
	*/
	u = ranf();
	w = u * (*(a + i) - aa);
	tt = (0.5 * w + aa) * w;
	goto S80;
S70:
	tt = u;
	ustar = ranf();
S80:
	if (ustar > tt) goto S50;
	u = ranf();
	if (ustar >= u) goto S70;
	ustar = ranf();
	goto S40;
S100:
	/*
									START TAIL
	*/
	i = 6;
	aa = *(a + 31);
	goto S120;
S110:
	aa += *(d + i - 1);
	i += 1;
S120:
	u += u;
	if (u < 1.0) goto S110;
	u -= 1.0;
S140:
	w = u * *(d + i - 1);
	tt = (0.5 * w + aa) * w;
	goto S160;
S150:
	tt = u;
S160:
	ustar = ranf();
	if (ustar > tt) goto S50;
	u = ranf();
	if (ustar >= u) goto S150;
	u = ranf();
	goto S140;
}
double snorm_mt(int tn)
/*
**********************************************************************


(STANDARD-)  N O R M A L  DISTRIBUTION


**********************************************************************
**********************************************************************

FOR DETAILS SEE:

AHRENS, J.H. AND DIETER, U.
EXTENSIONS OF FORSYTHE'S METHOD FOR RANDOM
SAMPLING FROM THE NORMAL DISTRIBUTION.
MATH. COMPUT., 27,124 (OCT. 1973), 927 - 937.

ALL STATEMENT NUMBERS CORRESPOND TO THE STEPS OF ALGORITHM 'FL'
(M=5) IN THE ABOVE PAPER     (SLIGHTLY MODIFIED IMPLEMENTATION)

Modified by Barry W. Brown, Feb 3, 1988 to use RANF instead of
SUNIF.  The argument IR thus goes away.

**********************************************************************
THE DEFINITIONS OF THE CONSTANTS A(K), D(K), T(K) AND
H(K) ARE ACCORDING TO THE ABOVEMENTIONED ARTICLE
*/
{
	static double a[32] = {
		0.0,3.917609E-2,7.841241E-2,0.11777,0.1573107,0.1970991,0.2372021,0.2776904,
			0.3186394,0.36013,0.4022501,0.4450965,0.4887764,0.5334097,0.5791322,
			0.626099,0.6744898,0.7245144,0.7764218,0.8305109,0.8871466,0.9467818,
			1.00999,1.077516,1.150349,1.229859,1.318011,1.417797,1.534121,1.67594,
			1.862732,2.153875
	};
	static double d[31] = {
		0.0,0.0,0.0,0.0,0.0,0.2636843,0.2425085,0.2255674,0.2116342,0.1999243,
			0.1899108,0.1812252,0.1736014,0.1668419,0.1607967,0.1553497,0.1504094,
			0.1459026,0.14177,0.1379632,0.1344418,0.1311722,0.128126,0.1252791,
			0.1226109,0.1201036,0.1177417,0.1155119,0.1134023,0.1114027,0.1095039
	};
	static double t[31] = {
		7.673828E-4,2.30687E-3,3.860618E-3,5.438454E-3,7.0507E-3,8.708396E-3,
			1.042357E-2,1.220953E-2,1.408125E-2,1.605579E-2,1.81529E-2,2.039573E-2,
			2.281177E-2,2.543407E-2,2.830296E-2,3.146822E-2,3.499233E-2,3.895483E-2,
			4.345878E-2,4.864035E-2,5.468334E-2,6.184222E-2,7.047983E-2,8.113195E-2,
			9.462444E-2,0.1123001,0.136498,0.1716886,0.2276241,0.330498,0.5847031
	};
	static double h[31] = {
		3.920617E-2,3.932705E-2,3.951E-2,3.975703E-2,4.007093E-2,4.045533E-2,
			4.091481E-2,4.145507E-2,4.208311E-2,4.280748E-2,4.363863E-2,4.458932E-2,
			4.567523E-2,4.691571E-2,4.833487E-2,4.996298E-2,5.183859E-2,5.401138E-2,
			5.654656E-2,5.95313E-2,6.308489E-2,6.737503E-2,7.264544E-2,7.926471E-2,
			8.781922E-2,9.930398E-2,0.11556,0.1404344,0.1836142,0.2790016,0.7010474
	};
	int32_t i;
	double snorm_mt, u, s, ustar, aa, w, y, tt;
	u = ranf_mt(tn);
	s = 0.0;
	if (u > 0.5) s = 1.0;
	u += (u - s);
	u = 32.0 * u;
	i = (int32_t)(u);
	if (i == 32) i = 31;
	if (i == 0) goto S100;
	/*
	START CENTER
	*/
	ustar = u - (double)i;
	aa = *(a + i - 1);
S40:
	if (ustar <= *(t + i - 1)) goto S60;
	w = (ustar - *(t + i - 1)) * *(h + i - 1);
S50:
	/*
	EXIT   (BOTH CASES)
	*/
	y = aa + w;
	snorm_mt = y;
	if (s == 1.0) snorm_mt = -y;
	return snorm_mt;
S60:
	/*
	CENTER CONTINUED
	*/
	u = ranf_mt(tn);
	w = u * (*(a + i) - aa);
	tt = (0.5 * w + aa) * w;
	goto S80;
S70:
	tt = u;
	ustar = ranf_mt(tn);
S80:
	if (ustar > tt) goto S50;
	u = ranf_mt(tn);
	if (ustar >= u) goto S70;
	ustar = ranf_mt(tn);
	goto S40;
S100:
	/*
	START TAIL
	*/
	i = 6;
	aa = *(a + 31);
	goto S120;
S110:
	aa += *(d + i - 1);
	i += 1;
S120:
	u += u;
	if (u < 1.0) goto S110;
	u -= 1.0;
S140:
	w = u * *(d + i - 1);
	tt = (0.5 * w + aa) * w;
	goto S160;
S150:
	tt = u;
S160:
	ustar = ranf_mt(tn);
	if (ustar > tt) goto S50;
	u = ranf_mt(tn);
	if (ustar >= u) goto S150;
	u = ranf_mt(tn);
	goto S140;
}
double fsign(double num, double sign)
/* Transfers sign of argument sign to argument num */
{
	if ((sign > 0.0f && num < 0.0f) || (sign < 0.0f && num>0.0f))
		return -num;
	else return num;
}
/*function gen_snorm
 * purpose: my own implementation of sampling from a uniform distribution, using Marsaglia polar method, but for multi-threading
 *
 * author: ggilani, date: 28/11/14
 */
double gen_norm_mt(double mu, double sd, int tn)
{
	double u, v, x, S;

	do
	{
		u = 2 * ranf_mt(tn) - 1; //u and v are uniform random numbers on the interval [-1,1]
		v = 2 * ranf_mt(tn) - 1;

		//calculate S=U^2+V^2
		S = u * u + v * v;
	} while (S >= 1 || S == 0);

	//calculate x,y - both of which are normally distributed
	x = u * sqrt((-2 * log(S)) / S);

	// This routine can be accelerated by storing the second normal value
	// and using it for the next call.
	//	y = v * sqrt((-2 * log(S)) / S);

	//return x
	return x * sd + mu;
}
/*function gen_gamma_mt
 * purpose: my own implementation of sampling from a gamma distribution, using Marsaglia-Tsang method, but for multi-threading
 *
 * author: ggilani, date: 01/12/14
 */
double gen_gamma_mt(double beta, double alpha, int tn)
{
	double d, c, u, v, z, f, alpha2, gamma;

	//error statment if either beta or alpha are <=0, as gamma distribution is undefined in this case
	if ((beta <= 0) || (alpha <= 0))
	{
		ERR_CRITICAL("Gamma distribution parameters in undefined range!\n");
	}

	//method is slightly different depending on whether alpha is greater than or equal to 1, or less than 1
	if (alpha >= 1)
	{
		d = alpha - (1.0 / 3.0);
		c = 1.0 / (3.0 * sqrt(d));
		do
		{
			//sample one random number from uniform distribution and one from standard normal distribution
			u = ranf_mt(tn);
			z = gen_norm_mt(0, 1, tn);
			v = 1 + z * c;
			v = v * v * v;
		} while ((z <= (-1.0 / c)) && (log(u) >= (0.5 * z * z + d - d * v + d * log(v))));
		//if beta is equal to 1, there is no scale. If beta is not equal to one, return scaled value
		if (beta != 1)
		{
			return (d * v) / beta;
		}
		else
		{
			return d * v;
		}
	}
	else
	{
		//if alpha is less than 1, initially sample from gamma(beta,alpha+1)
		alpha2 = alpha + 1;
		d = alpha2 - (1.0 / 3.0);
		c = 1.0 / (3.0 * sqrt(d));
		do
		{
			//sample one random number from uniform distribution and one from standard normal distribution
			u = ranf_mt(tn);
			z = gen_norm_mt(0, 1, tn);
			v = 1 + z * c;
			v = v * v * v;
		} while ((z <= (-1.0 / c)) && (log(u) >= (0.5 * z * z + d - d * v + d * log(v))));
		//if beta is equal to 1, there is no scale. If beta is not equal to one, return scaled value
		if (beta != 1)
		{
			gamma = (d * v) / beta;
		}
		else
		{
			gamma = d * v;
		}
		//now rescale again to take into account that alpha is less than 1
		f = pow(ranf_mt(tn), (1.0 / alpha));
		//return gamma scaled by f
		return gamma * f;
	}
}
/* function gen_lognormal(double mu, double sigma)
 * purpose: to generate samples from a lognormal distribution with parameters mu and sigma
 *
 * parameters:
 *  mean mu
 *  standard deviation sigma
 *
 * returns: double from the specified lognormal distribution
 *
 * author: ggilani, date: 09/02/17
 */
double gen_lognormal(double mu, double sigma)
{
	double randnorm, location, scale;

	randnorm = snorm();
	location = log(mu / sqrt(1 + ((sigma * sigma) / (mu * mu))));
	scale = sqrt(log(1 + ((sigma * sigma) / (mu * mu))));
	return exp(location + scale * randnorm);
}
void SampleWithoutReplacement(int tn, int k, int n)
{
	/* Based on algorithm SG of http://portal.acm.org/citation.cfm?id=214402
	ACM Transactions on Mathematical Software (TOMS) archive
	Volume 11 ,  Issue 2  (June 1985) table of contents
	Pages: 157 - 169
	Year of Publication: 1985
	ISSN:0098-3500
	*/

	double t, r, a, mu, f;
	int i, j, q, b;

	if (k < 3)
	{
		for (i = 0; i < k; i++)
		{
			do
			{
				SamplingQueue[tn][i] = (int)(ranf_mt(tn) * ((double)n));
// This original formulation is completely valid, but the PVS Studio analyzer
// notes this, so I am changing it just to get report-clean.
// "V1008 Consider inspecting the 'for' operator. No more than one iteration of the loop will be performed. Rand.cpp 2450"
//				for (j = q = 0; (j < i) && (!q); j++)
//					q = (SamplingQueue[tn][i] == SamplingQueue[tn][j]);
				j = q = 0;
				if (i == 1)
					q = (SamplingQueue[tn][i] == SamplingQueue[tn][j]);
			} while (q);
		}
		q = k;
	}
	else if (2 * k > n)
	{
		for (i = 0; i < n; i++)
			SamplingQueue[tn][i] = i;
		for (i = n; i > k; i--)
		{
			j = (int)(ranf_mt(tn) * ((double)i));
			if (j != i - 1)
			{
				b = SamplingQueue[tn][j];
				SamplingQueue[tn][j] = SamplingQueue[tn][i - 1];
				SamplingQueue[tn][i - 1] = b;
			}
		}
		q = k;
	}
	else if (4 * k > n)
	{
		for (i = 0; i < n; i++)
			SamplingQueue[tn][i] = i;
		for (i = 0; i < k; i++)
		{
			j = (int)(ranf_mt(tn) * ((double)(n - i)));
			if (j > 0)
			{
				b = SamplingQueue[tn][i];
				SamplingQueue[tn][i] = SamplingQueue[tn][i + j];
				SamplingQueue[tn][i + j] = b;
			}
		}
		q = k;
	}
	else
	{
		/* fprintf(stderr,"@%i %i:",k,n); */
		t = (double)k;
		r = sqrt(t);
		a = sqrt(log(1 + t / 2 * PI));
		a = a + a * a / (3 * r);
		mu = t + a * r;
		b = 2 * MAX_PLACE_SIZE; /* (int) (k+4*a*r); */
		f = -1 / (log(1 - mu / ((double)n)));
		i = q = 0;
		while (i <= n)
		{
			i += (int)ceil(-log(ranf_mt(tn)) * f);
			if (i <= n)
			{
				SamplingQueue[tn][q] = i - 1;
				q++;
				if (q >= b) i = q = 0;
			}
			else if (q < k)
				i = q = 0;
		}
	}
	/*	else
			{
			t=(double) (n-k);
			r=sqrt(t);
			a=sqrt(log(1+t/2*PI));
			a=a+a*a/(3*r);
			mu=t+a*r;
			b=2*MAX_PLACE_SIZE;
			f=-1/(log(1-mu/((double) n)));
			i=q=0;
			while(i<=n)
				{
				int i2=i+(int) ceil(-log(ranf_mt(tn))*f);
				i++;
				if(i2<=n)
					for(;(i<i2)&&(q<b);i++)
						{
						SamplingQueue[tn][q]=i-1;
						q++;
						}
				else
					{
					for(;(i<=n)&&(q<b);i++)
						{
						SamplingQueue[tn][q]=i-1;
						q++;
						}
					if(q<k) i=q=0;
					}
				if(q>=b) i=q=0;
				}
			}
	*/
	/*	if(k>2)
			{
			fprintf(stderr,"(%i) ",q);
			for(i=0;i<q;i++) fprintf(stderr,"%i ",SamplingQueue[tn][i]);
			fprintf(stderr,"\n");
			}
	*/	while (q > k)
	{
		i = (int)(ranf_mt(tn) * ((double)q));
		if (i < q - 1) SamplingQueue[tn][i] = SamplingQueue[tn][q - 1];
		q--;
	}

}

