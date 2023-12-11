/* Copyright (c) 2012 Massachusetts Institute of Technology
 * 
 * [Also included are functions derived from derfc in SLATEC
 *  (netlib.org/slatec), which "is in the public domain"
 *  and hence may be redistributed under these or any terms.]
 *
 * Permission is hereby granted, free of charge, to any person obtaining
 * a copy of this software and associated documentation files (the
 * "Software"), to deal in the Software without restriction, including
 * without limitation the rights to use, copy, modify, merge, publish,
 * distribute, sublicense, and/or sell copies of the Software, and to
 * permit persons to whom the Software is furnished to do so, subject to
 * the following conditions:
 * 
 * The above copyright notice and this permission notice shall be
 * included in all copies or substantial portions of the Software.
 * 
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND,
 * EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF
 * MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND
 * NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE
 * LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION
 * OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION
 * WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE. 
 */

// include this declaration in your code or header file to call:
#include <complex.h>
//extern std::double complex Faddeeva_w(std::double complex z,double relerr=0);

/* Available at: http://ab-initio.mit.edu/Faddeeva_w
 
 Compute the Faddeyeva function w(z) = exp(-z^2) * erfc(-i*z),
 to a desired relative accuracy relerr, for arbitrary complex
 arguments z.
 
 For sufficiently large |z|, we use a continued-fraction expansion
 for w(z) similar to those described in:
 
 Walter Gautschi, "Efficient computation of the complex error
 function," SIAM J. Numer. Anal. 7(1), pp. 187-198 (1970)
 
 G. P. M. Poppe and C. M. J. Wijers, "More efficient computation
 of the complex error function," ACM Trans. Math. Soft. 16(1),
 pp. 38-46 (1990).
 
 Unlike those papers, however, we switch to a completely different
 algorithm for smaller |z|:
 
 Mofreh R. Zaghloul and Ahmed N. Ali, "Algorithm 916: Computing the
 Faddeyeva and Voigt Functions," ACM Trans. Math. Soft. 38(2), 15
 (2011).
 
 (I initially used this algorithm for all z, but it turned out to be
 significantly slower than the continued-fraction expansion for
 larger |z|.  On the other hand, it is competitive for smaller |z|, 
 and is significantly more accurate than the Poppe & Wijers code
 in some regions, e.g. in the vicinity of z=1+1i.)
 
 Note that this is an INDEPENDENT RE-IMPLEMENTATION of these algorithms,
 based on the description in the papers ONLY.  In particular, I did
 not refer to the authors' Fortran or Matlab implementations, respectively,
 (which are under restrictive ACM copyright terms and therefore unusable
 in free/open-source software).
 
 Steven G. Johnson, Massachusetts Institute of Technology
 http://math.mit.edu/~stevenj
 October 2012.
 
 -- Note that Algorithm 916 assumes that the erfc(x) function, 
 or rather the scaled function erfcx(x) = exp(x*x)*erfc(x),
 is supplied for REAL arguments x.   I include an erfcx routine
 derived from the derfc function in SLATEC (modified to
 compute erfcx instead of erfc).
 
 A small test program is included the end, which checks
 the w(z) results against several known values.  To compile
 the test function, compile with -DFADDEEVA_W_TEST (that is,
 #define FADDEEVA_W_TEST).
 
 REVISION HISTORY:
 4 October 2012: Initial public release (SGJ)
 5 October 2012: Revised (SGJ) to fix spelling error,
 start summation for large x at round(x/a) (> 1)
 rather than ceil(x/a) as in the original
 paper, which should slightly improve performance
 (and, apparently, slightly improves accuracy)
 19 October 2012: Revised (SGJ) to fix bugs for large x, large -y,
 and 15<x<26. Performance improvements. Prototype
 now supplies default value for relerr.
 24 October 2012: Switch to continued-fraction expansion for
 sufficiently large z, for performance reasons.
 Also, avoid spurious overflow for |z| > 1e154.
 Set relerr argument to min(relerr,0.1).
 */

#include <float.h>
#include <math.h>

//using namespace std;

static double derfcx(double x);
#define erfcx(x) derfcx(x)

// return sinc(x) = sin(x)/x, given both x and sin(x) 
// [since we only use this in cases where sin(x) has already been computed]
static inline double sinc(double x, double sinx) { 
	return fabs(x) < 1e-4 ? 1 - (0.1666666666666666666667)*x*x : sinx / x; 
}

// sinh(x) via Taylor series, accurate to machine precision for |x| < 1e-2
static inline double sinh_taylor(double x) {
	return x * (1 + (x*x) * (0.1666666666666666666667
							 + 0.00833333333333333333333 * (x*x)));
}

static inline double sqr(double x) { return x*x; }

// precomputed table of expa2n2[n-1] = exp(-a2*n*n)
// for double-precision a2 = 0.26865... in Faddeeva_w, below.
static const double expa2n2[] = {
7.64405281671221563e-01,
3.41424527166548425e-01,
8.91072646929412548e-02,
1.35887299055460086e-02,
1.21085455253437481e-03,
6.30452613933449404e-05,
1.91805156577114683e-06,
3.40969447714832381e-08,
3.54175089099469393e-10,
2.14965079583260682e-12,
7.62368911833724354e-15,
1.57982797110681093e-17,
1.91294189103582677e-20,
1.35344656764205340e-23,
5.59535712428588720e-27,
1.35164257972401769e-30,
1.90784582843501167e-34,
1.57351920291442930e-38,
7.58312432328032845e-43,
2.13536275438697082e-47,
3.51352063787195769e-52,
3.37800830266396920e-57,
1.89769439468301000e-62,
6.22929926072668851e-68,
1.19481172006938722e-73,
1.33908181133005953e-79,
8.76924303483223939e-86,
3.35555576166254986e-92,
7.50264110688173024e-99,
9.80192200745410268e-106,
7.48265412822268959e-113,
3.33770122566809425e-120,
8.69934598159861140e-128,
1.32486951484088852e-135,
1.17898144201315253e-143,
6.13039120236180012e-152,
1.86258785950822098e-160,
3.30668408201432783e-169,
3.43017280887946235e-178,
2.07915397775808219e-187,
7.36384545323984966e-197,
1.52394760394085741e-206,
1.84281935046532100e-216,
1.30209553802992923e-226,
5.37588903521080531e-237,
1.29689584599763145e-247,
1.82813078022866562e-258,
1.50576355348684241e-269,
7.24692320799294194e-281,
2.03797051314726829e-292,
3.34880215927873807e-304,
0.0 // underflow (also prevents reads past array end, below)
};

double complex Faddeeva_w(double complex z, double relerr)
{
	double a, a2, c;
	if (relerr <= DBL_EPSILON) {
		relerr = DBL_EPSILON;
		a = 0.518321480430085929872; // pi / sqrt(-log(eps*0.5))
		c = 0.329973702884629072537; // (2/pi) * a;
		a2 = 0.268657157075235951582; // a^2
	}
	else {
		const double pi = 3.14159265358979323846264338327950288419716939937510582;
		if (relerr > 0.1) relerr = 0.1; // not sensible to compute < 1 digit
		a = pi / sqrt(-log(relerr*0.5));
		c = (2/pi)*a;
		a2 = a*a;
	}
	const double x = fabs(creal(z));
	const double y = cimag(z);
	
	double complex ret=0.+0.*I; // return value
	
	double sum1 = 0, sum2 = 0, sum3 = 0, sum4 = 0, sum5 = 0;
	
#define USE_CONTINUED_FRACTION 1 // 1 to use continued fraction for large |z|
	
#if USE_CONTINUED_FRACTION
	double ya;
	if ((ya = fabs(y)) > 15 || x > 6) { // continued fraction is faster
		/* Poppe & Wijers suggest using a number of terms
		 nu = 3 + 1442 / (26*rho + 77)
		 where rho = sqrt((x/x0)^2 + (y/y0)^2) where x0=6.3, y0=4.4.
		 (They only use this expansion for rho >= 1, but rho a little less
		 than 1 seems okay too.)
		 Instead, I did my own fit to a slightly different function
		 that avoids the hypotenuse calculation, using NLopt to minimize
		 the sum of the squares of the errors in nu with the constraint
		 that the estimated nu be >= minimum nu to attain machine precision.
		 I also separate the regions where nu == 2 and nu == 1. */
		const double ispi = 0.56418958354775628694807945156; // 1 / sqrt(pi)
		double xs = y < 0 ? -creal(z) : creal(z); // compute for -z if y < 0
		if (x + ya > 4000) { // nu <= 2
			if (x + ya > 1e7) { // nu == 1, w(z) = i/sqrt(pi) / z
				// scale to avoid overflow
				if (x > ya) {
					double yax = ya / xs; 
					double denom = ispi / (xs + yax*ya);
					ret = denom*yax + denom*I;
				}
				else {
					double xya = xs / ya;
					double denom = ispi / (xya*xs + ya);
					ret = (denom + denom*xya*I);
				}
			}
			else { // nu == 2, w(z) = i/sqrt(pi) * z / (z*z - 0.5)
				double dr = xs*xs - ya*ya - 0.5, di = 2*xs*ya;
				double denom = ispi / (dr*dr + di*di);
				ret = (denom * (xs*di-ya*dr) + denom * (xs*dr+ya*di)*I);
			}
		}
		else { // compute nu(z) estimate and do general continued fraction
			const double c0=3.9, c1=11.398, c2=0.08254, c3=0.1421, c4=0.2023; // fit
			double nu = floor(c0 + c1 / (c2*x + c3*ya + c4));
			double wr = xs, wi = ya;
			for (nu = 0.5 * (nu - 1); nu > 0.4; nu -= 0.5) {
				// w <- z - nu/w:
				double denom = nu / (wr*wr + wi*wi);
				wr = xs - wr * denom;
				wi = ya + wi * denom;
			}
			{ // w(z) = i/sqrt(pi) / w:
				double denom = ispi / (wr*wr + wi*wi);
				ret = (denom*wi + denom*wr*I);
			}
		}
		if (y < 0) {
			// use w(z) = 2.0*exp(-z*z) - w(-z), 
			// but be careful of overflow in exp(-z*z) 
			//                                = exp(-(xs*xs-ya*ya) -2*i*xs*ya) 
			return 2.0*cexp(((ya-xs)*(xs+ya) + 2*xs*y*I)) - ret;
		}
		else
			return ret;
	}
	else
#endif // USE_CONTINUED_FRACTION
		if (x == 0.0)
			return (erfcx(y) + x*I); // give correct sign of 0 in imag(w)
	
	/* Note: The test that seems to be suggested in the paper is x <
     sqrt(-log(DBL_MIN)), about 26.6, since otherwise exp(-x^2)
     underflows to zero and sum1,sum2,sum4 are zero.  However, long
     before this occurs, the sum1,sum2,sum4 contributions are
     negligible in double precision; I find that this happens for x >
     about 6, for all y.  On the other hand, I find that the case
     where we compute all of the sums is faster (at least with the
     precomputed expa2n2 table) until about x=10.  Furthermore, if we
     try to compute all of the sums for x > 20, I find that we
     sometimes run into numerical problems because underflow/overflow
     problems start to appear in the various coefficients of the sums,
     below.  Therefore, we use x < 10 here. */
#if USE_CONTINUED_FRACTION
		else { // x < 10 always here
#else
			else if (x < 10) {
#endif
				double prod2ax = 1, prodm2ax = 1;
				double expx2;
				int n;
				
				/* Somewhat ugly copy-and-paste duplication here, but I see significant
				 speedups from using the special-case code with the precomputed
				 exponential, and the x < 5e-4 special case is needed for accuracy. */
				
				if (relerr == DBL_EPSILON) { // use precomputed exp(-a2*(n*n)) table
					if (x < 5e-4) { // compute sum4 and sum5 together as sum5-sum4
						const double x2 = x*x;
						expx2 = 1 - x2 * (1 - 0.5*x2); // exp(-x*x) via Taylor
						// compute exp(2*a*x) and exp(-2*a*x) via Taylor, to double precision
						const double ax2 = 1.036642960860171859744*x; // 2*a*x
						const double exp2ax =
						1 + ax2 * (1 + ax2 * (0.5 + 0.166666666666666666667*ax2));
						const double expm2ax =
						1 - ax2 * (1 - ax2 * (0.5 - 0.166666666666666666667*ax2));
						for (n = 1; 1; ++n) {
							const double coef = expa2n2[n-1] * expx2 / (a2*(n*n) + y*y);
							prod2ax *= exp2ax;
							prodm2ax *= expm2ax;
							sum1 += coef;
							sum2 += coef * prodm2ax;
							sum3 += coef * prod2ax;
							
							// really = sum5 - sum4
							sum5 += coef * (2*a) * n * sinh_taylor((2*a)*n*x);
							
							// test convergence via sum3
							if (coef * prod2ax < relerr * sum3) break;
						}
					}
					else { // x > 5e-4, compute sum4 and sum5 separately
						expx2 = exp(-x*x);
						const double exp2ax = exp((2*a)*x), expm2ax = 1 / exp2ax;
						for (n = 1; 1; ++n) {
							const double coef = expa2n2[n-1] * expx2 / (a2*(n*n) + y*y);
							prod2ax *= exp2ax;
							prodm2ax *= expm2ax;
							sum1 += coef;
							sum2 += coef * prodm2ax;
							sum4 += (coef * prodm2ax) * (a*n);
							sum3 += coef * prod2ax;
							sum5 += (coef * prod2ax) * (a*n);
							// test convergence via sum5, since this sum has the slowest decay
							if ((coef * prod2ax) * (a*n) < relerr * sum5) break;
						}
					}
				}
				else { // relerr != DBL_EPSILON, compute exp(-a2*(n*n)) on the fly
					const double exp2ax = exp((2*a)*x), expm2ax = 1 / exp2ax;
					if (x < 5e-4) { // compute sum4 and sum5 together as sum5-sum4
						const double x2 = x*x;
						expx2 = 1 - x2 * (1 - 0.5*x2); // exp(-x*x) via Taylor
						for (n = 1; 1; ++n) {
							const double coef = exp(-a2*(n*n)) * expx2 / (a2*(n*n) + y*y);
							prod2ax *= exp2ax;
							prodm2ax *= expm2ax;
							sum1 += coef;
							sum2 += coef * prodm2ax;
							sum3 += coef * prod2ax;
							
							// really = sum5 - sum4
							sum5 += coef * (2*a) * n * sinh_taylor((2*a)*n*x);
							
							// test convergence via sum3
							if (coef * prod2ax < relerr * sum3) break;
						}
					}
					else { // x > 5e-4, compute sum4 and sum5 separately
						expx2 = exp(-x*x);
						for (n = 1; 1; ++n) {
							const double coef = exp(-a2*(n*n)) * expx2 / (a2*(n*n) + y*y);
							prod2ax *= exp2ax;
							prodm2ax *= expm2ax;
							sum1 += coef;
							sum2 += coef * prodm2ax;
							sum4 += (coef * prodm2ax) * (a*n);
							sum3 += coef * prod2ax;
							sum5 += (coef * prod2ax) * (a*n);
							// test convergence via sum5, since this sum has the slowest decay
							if ((coef * prod2ax) * (a*n) < relerr * sum5) break;
						}
					}
				}
				const double expx2erfcxy = // avoid spurious overflow for large negative y
				y > -6 // for y < -6, erfcx(y) = 2*exp(y*y) to double precision
				? expx2*erfcx(y) : 2*exp(y*y-x*x);
				if (y > 5) { // imaginary terms cancel
					const double sinxy = sin(x*y);
					ret = (expx2erfcxy - c*y*sum1) * cos(2*x*y)
					+ (c*x*expx2) * sinxy * sinc(x*y, sinxy);
				}
				else {
					double xs = creal(z);
					const double sinxy = sin(xs*y);
					const double sin2xy = sin(2*xs*y), cos2xy = cos(2*xs*y);
					const double coef1 = expx2erfcxy - c*y*sum1;
					const double coef2 = c*xs*expx2;
					ret = (coef1 * cos2xy + coef2 * sinxy * sinc(xs*y, sinxy) +
										  (coef2 * sinc(2*xs*y, sin2xy) - coef1 * sin2xy)*I);
				}
			}
#if !USE_CONTINUED_FRACTION
			// Originally used algorithm 916 for large x, but this is
			// unnecessary if the continued-fraction expansion is used for large |z|.
			else { // x large: only sum3 & sum5 contribute (see above note)
				if (y < 0) {
					/* erfcx(y) ~ 2*exp(y*y) + (< 1) if y < 0, so
					 erfcx(y)*exp(-x*x) ~ 2*exp(y*y-x*x) term may not be negligible
					 if y*y - x*x > -36 or so.  So, compute this term just in case. */
					ret = polar(2*exp(y*y-x*x), -2*real(z)*y);
				}
				// (round instead of ceil as in original paper; note that x/a > 1 here)
				double n0 = floor(x/a + 0.5); // sum in both directions, starting at n0
				double dx = a*n0 - x;
				sum3 = exp(-dx*dx) / (a2*(n0*n0) + y*y);
				sum5 = a*n0 * sum3;
				double exp1 = exp(4*a*dx), exp1dn = 1;
				int dn;
				for (dn = 1; n0 - dn > 0; ++dn) { // loop over n0-dn and n0+dn terms
					double np = n0 + dn, nm = n0 - dn;
					double tp = exp(-sqr(a*dn+dx));
					double tm = tp * (exp1dn *= exp1); // trick to get tm from tp
					tp /= (a2*(np*np) + y*y);
					tm /= (a2*(nm*nm) + y*y);
					sum3 += tp + tm;
					sum5 += a * (np * tp + nm * tm);
					if (a * (np * tp + nm * tm) < relerr * sum5) goto finish;
				}
				while (1) { // loop over n0+dn terms only (since n0-dn <= 0)
					double np = n0 + dn++;
					double tp = exp(-sqr(a*dn+dx)) / (a2*(np*np) + y*y);
					sum3 += tp;
					sum5 += a * np * tp;
					if (a * np * tp < relerr * sum5) goto finish;
				}
			}
		finish:
#endif // !USE_CONTINUED_FRACTION
			
			return ret + ((0.5*c)*y*(sum2+sum3) + (0.5*c)*copysign(sum5-sum4, creal(z))*I);
		}
	
	/////////////////////////////////////////////////////////////////////////
	
	/* erfcx function derived from derfc function in SLATEC, converted
	 to C by Steven G. Johnson (with help from f2c) and modified to
	 compute erfcx instead of erfc, October 2012.
	 
	 According to http://www.netlib.org/slatec/guide:
	 
	 The Library is in the public domain and distributed by the
	 Energy Science and Technology Software Center.
	 
	 Energy Science and Technology Software Center
	 P.O. Box 1020
	 Oak Ridge, TN  37831
	 
	 Telephone  615-576-2606
	 E-mail  estsc%a1.adonis.mrouter@zeus.osti.gov
	 */
	
#if 0 /* comment out for hard-coded IEEE values in derfc */
	/* BEGIN PROLOGUE  INITDS */
	/* PURPOSE  Determine the number of terms needed in an orthogonal */
	/*          polynomial series so that it meets a specified accuracy. */
	/* LIBRARY   SLATEC (FNLIB) */
	/* CATEGORY  C3A2 */
	/* TYPE      DOUBLE PRECISION (INITS-S, INITDS-D) */
	/* KEYWORDS  CHEBYSHEV, FNLIB, INITIALIZE, ORTHOGONAL POLYNOMIAL, */
	/*             ORTHOGONAL SERIES, SPECIAL FUNCTIONS */
	/* AUTHOR  Fullerton, W., (LANL) */
	static int initds(const double *os, int nos, double eta)
	{
		/* System generated locals */
		int i__1;
		
		/* Local variables */
		int i__, ii;
		double err;
		
		/* DESCRIPTION */
		/*  Initialize the orthogonal series, represented by the array OS, so */
		/*  that INITDS is the number of terms needed to insure the error is no */
		/*  larger than ETA.  Ordinarily, ETA will be chosen to be one-tenth */
		/*  machine precision. */
		
		/*             Input Arguments -- */
		/*   OS     double precision array of NOS coefficients in an orthogonal */
		/*          series. */
		/*   NOS    number of coefficients in OS. */
		/*   ETA    single precision scalar containing requested accuracy of */
		/*          series. */
		
		/* REVISION HISTORY  (YYMMDD) */
		/*   770601  DATE WRITTEN */
		/*   890531  Changed all specific intrinsics to generic.  (WRB) */
		/*   890831  Modified array declarations.  (WRB) */
		/*   891115  Modified error message.  (WRB) */
		/*   891115  REVISION DATE from Version 3.2 */
		/*   891214  Prologue converted to Version 4.0 format.  (BAB) */
		/*   900315  CALLs to XERROR changed to CALLs to XERMSG.  (THJ) */
		/*   121094  f2c cleanup, eta changed to double(SGJ) */
		/* END PROLOGUE  INITDS */
		/* FIRST EXECUTABLE STATEMENT  INITDS */
		/* Parameter adjustments */
		--os;
		
		/* Function Body */
		err = 0.;
		i__1 = nos;
		for (ii = 1; ii <= i__1; ++ii) {
			i__ = nos + 1 - ii;
			err += fabs(os[i__]);
			if (err > eta) {
				goto L20;
			}
			/* L10: */
		}
		
	L20:
		return i__;
	} /* initds */
#endif
	
	/* BEGIN PROLOGUE  DCSEVL */
	/* PURPOSE  Evaluate a Chebyshev series. */
	/* LIBRARY   SLATEC (FNLIB) */
	/* CATEGORY  C3A2 */
	/* TYPE      DOUBLE PRECISION (CSEVL-S, DCSEVL-D) */
	/* KEYWORDS  CHEBYSHEV SERIES, FNLIB, SPECIAL FUNCTIONS */
	/* AUTHOR  Fullerton, W., (LANL) */
	static double dcsevl_(double x, const double *cs, int n)
	{
		/* System generated locals */
		int i__1;
		
		/* Local variables */
		int i__;
		double b0, b1, b2;
		int ni;
		double twox;
		
		/* DESCRIPTION */
		/*  Evaluate the N-term Chebyshev series CS at X.  Adapted from */
		/*  a method presented in the paper by Broucke referenced below. */
		
		/*       Input Arguments -- */
		/*  X    value at which the series is to be evaluated. */
		/*  CS   array of N terms of a Chebyshev series.  In evaluating */
		/*       CS, only half the first coefficient is summed. */
		/*  N    number of terms in array CS. */
		
		/* REFERENCES  R. Broucke, Ten subroutines for the manipulation of */
		/*                 Chebyshev series, Algorithm 446, Communications of */
		/*                 the A.C.M. 16, (1973) pp. 254-256. */
		/*               L. Fox and I. B. Parker, Chebyshev Polynomials in */
		/*                 Numerical Analysis, Oxford University Press, 1968, */
		/*                 page 56. */
		/* ROUTINES CALLED  D1MACH, XERMSG */
		/* REVISION HISTORY  (YYMMDD) */
		/*   770401  DATE WRITTEN */
		/*   890831  Modified array declarations.  (WRB) */
		/*   890831  REVISION DATE from Version 3.2 */
		/*   891214  Prologue converted to Version 4.0 format.  (BAB) */
		/*   900315  CALLs to XERROR changed to CALLs to XERMSG.  (THJ) */
		/*   900329  Prologued revised extensively and code rewritten to allow */
		/*           X to be slightly outside interval (-1,+1).  (WRB) */
		/*   920501  Reformatted the REFERENCES section.  (WRB) */
		/*   121094  f2c cleanup, made thread-safe (SGJ) */
		/* END PROLOGUE  DCSEVL */
		/* Parameter adjustments */
		--cs;
		
		/* Function Body */
		b1 = 0.;
		b0 = 0.;
		twox = x * 2.;
		i__1 = n;
		for (i__ = 1; i__ <= i__1; ++i__) {
			b2 = b1;
			b1 = b0;
			ni = n + 1 - i__;
			b0 = twox * b1 - b2 + cs[ni];
			/* L10: */
		}
		
		return (b0 - b2) * .5;
	} /* dcsevl_ */
	
	/* BEGIN PROLOGUE  DERFCX (originally DERFC) */
	/* PURPOSE  Compute the scaled complementary error function. */
	/* LIBRARY   SLATEC (FNLIB) */
	/* CATEGORY  C8A, L5A1E */
	/* TYPE      DOUBLE PRECISION (ERFC-S, DERFC-D) */
	/* KEYWORDS  COMPLEMENTARY ERROR FUNCTION, ERFC, FNLIB, */
	/*             SPECIAL FUNCTIONS */
	/* AUTHOR  Fullerton, W., (LANL) */
	static double derfcx(double x)
	{
		/* Initialized data */
		
		const double erfcs[21] = { -.049046121234691808039984544033376,
			-.14226120510371364237824741899631,
			.010035582187599795575754676712933,
			-5.7687646997674847650827025509167e-4,
			2.7419931252196061034422160791471e-5,
			-1.1043175507344507604135381295905e-6,
			3.8488755420345036949961311498174e-8,
			-1.1808582533875466969631751801581e-9,
			3.2334215826050909646402930953354e-11,
			-7.9910159470045487581607374708595e-13,
			1.7990725113961455611967245486634e-14,
			-3.7186354878186926382316828209493e-16,
			7.1035990037142529711689908394666e-18,
			-1.2612455119155225832495424853333e-19,
			2.0916406941769294369170500266666e-21,
			-3.253973102931407298236416e-23,
			4.7668672097976748332373333333333e-25,
			-6.5980120782851343155199999999999e-27,
			8.6550114699637626197333333333333e-29,
			-1.0788925177498064213333333333333e-30,
	    1.2811883993017002666666666666666e-32 };
		const double erc2cs[49] = { -.06960134660230950112739150826197,
			-.04110133936262089348982212084666,
			.003914495866689626881561143705244,
			-4.906395650548979161280935450774e-4,
			7.157479001377036380760894141825e-5,
			-1.153071634131232833808232847912e-5,
			1.994670590201997635052314867709e-6,
			-3.642666471599222873936118430711e-7,
			6.944372610005012589931277214633e-8,
			-1.37122090210436601953460514121e-8,
			2.788389661007137131963860348087e-9,
			-5.814164724331161551864791050316e-10,
			1.23892049175275318118016881795e-10,
			-2.690639145306743432390424937889e-11,
			5.94261435084791098244470968384e-12,
			-1.33238673575811957928775442057e-12,
			3.028046806177132017173697243304e-13,
			-6.966648814941032588795867588954e-14,
			1.620854541053922969812893227628e-14,
			-3.809934465250491999876913057729e-15,
			9.040487815978831149368971012975e-16,
			-2.164006195089607347809812047003e-16,
			5.222102233995854984607980244172e-17,
			-1.26972960236455533637241552778e-17,
			3.109145504276197583836227412951e-18,
			-7.663762920320385524009566714811e-19,
			1.90081925136274520253692973329e-19,
			-4.742207279069039545225655999965e-20,
			1.189649200076528382880683078451e-20,
			-3.000035590325780256845271313066e-21,
			7.602993453043246173019385277098e-22,
			-1.93590944760687288156981104913e-22,
			4.951399124773337881000042386773e-23,
			-1.271807481336371879608621989888e-23,
			3.280049600469513043315841652053e-24,
			-8.492320176822896568924792422399e-25,
			2.206917892807560223519879987199e-25,
			-5.755617245696528498312819507199e-26,
			1.506191533639234250354144051199e-26,
			-3.954502959018796953104285695999e-27,
			1.041529704151500979984645051733e-27,
			-2.751487795278765079450178901333e-28,
			7.29005820549755740899770368e-29,
			-1.936939645915947804077501098666e-29,
			5.160357112051487298370054826666e-30,
			-1.3784193221930940993896448e-30,
			3.691326793107069042251093333333e-31,
			-9.909389590624365420653226666666e-32,
	    2.666491705195388413323946666666e-32 };
		const double erfccs[59] = { .0715179310202924774503697709496,
			-.0265324343376067157558893386681,
			.00171115397792085588332699194606,
			-1.63751663458517884163746404749e-4,
			1.98712935005520364995974806758e-5,
			-2.84371241276655508750175183152e-6,
			4.60616130896313036969379968464e-7,
			-8.22775302587920842057766536366e-8,
			1.59214187277090112989358340826e-8,
			-3.29507136225284321486631665072e-9,
			7.2234397604005554658126115389e-10,
			-1.66485581339872959344695966886e-10,
			4.01039258823766482077671768814e-11,
			-1.00481621442573113272170176283e-11,
			2.60827591330033380859341009439e-12,
			-6.99111056040402486557697812476e-13,
			1.92949233326170708624205749803e-13,
			-5.47013118875433106490125085271e-14,
			1.58966330976269744839084032762e-14,
			-4.7268939801975548392036958429e-15,
			1.4358733767849847867287399784e-15,
			-4.44951056181735839417250062829e-16,
			1.40481088476823343737305537466e-16,
			-4.51381838776421089625963281623e-17,
			1.47452154104513307787018713262e-17,
			-4.89262140694577615436841552532e-18,
			1.64761214141064673895301522827e-18,
			-5.62681717632940809299928521323e-19,
			1.94744338223207851429197867821e-19,
			-6.82630564294842072956664144723e-20,
			2.42198888729864924018301125438e-20,
			-8.69341413350307042563800861857e-21,
			3.15518034622808557122363401262e-21,
			-1.15737232404960874261239486742e-21,
			4.28894716160565394623737097442e-22,
			-1.60503074205761685005737770964e-22,
			6.06329875745380264495069923027e-23,
			-2.31140425169795849098840801367e-23,
			8.88877854066188552554702955697e-24,
			-3.44726057665137652230718495566e-24,
			1.34786546020696506827582774181e-24,
			-5.31179407112502173645873201807e-25,
			2.10934105861978316828954734537e-25,
			-8.43836558792378911598133256738e-26,
			3.39998252494520890627359576337e-26,
			-1.3794523880732420900223837711e-26,
			5.63449031183325261513392634811e-27,
			-2.316490434477065448234277527e-27,
			9.58446284460181015263158381226e-28,
			-3.99072288033010972624224850193e-28,
			1.67212922594447736017228709669e-28,
			-7.04599152276601385638803782587e-29,
			2.97976840286420635412357989444e-29,
			-1.26252246646061929722422632994e-29,
			5.39543870454248793985299653154e-30,
			-2.38099288253145918675346190062e-30,
			1.0990528301027615735972668375e-30,
			-4.86771374164496572732518677435e-31,
	    1.52587726411035756763200828211e-31 };
		const double sqrtpi = 1.77245385090551602729816748334115;
		
		/* System generated locals */
		double d__1;
		
		/* Local variables */
		double y;
		
		/* DERFCX(X) calculates the double precision scaled complementary
		 * error function [DERFCX(X) = EXP(X^2) DERFC(X)] for double precision
		 * argument X.  Original SLATEC function DERFC modified by SGJ in 2012
		 * to compute DERFCX. */
		
		/* Series for ERF        on the interval  0.          to  1.00000E+00 */
		/*                                        with weighted Error   1.28E-32 */
		/*                                         log weighted Error  31.89 */
		/*                               significant figures required  31.05 */
		/*                                    decimal places required  32.55 */
		
		/* Series for ERC2       on the interval  2.50000E-01 to  1.00000E+00 */
		/*                                        with weighted Error   2.67E-32 */
		/*                                         log weighted Error  31.57 */
		/*                               significant figures required  30.31 */
		/*                                    decimal places required  32.42 */
		
		/* Series for ERFC       on the interval  0.          to  2.50000E-01 */
		/*                                        with weighted error   1.53E-31 */
		/*                                         log weighted error  30.82 */
		/*                               significant figures required  29.47 */
		/*                                    decimal places required  31.70 */
		
		/* REFERENCES  (NONE) */
		/* ROUTINES CALLED  D1MACH, DCSEVL, INITDS, XERMSG */
		/* REVISION HISTORY  (YYMMDD) */
		/*   770701  DATE WRITTEN */
		/*   890531  Changed all specific intrinsics to generic.  (WRB) */
		/*   890531  REVISION DATE from Version 3.2 */
		/*   891214  Prologue converted to Version 4.0 format.  (BAB) */
		/*   900315  CALLs to XERROR changed to CALLs to XERMSG.  (THJ) */
		/*   920618  Removed space from variable names.  (RWC, WRB) */
		/*   121094  f2c cleanup, made thread-safe, converted erfc -> erfcx (SGJ) */
		/* END PROLOGUE  DERFCX */
		/* FIRST EXECUTABLE STATEMENT  DERFCX */
		
#if 0
		/* Compute limits on the fly, which is somewhat wasteful since
		 they will be the same for every call, but the alternative (as
		 in the original Fortran code) of using static variables is not
		 thread-safe.  (We could use thread-local static variables, but
		 that is nonstandard.)  In practice, all machines nowadays use
		 IEEE-754 arithmetic, so we should be able to just hard-code these. */
		double xsml;
		int nterf;
		double sqeps;
		int nterc2;
		int nterfc;
		{
			double d1mach3 = DBL_EPSILON/FLT_RADIX;
			double eta = d1mach3 * .1;
			nterf = initds(erfcs, 21, eta);
			nterfc = initds(erfccs, 59, eta);
			nterc2 = initds(erc2cs, 49, eta);
			
			xsml = -sqrt(-log(sqrtpi * DBL_EPSILON/FLT_RADIX));
			sqeps = sqrt(DBL_EPSILON/FLT_RADIX * 2.);
		}
#else /* hard-code IEEE values, so as to be re-entrant */
		const int nterf = 12;
		const int nterfc = 25;
		const int nterc2 = 24;
		const double xsml = -6.01368735691775061802858;
		const double sqeps = 1.490116119384765625e-08;
#endif
		
		if (x > xsml) {
			goto L20;
		}
		
		/* ERFC(X) = 1.0 - ERF(X)  FOR  X .LT. XSML */
		
		return 2. 
		/* SGJ: */ * exp(x*x);
		
	L20:
		/* SGJ: removed x > xmax test, in which case derfc returned 0 (at L40).
		 This was to account for underflow, which no longer occurs for erfcx.
		 {
		 double txmax;
		 txmax = sqrt(-log(sqrtpi * DBL_MIN));
		 xmax = txmax - log(txmax) * .5 / txmax - .01;
		 if (x > xmax) {
		 goto L40;
		 }
		 } 
		 */
		y = fabs(x);
		if (y > 1.) {
			goto L30;
		}
		
		/* ERFC(X) = 1.0 - ERF(X)  FOR ABS(X) .LE. 1.0 */
		
		if (y < sqeps) {
			return 1. - x * 2. / sqrtpi; /* note: exp(x*x) = 1 to mach. prec. */
		}
		else {
			d__1 = x * 2. * x - 1.;
			return (1. - x * (dcsevl_(d__1, erfcs, nterf) + 1.))
			/* SGJ: */ * exp(x*x); 
			/* TODO: should be more efficient to use a Chebyshev approximation
			 for erfcx directly, to avoid computation of exp(x*x) */
		}
		
		/* ERFC(X) = 1.0 - ERF(X)  FOR  1.0 .LT. ABS(X) .LE. XMAX */
		
	L30:
		y *= y;
		{
			double ret;
			if (y <= 4.) {
				d__1 = (8. / y - 5.) / 3.;
				ret = (dcsevl_(d__1, erc2cs, nterc2) + .5) / fabs(x);
				/* SGJ: removed exp(-y) factor from erfc */
			}
			else {
				d__1 = 8. / y - 1.;
				ret = (dcsevl_(d__1, erfccs, nterfc) + .5) / fabs(x);
				/* SGJ: removed exp(-y) factor from erfc */
			}
			if (x < 0.) {
				ret = 2. * exp(y) - ret; /* SGJ: added exp(y) factor */
			}
			return ret;
		}
	} /* derfcx */
	
	/////////////////////////////////////////////////////////////////////////
	
	// Compile with -DFADDEEVA_W_TEST to compile a little test program
#ifdef FADDEEVA_W_TEST
	
#include <cstdio>
	
#define NTST 22
	
	int main(void) {
		double complex z[NTST] = {
			double complex(624.2,-0.26123),
			double complex(-0.4,3.),
			double complex(0.6,2.),
			double complex(-1.,1.),
			double complex(-1.,-9.),
			double complex(-1.,9.),
			double complex(-0.0000000234545,1.1234),
			double complex(-3.,5.1),
			double complex(-53,30.1),
			double complex(0.0,0.12345),
			double complex(11,1),
			double complex(-22,-2),
			double complex(9,-28),
			double complex(21,-33),
			double complex(1e5,1e5),
			double complex(1e14,1e14),
			double complex(-3001,-1000),
			double complex(1e160,-1e159),
			double complex(-6.01,0.01),
			double complex(-0.7,-0.7),
			double complex(2.611780000000000e+01, 4.540909610972489e+03),
			double complex(0.8e7,0.3e7)
		};
		double complex w[NTST] = { // w(z), computed with WolframAlpha
			double complex(-3.78270245518980507452677445620103199303131110e-7,
							0.000903861276433172057331093754199933411710053155),
			double complex(0.1764906227004816847297495349730234591778719532788,
							-0.02146550539468457616788719893991501311573031095617),
			double complex(0.2410250715772692146133539023007113781272362309451,
							0.06087579663428089745895459735240964093522265589350),
			double complex(0.30474420525691259245713884106959496013413834051768,
							-0.20821893820283162728743734725471561394145872072738),
			double complex(7.317131068972378096865595229600561710140617977e34,
							8.321873499714402777186848353320412813066170427e34),
			double complex(0.0615698507236323685519612934241429530190806818395,
							-0.00676005783716575013073036218018565206070072304635),
			double complex(0.3960793007699874918961319170187598400134746631,
							-5.593152259116644920546186222529802777409274656e-9),
			double complex(0.08217199226739447943295069917990417630675021771804,
							-0.04701291087643609891018366143118110965272615832184),
			double complex(0.00457246000350281640952328010227885008541748668738,
							-0.00804900791411691821818731763401840373998654987934),
			double complex(0.8746342859608052666092782112565360755791467973338452,
							0.),
			double complex(0.00468190164965444174367477874864366058339647648741,
							0.0510735563901306197993676329845149741675029197050),
			double complex(-0.0023193175200187620902125853834909543869428763219,
							-0.025460054739731556004902057663500272721780776336),
			double complex(9.11463368405637174660562096516414499772662584e304,
							3.97101807145263333769664875189354358563218932e305),
			double complex(-4.4927207857715598976165541011143706155432296e281,
							-2.8019591213423077494444700357168707775769028e281),
			double complex(2.820947917809305132678577516325951485807107151e-6,
							2.820947917668257736791638444590253942253354058e-6),
			double complex(2.82094791773878143474039725787438662716372268e-15,
							2.82094791773878143474039725773333923127678361e-15),
			double complex(-0.0000563851289696244350147899376081488003110150498,
							-0.000169211755126812174631861529808288295454992688),
			double complex(-5.586035480670854326218608431294778077663867e-162,
							5.586035480670854326218608431294778077663867e-161),
			double complex(0.00016318325137140451888255634399123461580248456,
							-0.095232456573009287370728788146686162555021209999),
			double complex(0.69504753678406939989115375989939096800793577783885,
							-1.8916411171103639136680830887017670616339912024317),
			double complex(0.0001242418269653279656612334210746733213167234822,
							7.145975826320186888508563111992099992116786763e-7),
			double complex(2.318587329648353318615800865959225429377529825e-8,
							6.182899545728857485721417893323317843200933380e-8)
			
		};
		for (int i = 0; i < NTST; ++i) {
			double complex fw = Faddeeva_w(z[i],0.);
			double err = abs((fw - w[i]) / w[i]);
			printf("w(%g%+gi) = %g%+gi (vs. %g%+gi), rel. err. = %0.2g\n",
				   real(z[i]),imag(z[i]), real(fw),imag(fw), real(w[i]),imag(w[i]),
				   err);
			if (err > 1e-13) {
				printf("FAILURE -- relative error %g too large!\n", err);
				return 1;
			}
		}
		printf("SUCCESS\n");
		return 0;
	}
	
#endif