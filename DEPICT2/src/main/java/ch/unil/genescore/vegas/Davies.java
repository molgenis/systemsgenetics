/*******************************************************************************
 * Copyright (c) 2015 David Lamparter, Daniel Marbach
 *
 * Permission is hereby granted, free of charge, to any person obtaining a copy
 * of this software and associated documentation files (the "Software"), to deal
 * in the Software without restriction, including without limitation the rights
 * to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
 * copies of the Software, and to permit persons to whom the Software is
 * furnished to do so, subject to the following conditions:
 *
 * The above copyright notice and this permission notice shall be included in
 * all copies or substantial portions of the Software.
 *
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
 * IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
 * FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
 * AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
 * LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
 * OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN
 * THE SOFTWARE.
 *******************************************************************************/
package ch.unil.genescore.vegas;


/**
 * Refactored from R package CompQuadForm by Pierre Lafaye de Micheaux and Pierre Duchesne. See bottom of file for original code.
 * 
 * Distribution function (survival function in fact) of quadratic forms in normal variables using Davies algorithm.
 *
 */
public class Davies implements WeightedChisquareAlgorithm {

	// ARGUMENTS
	private final double log28 = .0866;  /*  log(2.0) / 8.0  */

	private double sigsq, lmax, lmin, mean, c;
	private double intl, ersm;
	private int count, r, lim;  
	private boolean ndtsrt, fail;
	private int n[], th[]; 
	private double lb[], nc[];
	double sigma_ = 0;

	// //	void  qfc(double* lb1, double* nc1, int* n1, int *r1, double *sigma, double *c1, int *lim1, double *acc, double* trace, int* ifault, double *res)

	//double[] lb1;
	//double[] nc1;
	//int[] n1;
	//int[] r1;
	//double[] sigma;
	//double[] c1;
	//int[] lim1;
	//double[] acc;
	double[] trace;
	
	/** Value point at which the survival function is to be evaluated */
	private double q_ = -1;
	/** Distinct non-zero characteristic roots of A\Sigma */
	private double[] lambda_ = null;
	/** Respective orders of multiplicity of the lambdas (degrees of the chi2 variables) */
	private int[] h_ = null;
	/** Non-centrality parameters of the chi2 variables */
	private double[] delta_ = null;
	/** Error bound. R documentation: error bound. Suitable values for acc range from 0.001 to 0.00005 */
	private double acc_ = -1;
	/**
	 * Maximum number of integration terms. Realistic values for lim range from 1000 if the procedure is 
	 * to be called repeatedly up to 50 000 if it is to be called only occasionally 
	 */
	private int lim_ = -1;
	
	/** The result */
	private double qfval = -1;
	/** Error code */
	private int ifault = 0;
	
	public class NonCovergeneceException extends RuntimeException {
		private static final long serialVersionUID = 7845175270361304271L;
	};
	

	// ============================================================================
	// PUBLIC METHODS
	
//	davies <- function (q, lambda, h = rep(1, length(lambda)), delta = rep(0, 
//		    length(lambda)), sigma = 0, lim = 10000, acc = 1e-04) 
//		{
//		    r <- length(lambda)
//		    if (length(h) != r) 
//		        stop("lambda and h should have the same length!")
//		    if (length(delta) != r) 
//		        stop("lambda and delta should have the same length!")
//		    out <- .C("qfc", lambdas = as.double(lambda), noncentral = as.double(delta), 
//		        df = as.integer(h), r = as.integer(r), sigma = as.double(sigma), 
//		        q = as.double(q), lim = as.integer(lim), acc = as.double(acc), 
//		        trace = as.double(rep(0, 7)), ifault = as.integer(0), 
//		        res = as.double(0), PACKAGE = "CompQuadForm")
//		    out$res <- 1 - out$res
//		    return(list(trace = out$trace, ifault = out$ifault, Qq = out$res))
//		}
//		<environment: namespace:CompQuadForm>

	/** Constructor */
	public Davies(double[] lambda) {

		lambda_ = lambda;
		
		h_ = new int[lambda.length];
		for (int i=0; i<lambda.length; i++)
			h_[i] = 1;
		
		delta_ = new double[lambda.length];
		for (int i=0; i<lambda.length; i++)
			delta_[i] = 0;
		
		acc_ = 1e-20;
		lim_ = 50000;
	}

		
	// ----------------------------------------------------------------------------

	/** Compute P(Q > q) */
	public double probQsupx(double q) {

		// initialize
		q_ = q;
		ifault = 0;
		qfval = -41;
		
		trace = new double[7];
		for (int i=0; i<7; i++)
			trace[i] = 0.0;
		
		// compute
		qfc();
		qfval = 1 - qfval;
		return qfval;
	}

	public int getIfault() {
		return ifault;
	}
	
	// ============================================================================
	// PRIVATE FUNCTIONS


	/* to avoid underflows  */
	private double exp1(double x) { 
		return x < -50.0 ? 0.0 : Math.exp(x); 
	}

	/*  count number of calls to errbd, truncation, cfe */
	private void counter() {
	      
		count = count + 1;
		if ( count > lim ) {
			ifault = 4;
			throw new NonCovergeneceException();
		}
	}

	private double square(double x)  { return x*x; }

	private double cube(double x)  { return x*x*x; }

	/* if (first) log(1 + x) ; else  log(1 + x) - x */
	private double log1(double x, boolean first) {

		if (Math.abs(x) > 0.1)
		{
			return (first ? Math.log(1.0 + x) : (Math.log(1.0 + x) - x));
		}
		else
		{
			double s, s1, term, y, k;
			y = x / (2.0 + x);  term = 2.0 * cube(y);  k = 3.0;
			s = (first ? 2.0 : - x) * y;
			y = square(y);
			for (s1 = s + term / k; s1 != s; s1 = s + term / k)
			{ k = k + 2.0; term = term * y; s = s1; }
			return s;
		}
	}

	/* find order of absolute values of lb */
	private void order() {
	      
		int j, k; double lj;
	      //extern double *lb; extern int *th; extern int r; extern BOOL ndtsrt;
	      for ( j=0; j<r; j++ )
	      {
	         lj = Math.abs(lb[j]);
	         boolean gotol1 = false;
	         for (k = j-1; k>=0; k--)
	         {
	            if ( lj > Math.abs(lb[th[k]]) )  th[k + 1] = th[k];
	            //else goto l1;
	            else {
	            	gotol1 = true;
	            	break;
	            }
	         }
	         if (!gotol1)
	        	 k = -1;
	      //l1 :
	         th[k + 1] = j;
	      }
	      ndtsrt = false;
	   }


	   private double errbd(double u, double[] cx)
	   /*  find bound on tail probability using mgf, cutoff
	      point returned to *cx */
	   {
	      double sum1, lj, ncj, x, y, xconst; int j, nj;
	      //extern double sigsq,*lb,*nc; extern int *n; extern int r;
	      counter();
	      xconst = u * sigsq;  sum1 = u * xconst;  u = 2.0 * u;
	      for (j=r-1; j>=0; j--)
	      {
	         nj = n[j]; lj = lb[j]; ncj = nc[j];
	         x = u * lj; y = 1.0 - x;
	         xconst = xconst + lj * (ncj / y + nj) / y;
	         sum1 = sum1 + ncj * square(x / y)
	            + nj * (square(x) / y + log1(-x, false ));
	      }
	      cx[0] = xconst; return exp1(-0.5 * sum1);
	   }

	   private double  ctff(double accx, double[] upn)
	   /*  find ctff so that p(qf > ctff) < accx  if (upn > 0,
	       p(qf < ctff) < accx otherwise */
	   {
	      double u1, u2, u, rb, c1;// c2;
	      double[] c2 = new double[1];
	      double[] xconst = new double[1];
	      //extern double lmin,lmax,mean;
	      u2 = upn[0];   u1 = 0.0;  c1 = mean;
	      rb = 2.0 * ((u2 > 0.0) ? lmax : lmin);
	      for (u = u2 / (1.0 + u2 * rb); errbd(u, c2) > accx; 
	         u = u2 / (1.0 + u2 * rb))
	      {
	         u1 = u2;  c1 = c2[0];  u2 = 2.0 * u2;
	      }
	      for (u = (c1 - mean) / (c2[0] - mean); u < 0.9;
	         u = (c1 - mean) / (c2[0] - mean))
	      {
	         u = (u1 + u2) / 2.0;
	         if (errbd(u / (1.0 + u * rb), xconst) > accx)
	            {  u1 = u; c1 = xconst[0];  }
	         else
	            {  u2 = u;  c2[0] = xconst[0]; }
	      }
	      upn[0] = u2; return c2[0];
	   }

	   private double truncation(double u, double tausq)
	   /* bound integration error due to truncation at u */
	   {
	      double sum1, sum2, prod1, prod2, prod3, lj, ncj,
	         x, y, err1, err2;
	      int j, nj, s;
	      //extern double sigsq,*lb,*nc; extern int *n; extern int r;

	      counter();
	      sum1  = 0.0; prod2 = 0.0;  prod3 = 0.0;  s = 0;
	      sum2 = (sigsq + tausq) * square(u); prod1 = 2.0 * sum2;
	      u = 2.0 * u;
	      for (j=0; j<r; j++ )
	      {
	         lj = lb[j];  ncj = nc[j]; nj = n[j];
	         x = square(u * lj);
	         sum1 = sum1 + ncj * x / (1.0 + x);
	         if (x > 1.0)
	         {
	            prod2 = prod2 + nj * Math.log(x);
	            prod3 = prod3 + nj * log1(x, true );
	            s = s + nj;
	         }
	         else  prod1 = prod1 + nj * log1(x, true );
	      }
	      sum1 = 0.5 * sum1;
	      prod2 = prod1 + prod2;  prod3 = prod1 + prod3;
	      x = exp1(-sum1 - 0.25 * prod2) / Math.PI;
	      y = exp1(-sum1 - 0.25 * prod3) / Math.PI;
	      err1 =  ( s  ==  0 )  ? 1.0 : x * 2.0 / s;
	      err2 =  ( prod3 > 1.0 )  ? 2.5 * y : 1.0;
	      if (err2 < err1) err1 = err2;
	      x = 0.5 * sum2;
	      err2 =  ( x  <=  y )  ? 1.0  : y / x;
	      return  ( err1 < err2 )  ? err1  :  err2;
	   }

	   private void findu(double[] utx, double accx)
	   /*  find u such that truncation(u) < accx and truncation(u / 1.2) > accx */
	   {
	      double u, ut; int i;
	      final double divis[]= { 2.0,1.4,1.2,1.1 };
	      ut = utx[0]; u = ut / 4.0;
	      if ( truncation(u, 0.0) > accx )
	      {
	         for ( u = ut; truncation(u, 0.0) > accx; u = ut) ut = ut * 4.0;
	      }
	      else
	      {
	         ut = u;
	         for ( u = u / 4.0; truncation(u, 0.0) <=  accx; u = u / 4.0 )
	         ut = u;
	      }
	      for ( i=0;i<4;i++)
	         { u = ut/divis[i]; if ( truncation(u, 0.0)  <=  accx )  ut = u; }
	      utx[0] = ut;
	   }


	   private void integrate(int nterm, double interv, double tausq, boolean mainx)
	   /*  carry out integration with nterm terms, at stepsize
	      interv.  if (! mainx) multiply integrand by
	         1.0-exp(-0.5*tausq*u^2) */
	   {
	      double inpi, u, sum1, sum2, sum3, x, y, z;
	      int k, j, nj;
	      //extern double intl,ersm; extern double sigsq,c;
	      //extern int *n; extern double *lb,*nc; extern int r;
	      inpi = interv / Math.PI;
	      for ( k = nterm; k>=0; k--)
	      {
	         u = (k + 0.5) * interv;
	         sum1 = - 2.0 * u * c;  sum2 = Math.abs(sum1);
	         sum3 = - 0.5 * sigsq * square(u);
	         for ( j = r-1; j>=0; j--)
	         {
	            nj = n[j];  x = 2.0 * lb[j] * u;  y = square(x);
	            sum3 = sum3 - 0.25 * nj * log1(y, true );
	            y = nc[j] * x / (1.0 + y);
	            z = nj * Math.atan(x) + y;
	            sum1 = sum1 + z;   sum2 = sum2 + Math.abs(z);
	            sum3 = sum3 - 0.5 * x * y;
	         }
	         x = inpi * exp1(sum3) / u;
		 if ( !  mainx )
	         x = x * (1.0 - exp1(-0.5 * tausq * square(u)));
	         sum1 = Math.sin(0.5 * sum1) * x;  sum2 = 0.5 * sum2 * x;
	         intl = intl + sum1; ersm = ersm + sum2;
	      }
	   }

	   private double cfe(double x)
	   /*  coef of tausq in error when convergence factor of
	      exp1(-0.5*tausq*u^2) is used when df is evaluated at x */
	   {
	      double axl, axl1, axl2, sxl, sum1, lj; int j, k, t;
	      //extern BOOL ndtsrt,fail; extern int *th,*n; extern double *lb,*nc;
	      //extern int r;
	      counter();
	      if (ndtsrt) order();
	      axl = Math.abs(x);  sxl = (x>0.0) ? 1.0 : -1.0;  sum1 = 0.0;
	      for ( j = r-1; j>=0; j-- ) { 
	    	  t = th[j];
	    	  if ( lb[t] * sxl > 0.0 )
	    	  {
	    		  lj = Math.abs(lb[t]);
	    		  axl1 = axl - lj * (n[t] + nc[t]);  axl2 = lj / log28;
	    		  if ( axl1 > axl2 )  
	    			  axl = axl1; 
	    		  else
	    		  {
	    			  if ( axl > axl2 )
	    				  axl = axl2;
	    			  sum1 = (axl - axl1) / lj;
	    			  for ( k = j-1; k>=0; k--)
	    				  sum1 = sum1 + (n[th[k]] + nc[th[k]]);
	    			  //goto  l;
	    			  break;
	    		  }
	    	  }
	      }
	      //l:
	      if (sum1 > 100.0) { 
	    	  fail = true; 
	    	  return 1.0; 
	      } else
	    	  return Math.pow(2.0,(sum1 / 4.0)) / (Math.PI * square(axl));
	   }


	private void  qfc()

	/*  distribution function of a linear combination of non-central
	   chi-squared random variables :

	input:
	   lb[j]            coefficient of j-th chi-squared variable
	   nc[j]            non-centrality parameter
	   n[j]             degrees of freedom
	   j = 0, 2 ... r-1
	   sigma            coefficient of standard normal variable
	   c1                point at which df is to be evaluated
	   lim              maximum number of terms in integration
	   acc              maximum error

	output:
	   ifault = 1       required accuracy NOT achieved
	            2       round-off error possibly significant
	            3       invalid parameters
	            4       unable to locate integration parameters
	            5       out of memory

	   trace[0]         absolute sum
	   trace[1]         total number of integration terms
	   trace[2]         number of integrations
	   trace[3]         integration interval in final integration
	   trace[4]         truncation point in initial integration
	   trace[5]         s.d. of initial convergence factor
	   trace[6]         cycles to locate integration parameters     */

	{
	      int j, nj, nt, ntm;  double acc1, almx, xlim, xnt, xntm;
	      double tausq, sd, intv, intv1, x, d1, d2, lj, ncj;
	      double[] utx = new double[1];
	      double[] up = new double[1];
	      double[] un = new double[1];
	      //extern double sigsq, lmax, lmin, mean;
	      //extern double intl,ersm;
	      //extern int r,lim; extern double c;
	      //extern int *n,*th; extern double *lb,*nc;
	      qfval = -1.0;
	      final int rats[]={1,2,4,8};

	      //if (setjmp(env) != 0) { *ifault=4; goto endofproc; }
	      r=lambda_.length; lim=lim_; c=q_;
	      n=h_; lb=lambda_; nc=delta_;
	      for ( j = 0; j<7; j++ )  trace[j] = 0.0;
	      ifault = 0; count = 0;
	      intl = 0.0; ersm = 0.0;
	      qfval = -1.0; acc1 = acc_; ndtsrt = true;  fail = false;
	      xlim = (double)lim;
	      th= new int[r]; //(int*)malloc(r*(sizeof(int)));

	      /* find mean, sd, max and min of lb,
	         check that parameter values are valid */
	      sigsq = square(sigma_); sd = sigsq;
	      lmax = 0.0; lmin = 0.0; mean = 0.0;
	      for (j=0; j<r; j++ )
	      {
	         nj = n[j];  lj = lb[j];  ncj = nc[j];
	         if ( nj < 0  ||  ncj < 0.0 ) { ifault = 3;  return; } //goto  endofproc;  }
	         sd  = sd  + square(lj) * (2 * nj + 4.0 * ncj);
	         mean = mean + lj * (nj + ncj);
	         if (lmax < lj) lmax = lj ; else if (lmin > lj) lmin = lj;
	      }
	      if ( sd == 0.0  )
	      {  qfval = (c > 0.0) ? 1.0 : 0.0; return; } //goto  endofproc;  }
	      if ( lmin == 0.0 && lmax == 0.0 && sigma_ == 0.0 )
	         { ifault = 3;  return; } //goto  endofproc;  }
	      sd = Math.sqrt(sd);
	      almx = (lmax < - lmin) ? - lmin : lmax;

	      /* starting values for findu, ctff */
	      utx[0] = 16.0 / sd;  up[0] = 4.5 / sd;  un[0] = - up[0];
	      /* truncation point with no convergence factor */
	      findu(utx, .5 * acc1);
	      /* does convergence factor help */
	      if (c != 0.0  && (almx > 0.07 * sd))
	      {
	         tausq = .25 * acc1 / cfe(c);
	         if (fail) fail = false ;
	         else if (truncation(utx[0], tausq) < .2 * acc1)
	         {
	            sigsq = sigsq + tausq;
	            findu(utx, .25 * acc1);
	            trace[5] = Math.sqrt(tausq);
	         }
	      }
	      trace[4] = utx[0];  acc1 = 0.5 * acc1;
	      
	      /* find RANGE of distribution, quit if outside this */
	   //l1:
	      while (true) {
	    	  d1 = ctff(acc1, up) - c;
	    	  if (d1 < 0.0) { qfval = 1.0; return; }
	    	  d2 = c - ctff(acc1, un);
	    	  if (d2 < 0.0) { qfval = 0.0; return; }
	    	  /* find integration interval */
	    	  intv = 2.0 * Math.PI / ((d1 > d2) ? d1 : d2);
	    	  /* calculate number of terms required for main and
	         auxillary integrations */
	    	  xnt = utx[0] / intv;  xntm = 3.0 / Math.sqrt(acc1);
	    	  if (!(xnt > xntm * 1.5))
	    		  break;

	    	  /* parameters for auxillary integration */
	    	  if (xntm > xlim) { ifault = 1; return; }
	    	  ntm = (int)Math.floor(xntm+0.5);
	    	  intv1 = utx[0] / ntm;  x = 2.0 * Math.PI / intv1;
	    	  if (x <= Math.abs(c)) break;
	    	  /* calculate convergence factor */
	    	  tausq = .33 * acc1 / (1.1 * (cfe(c - x) + cfe(c + x)));
	    	  if (fail) break;
	    	  acc1 = .67 * acc1;
	    	  /* auxillary integration */
	    	  integrate(ntm, intv1, tausq, false );
	    	  xlim = xlim - xntm;  sigsq = sigsq + tausq;
	    	  trace[2] = trace[2] + 1; trace[1] = trace[1] + ntm + 1;
	    	  /* find truncation point with new convergence factor */
	    	  findu(utx, .25 * acc1);  acc1 = 0.75 * acc1;
	    	  //goto l1;
	      }

	      /* main integration */
//	   l2:
	      trace[3] = intv;
	      if (xnt > xlim) { ifault = 1; return; }
	      nt = (int)Math.floor(xnt+0.5);
	      integrate(nt, intv, 0.0, true );
	      trace[2] = trace[2] + 1; trace[1] = trace[1] + nt + 1;
	      qfval = 0.5 - intl;
	      trace[0] = ersm;

	      /* test whether round-off error could be significant
	         allow for radix 8 or 16 machines */
	      up[0]=ersm; x = up[0] + acc_ / 10.0;
	      for (j=0;j<4;j++) { if (rats[j] * x == rats[j] * up[0]) ifault = 2; }

	   //endofproc :
	      //trace[6] = (double)count;
	      //res[0] = qfval;
	      //return;
	}


	}




	
	
	
	
	
//	#define UseDouble 0             /* all floating point double */
//
//	#include <stdio.h>
//	#include <stdlib.h>
//	#include <math.h>
//	#include <setjmp.h>
//
//	#define TRUE  1
//	#define FALSE 0
//	typedef int BOOL;
//
//
//	#define pi 3.14159265358979
//	#define log28 .0866  /*  log(2.0) / 8.0  */
//
//	extern "C" {
//
//	static double sigsq, lmax, lmin, mean, c;
//	static double intl, ersm;
//	static int count, r, lim;  static BOOL ndtsrt, fail;
//	static int *n,*th; static double *lb,*nc;
//	static jmp_buf env;
//
//
//
//	static double exp1(double x)               /* to avoid underflows  */
//	{ return x < -50.0 ? 0.0 : exp(x); }
//
//	   static void counter(void)
//	   /*  count number of calls to errbd, truncation, cfe */
//	   {
//	      extern int count,lim;
//	      count = count + 1;
//	      if ( count > lim ) longjmp(env,1);
//	   }
//
//	   static double square(double x)  { return x*x; }
//
//	   static double cube(double x)  { return x*x*x; }
//
//	   static double  log1(double x, BOOL first)
//	   /* if (first) log(1 + x) ; else  log(1 + x) - x */
//	   {
//	      if (fabs(x) > 0.1)
//	      {
//	         return (first ? log(1.0 + x) : (log(1.0 + x) - x));
//	      }
//	      else
//	      {
//	         double s, s1, term, y, k;
//	         y = x / (2.0 + x);  term = 2.0 * cube(y);  k = 3.0;
//	         s = (first ? 2.0 : - x) * y;
//	         y = square(y);
//	         for (s1 = s + term / k; s1 != s; s1 = s + term / k)
//	         { k = k + 2.0; term = term * y; s = s1; }
//	         return s;
//	      }
//	  }
//
//	  static void order(void)
//	  /* find order of absolute values of lb */
//	  {
//	      int j, k; double lj;
//	      extern double *lb; extern int *th; extern int r; extern BOOL ndtsrt;
//	      for ( j=0; j<r; j++ )
//	      {
//	         lj = fabs(lb[j]);
//	         for (k = j-1; k>=0; k--)
//	         {
//	            if ( lj > fabs(lb[th[k]]) )  th[k + 1] = th[k];
//	            else goto l1;
//	         }
//	         k = -1;
//	      l1 :
//	         th[k + 1] = j;
//	      }
//	      ndtsrt = FALSE;
//	   }
//
//
//	   static double   errbd(double u, double* cx)
//	   /*  find bound on tail probability using mgf, cutoff
//	      point returned to *cx */
//	   {
//	      double sum1, lj, ncj, x, y, xconst; int j, nj;
//	      extern double sigsq,*lb,*nc; extern int *n; extern int r;
//	      counter();
//	      xconst = u * sigsq;  sum1 = u * xconst;  u = 2.0 * u;
//	      for (j=r-1; j>=0; j--)
//	      {
//	         nj = n[j]; lj = lb[j]; ncj = nc[j];
//	         x = u * lj; y = 1.0 - x;
//	         xconst = xconst + lj * (ncj / y + nj) / y;
//	         sum1 = sum1 + ncj * square(x / y)
//	            + nj * (square(x) / y + log1(-x, FALSE ));
//	      }
//	      *cx = xconst; return exp1(-0.5 * sum1);
//	   }
//
//	   static double  ctff(double accx, double* upn)
//	   /*  find ctff so that p(qf > ctff) < accx  if (upn > 0,
//	       p(qf < ctff) < accx otherwise */
//	   {
//	      double u1, u2, u, rb, xconst, c1, c2;
//	      extern double lmin,lmax,mean;
//	      u2 = *upn;   u1 = 0.0;  c1 = mean;
//	      rb = 2.0 * ((u2 > 0.0) ? lmax : lmin);
//	      for (u = u2 / (1.0 + u2 * rb); errbd(u, &c2) > accx; 
//	         u = u2 / (1.0 + u2 * rb))
//	      {
//	         u1 = u2;  c1 = c2;  u2 = 2.0 * u2;
//	      }
//	      for (u = (c1 - mean) / (c2 - mean); u < 0.9;
//	         u = (c1 - mean) / (c2 - mean))
//	      {
//	         u = (u1 + u2) / 2.0;
//	         if (errbd(u / (1.0 + u * rb), &xconst) > accx)
//	            {  u1 = u; c1 = xconst;  }
//	         else
//	            {  u2 = u;  c2 = xconst; }
//	      }
//	      *upn = u2; return c2;
//	   }
//
//	   static double truncation(double u, double tausq)
//	   /* bound integration error due to truncation at u */
//	   {
//	      double sum1, sum2, prod1, prod2, prod3, lj, ncj,
//	         x, y, err1, err2;
//	      int j, nj, s;
//	      extern double sigsq,*lb,*nc; extern int *n; extern int r;
//
//	      counter();
//	      sum1  = 0.0; prod2 = 0.0;  prod3 = 0.0;  s = 0;
//	      sum2 = (sigsq + tausq) * square(u); prod1 = 2.0 * sum2;
//	      u = 2.0 * u;
//	      for (j=0; j<r; j++ )
//	      {
//	         lj = lb[j];  ncj = nc[j]; nj = n[j];
//	         x = square(u * lj);
//	         sum1 = sum1 + ncj * x / (1.0 + x);
//	         if (x > 1.0)
//	         {
//	            prod2 = prod2 + nj * log(x);
//	            prod3 = prod3 + nj * log1(x, TRUE );
//	            s = s + nj;
//	         }
//	         else  prod1 = prod1 + nj * log1(x, TRUE );
//	      }
//	      sum1 = 0.5 * sum1;
//	      prod2 = prod1 + prod2;  prod3 = prod1 + prod3;
//	      x = exp1(-sum1 - 0.25 * prod2) / pi;
//	      y = exp1(-sum1 - 0.25 * prod3) / pi;
//	      err1 =  ( s  ==  0 )  ? 1.0 : x * 2.0 / s;
//	      err2 =  ( prod3 > 1.0 )  ? 2.5 * y : 1.0;
//	      if (err2 < err1) err1 = err2;
//	      x = 0.5 * sum2;
//	      err2 =  ( x  <=  y )  ? 1.0  : y / x;
//	      return  ( err1 < err2 )  ? err1  :  err2;
//	   }
//
//	   static void findu(double* utx, double accx)
//	   /*  find u such that truncation(u) < accx and truncation(u / 1.2) > accx */
//	   {
//	      double u, ut; int i;
//	      static double divis[]={2.0,1.4,1.2,1.1};
//	      ut = *utx; u = ut / 4.0;
//	      if ( truncation(u, 0.0) > accx )
//	      {
//	         for ( u = ut; truncation(u, 0.0) > accx; u = ut) ut = ut * 4.0;
//	      }
//	      else
//	      {
//	         ut = u;
//	         for ( u = u / 4.0; truncation(u, 0.0) <=  accx; u = u / 4.0 )
//	         ut = u;
//	      }
//	      for ( i=0;i<4;i++)
//	         { u = ut/divis[i]; if ( truncation(u, 0.0)  <=  accx )  ut = u; }
//	      *utx = ut;
//	   }
//
//
//	   static void integrate(int nterm, double interv, double tausq, BOOL mainx)
//	   /*  carry out integration with nterm terms, at stepsize
//	      interv.  if (! mainx) multiply integrand by
//	         1.0-exp(-0.5*tausq*u^2) */
//	   {
//	      double inpi, u, sum1, sum2, sum3, x, y, z;
//	      int k, j, nj;
//	      extern double intl,ersm; extern double sigsq,c;
//	      extern int *n; extern double *lb,*nc; extern int r;
//	      inpi = interv / pi;
//	      for ( k = nterm; k>=0; k--)
//	      {
//	         u = (k + 0.5) * interv;
//	         sum1 = - 2.0 * u * c;  sum2 = fabs(sum1);
//	         sum3 = - 0.5 * sigsq * square(u);
//	         for ( j = r-1; j>=0; j--)
//	         {
//	            nj = n[j];  x = 2.0 * lb[j] * u;  y = square(x);
//	            sum3 = sum3 - 0.25 * nj * log1(y, TRUE );
//	            y = nc[j] * x / (1.0 + y);
//	            z = nj * atan(x) + y;
//	            sum1 = sum1 + z;   sum2 = sum2 + fabs(z);
//	            sum3 = sum3 - 0.5 * x * y;
//	         }
//	         x = inpi * exp1(sum3) / u;
//		 if ( !  mainx )
//	         x = x * (1.0 - exp1(-0.5 * tausq * square(u)));
//	         sum1 = sin(0.5 * sum1) * x;  sum2 = 0.5 * sum2 * x;
//	         intl = intl + sum1; ersm = ersm + sum2;
//	      }
//	   }
//
//	   static double cfe(double x)
//	   /*  coef of tausq in error when convergence factor of
//	      exp1(-0.5*tausq*u^2) is used when df is evaluated at x */
//	   {
//	      double axl, axl1, axl2, sxl, sum1, lj; int j, k, t;
//	      extern BOOL ndtsrt,fail; extern int *th,*n; extern double *lb,*nc;
//	      extern int r;
//	      counter();
//	      if (ndtsrt) order();
//	      axl = fabs(x);  sxl = (x>0.0) ? 1.0 : -1.0;  sum1 = 0.0;
//	      for ( j = r-1; j>=0; j-- )
//	      { t = th[j];
//	         if ( lb[t] * sxl > 0.0 )
//	         {
//	            lj = fabs(lb[t]);
//	            axl1 = axl - lj * (n[t] + nc[t]);  axl2 = lj / log28;
//	            if ( axl1 > axl2 )  axl = axl1  ; else
//	            {
//	               if ( axl > axl2 )  axl = axl2;
//	               sum1 = (axl - axl1) / lj;
//	               for ( k = j-1; k>=0; k--)
//	               sum1 = sum1 + (n[th[k]] + nc[th[k]]);
//	               goto  l;
//	            }
//	         }
//	      }
//	   l:
//	      if (sum1 > 100.0)
//	      { fail = TRUE; return 1.0; } else
//	      return pow(2.0,(sum1 / 4.0)) / (pi * square(axl));
//	   }
//
//
//	void  qfc(double* lb1, double* nc1, int* n1, int *r1, double *sigma, double *c1, int *lim1, double *acc, double* trace, int* ifault, double *res)
//
//	/*  distribution function of a linear combination of non-central
//	   chi-squared random variables :
//
//	input:
//	   lb[j]            coefficient of j-th chi-squared variable
//	   nc[j]            non-centrality parameter
//	   n[j]             degrees of freedom
//	   j = 0, 2 ... r-1
//	   sigma            coefficient of standard normal variable
//	   c1                point at which df is to be evaluated
//	   lim              maximum number of terms in integration
//	   acc              maximum error
//
//	output:
//	   ifault = 1       required accuracy NOT achieved
//	            2       round-off error possibly significant
//	            3       invalid parameters
//	            4       unable to locate integration parameters
//	            5       out of memory
//
//	   trace[0]         absolute sum
//	   trace[1]         total number of integration terms
//	   trace[2]         number of integrations
//	   trace[3]         integration interval in final integration
//	   trace[4]         truncation point in initial integration
//	   trace[5]         s.d. of initial convergence factor
//	   trace[6]         cycles to locate integration parameters     */
//
//	{
//	      int j, nj, nt, ntm;  double acc1, almx, xlim, xnt, xntm;
//	      double utx, tausq, sd, intv, intv1, x, up, un, d1, d2, lj, ncj;
//	      extern double sigsq, lmax, lmin, mean;
//	      extern double intl,ersm;
//	      extern int r,lim; extern double c;
//	      extern int *n,*th; extern double *lb,*nc;
//	      double qfval = -1.0;
//	      static int rats[]={1,2,4,8};
//
//	      if (setjmp(env) != 0) { *ifault=4; goto endofproc; }
//	      r=r1[0]; lim=lim1[0]; c=c1[0];
//	      n=n1; lb=lb1; nc=nc1;
//	      for ( j = 0; j<7; j++ )  trace[j] = 0.0;
//	      *ifault = 0; count = 0;
//	      intl = 0.0; ersm = 0.0;
//	      qfval = -1.0; acc1 = acc[0]; ndtsrt = TRUE;  fail = FALSE;
//	      xlim = (double)lim;
//	      th=(int*)malloc(r*(sizeof(int)));
//	      if (! th) { *ifault=5;  goto  endofproc; } 
//
//	      /* find mean, sd, max and min of lb,
//	         check that parameter values are valid */
//	      sigsq = square(sigma[0]); sd = sigsq;
//	      lmax = 0.0; lmin = 0.0; mean = 0.0;
//	      for (j=0; j<r; j++ )
//	      {
//	         nj = n[j];  lj = lb[j];  ncj = nc[j];
//	         if ( nj < 0  ||  ncj < 0.0 ) { *ifault = 3;  goto  endofproc;  }
//	         sd  = sd  + square(lj) * (2 * nj + 4.0 * ncj);
//	         mean = mean + lj * (nj + ncj);
//	         if (lmax < lj) lmax = lj ; else if (lmin > lj) lmin = lj;
//	      }
//	      if ( sd == 0.0  )
//	      {  qfval = (c > 0.0) ? 1.0 : 0.0; goto  endofproc;  }
//	      if ( lmin == 0.0 && lmax == 0.0 && sigma[0] == 0.0 )
//	         { *ifault = 3;  goto  endofproc;  }
//	      sd = sqrt(sd);
//	      almx = (lmax < - lmin) ? - lmin : lmax;
//
//	      /* starting values for findu, ctff */
//	      utx = 16.0 / sd;  up = 4.5 / sd;  un = - up;
//	      /* truncation point with no convergence factor */
//	      findu(&utx, .5 * acc1);
//	      /* does convergence factor help */
//	      if (c != 0.0  && (almx > 0.07 * sd))
//	      {
//	         tausq = .25 * acc1 / cfe(c);
//	         if (fail) fail = FALSE ;
//	         else if (truncation(utx, tausq) < .2 * acc1)
//	         {
//	            sigsq = sigsq + tausq;
//	            findu(&utx, .25 * acc1);
//	            trace[5] = sqrt(tausq);
//	         }
//	      }
//	      trace[4] = utx;  acc1 = 0.5 * acc1;
//
//	      /* find RANGE of distribution, quit if outside this */
//	   l1:
//	      d1 = ctff(acc1, &up) - c;
//	      if (d1 < 0.0) { qfval = 1.0; goto endofproc; }
//	      d2 = c - ctff(acc1, &un);
//	      if (d2 < 0.0) { qfval = 0.0; goto endofproc; }
//	      /* find integration interval */
//	      intv = 2.0 * pi / ((d1 > d2) ? d1 : d2);
//	      /* calculate number of terms required for main and
//	         auxillary integrations */
//	      xnt = utx / intv;  xntm = 3.0 / sqrt(acc1);
//	      if (xnt > xntm * 1.5)
//	      {
//	         /* parameters for auxillary integration */
//	         if (xntm > xlim) { *ifault = 1; goto endofproc; }
//	         ntm = (int)floor(xntm+0.5);
//	         intv1 = utx / ntm;  x = 2.0 * pi / intv1;
//	         if (x <= fabs(c)) goto l2;
//	         /* calculate convergence factor */
//	         tausq = .33 * acc1 / (1.1 * (cfe(c - x) + cfe(c + x)));
//	         if (fail) goto l2;
//	         acc1 = .67 * acc1;
//	         /* auxillary integration */
//	         integrate(ntm, intv1, tausq, FALSE );
//	         xlim = xlim - xntm;  sigsq = sigsq + tausq;
//	         trace[2] = trace[2] + 1; trace[1] = trace[1] + ntm + 1;
//	         /* find truncation point with new convergence factor */
//	         findu(&utx, .25 * acc1);  acc1 = 0.75 * acc1;
//	         goto l1;
//	      }
//
//	      /* main integration */
//	   l2:
//	      trace[3] = intv;
//	      if (xnt > xlim) { *ifault = 1; goto endofproc; }
//	      nt = (int)floor(xnt+0.5);
//	      integrate(nt, intv, 0.0, TRUE );
//	      trace[2] = trace[2] + 1; trace[1] = trace[1] + nt + 1;
//	      qfval = 0.5 - intl;
//	      trace[0] = ersm;
//
//	      /* test whether round-off error could be significant
//	         allow for radix 8 or 16 machines */
//	      up=ersm; x = up + acc[0] / 10.0;
//	      for (j=0;j<4;j++) { if (rats[j] * x == rats[j] * up) *ifault = 2; }
//
//	   endofproc :
//	      free((char*)th);
//	      trace[6] = (double)count;
//	      res[0] = qfval;
//	      return;
//	}
//
//
//	}



			


