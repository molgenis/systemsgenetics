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

import org.apache.commons.math3.distribution.NormalDistribution;


/**
 * Refactored from R package CompQuadForm by Pierre Lafaye de Micheaux and Pierre Duchesne. See bottom of file for original code.
 * 
 * Distribution function (survival function in fact) of quadratic forms in normal variables using Fare- brothers algorithm.
 *
 * Algorithm AS 204 Appl. Statist. (1984) Vol. 33, No.3
 * ruben evaluates the probability that a positive definite quadratic form in Normal variates is less than a given value
 */
public class Farebrother implements WeightedChisquareAlgorithm  {

	// ARGUMENTS
	/** Value point at which the survival function is to be evaluated */
	private double q_ = -1;
	/** Distinct non-zero characteristic roots of A\Sigma */
	private double[] lambda_ = null;
	/** Respective orders of multiplicity of the lambdas (degrees of the chi2 variables) */
	private int[] h_ = null;
	/** Non-centrality parameters of the chi2 variables */
	private double[] delta_ = null;
	/** Accuracy requested */
	private double eps_ = -1;
	/** Max number of iterations */
	private int maxit_ = -1;
	/** if mode>0 then beta = mode * lambda_min otherwise beta = beta*B = 2/(1/lambda_min + 1/lambda_max) */
	private double mode_ = 1;
	
	/** Normal distribution */
	private NormalDistribution normal_ = new NormalDistribution(0.0, 1.0);
	
	/** The result */
	private double res_ = -1;
	/** Error code */
	private int ifault_ = 0;
	
	
	// ============================================================================
	// PUBLIC METHODS
	
	/** Constructor */
	public Farebrother(double[] lambda) {

		lambda_ = lambda;
		
		h_ = new int[lambda.length];
		for (int i=0; i<lambda.length; i++)
			h_[i] = 1;
		
		delta_ = new double[lambda.length];
		for (int i=0; i<lambda.length; i++)
			delta_[i] = 0;
		
		eps_ = 1e-16;
		maxit_ = 1000000;
		mode_ = -1.0;
	}

		
	// ----------------------------------------------------------------------------

	/** Compute P(Q > q) */
	public double probQsupx(double q) {

		// initialize
		q_ = q;
		ifault_ = 0;
		res_ = -41;
		
		// compute
		ruben();
		res_ = 1 - res_;
		if(res_<1E-15){
			res_=1E-15;
		}
		return res_;
	}

	
	// ============================================================================
	// PRIVATE FUNCTIONS

	/**
	 * ruben.cpp:
	 * Algorithm AS 204 Appl. Statist. (1984) Vol. 33, No.3
	 * ruben evaluates the probability that a positive definite quadratic form in Normal variates is less than a given value
	 */
	private void ruben() {
		
		ifault_ = -41; // I added this initialization --daniel
		
		//void ruben(double *lambda, int *mult, double *delta, int *n, double *c, double *mode, int *maxit, double *eps, double *dnsty, int *ifault, double *res) {
		// Initialized at 0 in the R code
		double dnsty = 0.0;
			    
		int i,k,m,j;
		//double pnorm(double q, double mean, double sd, int lower_tail, int log_p);
		double ao, aoinv, z, bbeta, eps2, hold, hold2, sum, sum1, dans, lans, pans, prbty, tol;
		double[] gamma = new double[lambda_.length];
		double[] theta = new double[lambda_.length];
		double[] a = new double[maxit_];
		double[] b = new double[maxit_];
		
		hold = 0;//TMP
			    
		if ((lambda_.length<1) || (q_<=0) || (maxit_ <1) || (eps_<=0.0)) {
			res_ = -2.0;
			ifault_ = 2;
			return;
		}

		tol = -200.0;

		// Preliminaries
		sum = lambda_[0];
		bbeta = sum;
			      
		for (i=1;i<=lambda_.length;i++) {
			hold = lambda_[i-1];
			if ((hold<=0.0) || (h_[i-1]<1) || (delta_[i-1]<0.0)) {
				res_ = -7.0;
				ifault_ = -i;
				return;
			}	
			if (bbeta > hold) bbeta = hold; // calcul du max des lambdas
			if (sum < hold) sum = hold;    // calcul du min des lambdas
		}
			  
		if (mode_ > 0.0) {
			// if ((2.0/(1.0/bbeta+1.0/sum))>1.8*sum) bbeta = sum; // comme dans NAG : methode avec betaA
			bbeta = mode_*bbeta;
		} else {
			bbeta = 2.0/(1.0/bbeta+1.0/sum);  // methode avec betaB
		}

		k = 0;
		sum = 1.0;
		sum1 = 0.0;
		for (i=1;i<=lambda_.length;i++) {
			hold = bbeta/lambda_[i-1];
			gamma[i-1] = 1.0 - hold;
			sum = sum*Math.pow(hold,h_[i-1]); //???? pas sur ..
			sum1 = sum1 + delta_[i-1];
			k = k + h_[i-1];
			theta[i-1] = 1.0;
		}
		
//		System.out.println("sum: " + sum);
//		System.out.println("sum1: " +  sum1);
//		System.out.println("hold: " + hold);
			    
		ao = Math.exp(0.5*(Math.log(sum)-sum1));
		
//		System.out.println("ao:" + ao);
		
		if (ao <= 0.0) {
			res_ = 0.0;
			dnsty = 0.0;
			ifault_ = 1; // returns after this (the rest of the function is in the else)
		} else { // evaluate probability and density of chi-squared on k degrees of freedom. The constant 0.22579135264473 is ln(sqrt(pi/2))
			z = q_/bbeta;
			      
			if ((k%2) == 0) { // k est un entier donc on regarde si k est divisible par 2: k == (k/2)*k 
				i = 2;
				lans = -0.5*z;
				dans = Math.exp(lans);
				pans = 1.0 - dans;
			} else {
				i = 1;
				lans = -0.5*(z+Math.log(z)) - 0.22579135264473;
				dans = Math.exp(lans);
				//pans = pnorm(Math.sqrt(z),0.0,1.0,1,0) - pnorm(-Math.sqrt(z),0.0,1.0,1,0); 
				pans = normal_.cumulativeProbability(Math.sqrt(z)) - normal_.cumulativeProbability(-Math.sqrt(z)); 
			}
			      
			k = k-2;
			for (j=i;j<=k;j=j+2) {
				if (lans < tol) {
					lans = lans + Math.log(z/(double)j);
					dans = Math.exp(lans);
				} else {
					dans = dans*z/(double)j;
				}
				pans = pans -dans;
			}
			      
			// evaluate successive terms of expansion
			      
			prbty = pans;
			dnsty = dans;
			eps2 = eps_/ao;
			aoinv = 1.0/ao;
			sum = aoinv - 1.0;

			for (m=1;m<=maxit_;m++) {
				sum1 = 0.0;
				for (i=1;i<=lambda_.length;i++) {
					hold = theta[i-1];
					hold2 = hold*gamma[i-1];
					theta[i-1] = hold2;
					sum1 = sum1 + hold2*h_[i-1]+m*delta_[i-1]*(hold-hold2);
				}
				sum1 = 0.5*sum1;
				b[m-1] = sum1;
				for (i=m-1;i>=1;i--) {
					sum1 = sum1 + b[i-1]*a[m-i-1]; 
				}
				sum1 = sum1/(double)m;
				a[m-1] = sum1;
				k = k + 2;
				if (lans < tol) {
					lans = lans + Math.log(z/(double)k);
					dans = Math.exp(lans);
				} else {
					dans = dans*z/(double)k;
				}
				pans = pans - dans;
				sum = sum - sum1;
				dnsty = dnsty + dans*sum1;
				sum1 = pans*sum1;
				prbty = prbty + sum1;
				if (prbty<(-aoinv)) {
					res_ = -3.0;
					ifault_ = 3;
					return;
				}
				
//				if(m == 1 || m % 1000 == 0){
//					System.out.println(sum + "\t" + sum1 + "\t" + lans + "\t"  + dans + "\t" + pans + "\t" + dnsty  + "\t" + prbty);
//				}
				
				if (Math.abs(pans*sum) < eps2) {
					if (Math.abs(sum1) < eps2) {
						ifault_ = 0; // this is overwritten below (ifault_=4) -- I now changed the code below
						// I COMMENTED THIS SO WE CAN See IF THERE WAS CONVERGENCE OR NOT -daniel
						//m = maxit_+1; // and WHY would you do that?
						break;
					}
				}
				
			}

			if (m > maxit_)  // I ADDED THIS IF, OTHERWISE IT MAKES NO SENSE -daniel
				ifault_ = 4;
			// Check if I understood correctly
			if(!(ifault_ == 0 || ifault_ == 4)){
				throw new RuntimeException();
			}
			
//			System.out.println("ao : " + ao);
//			System.out.println("prbty: " + prbty);
//			
			dnsty = ao*dnsty/(bbeta+bbeta);
			prbty = ao*prbty;
			
			// With my edits above, this now makes a bit mores sense
			if (prbty<0.0 || prbty>1.0) {ifault_ = ifault_ + 5; // why the ... would they write it like this? I.e., ifault_ = 9
			} else {
				if (dnsty<0.0) ifault_ = ifault_ + 6;
			}
			res_ = prbty;
		}

		return;
	}
	
	
	// ============================================================================
	// GETTERS AND SETTERS
	
	public int getIfault() { return ifault_; }

}


//	Original code from the R-package
//  ================================

//farebrother <- function(q,lambda,h = rep(1,length(lambda)),delta = rep(0,length(lambda)),maxit=100000,eps=10^(-10),mode=1) {
//
//
//	  r <- length(lambda)
//	  if (length(h) != r) stop("lambda and h should have the same length!")
//	  if (length(delta) != r) stop("lambda and delta should have the same length!")
//
//	#  mode <- 1 # ??
//
//	  dnsty <- 0.0
//	  ifault <- 0
//	  res <- 0 
//	  
//
//	  out <- .C("ruben",lambda=as.double(lambda),h=as.integer(h),delta=as.double(delta),r=as.integer(r),q=as.double(q),mode=as.double(mode),maxit=as.integer(maxit),eps=as.double(eps),dnsty=as.double(dnsty),ifault=as.integer(ifault),res=as.double(res),PACKAGE="CompQuadForm")
//
//	  out$res <- 1 - out$res
//
//	  return(out)
//
//	}
//
//
//
//#include <R.h>
//#include "Rmath.h"
//
//extern "C" {
//
////Algorithm AS 204 Appl. Statist. (1984) Vol. 33, No.3
////ruben evaluates the probability that a positive definite quadratic form in Normal variates is less than a given value
// 
// void ruben(double *lambda, int *mult, double *delta, int *n, double *c, double *mode, int *maxit, double *eps, double *dnsty, int *ifault, double *res) {
//   
//   int i,k,m,j;
//   double pnorm(double q, double mean, double sd, int lower_tail, int log_p);
//   double ao, aoinv, z, bbeta, eps2, hold, hold2, sum, sum1, dans, lans, pans, prbty, tol;
//   double *gamma, *theta, *a, *b;
//   gamma = new double[n[0]];
//   theta = new double[n[0]];
//   a = new double[maxit[0]];
//   b = new double[maxit[0]];
//   
//   
//   
//   if ((n[0]<1) || (c[0]<=0) || (maxit[0] <1) || (eps[0]<=0.0)) {
//     res[0] = -2.0;
//     ifault[0] = 2;
//     delete[] gamma;
//     delete[] theta;
//     delete[] a;
//     delete[] b;
//     return;
//   } else {
//     tol = -200.0;
//   
//     // Preliminaries
//     sum = lambda[0];
//     bbeta = sum;
//     
//     for (i=1;i<=n[0];i++) {
//	hold = lambda[i-1];
//	if ((hold<=0.0) || (mult[i-1]<1) || (delta[i-1]<0.0)) {
//	  res[0] = -7.0;
//	  ifault[0] = -i;
//	  delete[] gamma;
//	  delete[] theta;
//	  delete[] a;
//	  delete[] b;
//	  return;
//	}	
//	if (bbeta > hold) bbeta = hold; // calcul du max des lambdas
//	if (sum < hold) sum = hold;    // calcul du min des lambdas
//     }
//     
// 
//   if (mode[0] > 0.0) {
//     // if ((2.0/(1.0/bbeta+1.0/sum))>1.8*sum) bbeta = sum; // comme dans NAG : methode avec betaA
//           bbeta = mode[0]*bbeta;
//   } else {
//     bbeta = 2.0/(1.0/bbeta+1.0/sum);  // methode avec betaB
//   }
//
//   k = 0;
//   sum = 1.0;
//   sum1 = 0.0;
//   for (i=1;i<=n[0];i++) {
//     hold = bbeta/lambda[i-1];
//     gamma[i-1] = 1.0 - hold;
//     sum = sum*R_pow(hold,mult[i-1]); //???? pas sur ..
//     sum1 = sum1 + delta[i-1];
//     k = k + mult[i-1];
//     theta[i-1] = 1.0;
//   }
//   
//   ao = exp(0.5*(log(sum)-sum1));
//   if (ao <= 0.0) {
//     res[0] = 0.0;
//     dnsty[0] = 0.0;
//     ifault[0] = 1;
//   } else { // evaluate probability and density of chi-squared on k degrees of freedom. The constant 0.22579135264473 is ln(sqrt(pi/2))
//     z = c[0]/bbeta;
//     
//     if ((k%2) == 0) { // k est un entier donc on regarde si k est divisible par 2: k == (k/2)*k 
//	i = 2;
//	lans = -0.5*z;
//	dans = exp(lans);
//	pans = 1.0 - dans;
//     } else {
//	i = 1;
//	lans = -0.5*(z+log(z)) - 0.22579135264473;
//	dans = exp(lans);
//	pans = pnorm(sqrt(z),0.0,1.0,1,0) - pnorm(-sqrt(z),0.0,1.0,1,0); 
//     }
//     
//     k = k-2;
//     for (j=i;j<=k;j=j+2) {
//	if (lans < tol) {
//	  lans = lans + log(z/(double)j);
//	  dans = exp(lans);
//	} else {
//	  dans = dans*z/(double)j;
//	}
//	pans = pans -dans;
//     }
//     
//     // evaluate successive terms of expansion
//     
//     prbty = pans;
//     dnsty[0] = dans;
//     eps2 = eps[0]/ao;
//     aoinv = 1.0/ao;
//     sum = aoinv - 1.0;
//   
//
//   for (m=1;m<=maxit[0];m++) {
//     sum1 = 0.0;
//     for (i=1;i<=n[0];i++) {
//	hold = theta[i-1];
//	hold2 = hold*gamma[i-1];
//	theta[i-1] = hold2;
//	sum1 = sum1 + hold2*mult[i-1]+m*delta[i-1]*(hold-hold2);
//     }
//     sum1 = 0.5*sum1;
//     b[m-1] = sum1;
//     for (i=m-1;i>=1;i--) {
//	sum1 = sum1 + b[i-1]*a[m-i-1]; 
//     }
//     sum1 = sum1/(double)m;
//     a[m-1] = sum1;
//     k = k + 2;
//     if (lans < tol) {
//	lans = lans + log(z/(double)k);
//	dans = exp(lans);
//     } else {
//	dans = dans*z/(double)k;
//     }
//     pans = pans - dans;
//     sum = sum - sum1;
//     dnsty[0] = dnsty[0] + dans*sum1;
//     sum1 = pans*sum1;
//     prbty = prbty + sum1;
//     if (prbty<(-aoinv)) {
//	res[0] = -3.0;
//	ifault[0] = 3;
//	return;
//     }
//     if (fabs(pans*sum) < eps2) {
//	if (fabs(sum1) < eps2) {
//	  ifault[0] = 0;
//	  
//	  m = maxit[0]+1;
//	  break;
//
//	}
//     }
//   }
//
//   ifault[0] = 4;
//   dnsty[0] = ao*dnsty[0]/(bbeta+bbeta);
//   prbty = ao*prbty;
//   if (prbty<0.0 || prbty>1.0) {ifault[0] = ifault[0] + 5;
//   } else {
//     if (dnsty[0]<0.0) ifault[0] = ifault[0] + 6;
//   }
//   res[0] = prbty;
//   }
//
//   delete[] gamma;
//   delete[] theta;
//   delete[] a;
//   delete[] b;
//   return;
//   }
//   
// }
//}



			


