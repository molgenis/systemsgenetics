//******************************************************************************
//
// File:    NonNegativeLeastSquares.java
// Package: edu.rit.numeric
// Unit:    Class edu.rit.numeric.NonNegativeLeastSquares
//
// This Java source file is copyright (C) 2005 by Alan Kaminsky. All rights
// reserved. For further information, contact the author, Alan Kaminsky, at
// ark@cs.rit.edu.
//
// This Java source file is part of the Parallel Java Library ("PJ"). PJ is free
// software; you can redistribute it and/or modify it under the terms of the GNU
// General Public License as published by the Free Software Foundation; either
// version 3 of the License, or (at your option) any later version.
//
// PJ is distributed in the hope that it will be useful, but WITHOUT ANY
// WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR
// A PARTICULAR PURPOSE. See the GNU General Public License for more details.
//
// Linking this library statically or dynamically with other modules is making a
// combined work based on this library. Thus, the terms and conditions of the
// GNU General Public License cover the whole combination.
//
// As a special exception, the copyright holders of this library give you
// permission to link this library with independent modules to produce an
// executable, regardless of the license terms of these independent modules, and
// to copy and distribute the resulting executable under terms of your choice,
// provided that you also meet, for each linked independent module, the terms
// and conditions of the license of that module. An independent module is a
// module which is not derived from or based on this library. If you modify this
// library, you may extend this exception to your version of the library, but
// you are not obligated to do so. If you do not wish to do so, delete this
// exception statement from your version.
//
// A copy of the GNU General Public License is provided in the file gpl.txt. You
// may also obtain a copy of the GNU General Public License on the World Wide
// Web at http://www.gnu.org/licenses/gpl.html.
//
//******************************************************************************

package main.java.decon_eQTL;

import org.apache.commons.math3.exception.MathIllegalArgumentException;

/**
 * 
 * This code is copied from the edu.rit.numeric package. Downloaded the source
 * from: https://www.cs.rit.edu/~ark/pj.shtml#download
 * Manual: https://www.cs.rit.edu/~ark/pj/doc/edu/rit/numeric/NonNegativeLeastSquares.html
 * 
 * I (Niek de Klein, 2017) made some adjustments for inputting data to solve() method and outputting some other data.
 * 
 * Also extended from math.commons AbstractMultipleLinearRegression
 * 
 * Most important change: original code replaces the A matrix and b vector given as input
 * to the solve() method by their orthogonals (see solve() documentation). However, here
 * they get cloned first, so that the input data given to solve() keeps its original values.
 * This is because I do not use the orthogonal data outside of this class, and I do want to
 * keep my original values.
 * 
 * Class NonNegativeLeastSquares provides a method for solving a least squares
 * minimization problem with nonnegativity constraints. The solve()
 * method finds an approximate solution to the linear system of equations
 * Ax = b, such that
 * ||Ax&nbsp;-&nbsp;b||<SUP>2</SUP> is minimized, and such that
 * x &gt;= 0. The inputs to and outputs from the solve()
 * method are stored in the fields of an instance of class
 * NonNegativeLeastSquares.
 * <P>
 * The Java code is a translation of the Fortran subroutine NNLS from
 * Charles L. Lawson and Richard J. Hanson, Solving Least Squares
 * Problems (Society for Industrial and Applied Mathematics, 1995), page
 * 161.
 *
 * @author  Alan Kaminsky (modified by Niek de Klein)
 * @version 14-Dec-2017
 */
public class NonNegativeLeastSquares
{

	/**
	 * The number of rows, typically the number of input data points, in the
	 * least squares problem.
	 */
	private int M;

	/**
	 * The number of columns, typically the number of output parameters, in the
	 * least squares problem.
	 */
	private int N;

	
	// initial y
	private double[] measuredValues;
	// initial x
	private double[][] observedValues;
	// predicted values for y
	private double[] predictedValues;
	/**
	 * The N-element x vector for the least squares problem. On
	 * output from the solve() method, x contains the solution
	 * vector x.
	 */
	private double[] x;

	/**
	 * The N-element index vector. On output from the solve()
	 * method: index[0] through index[nsetp-1] contain the
	 * indexes of the elements in x that are in set P, the set of
	 * positive values; that is, the elements that are not forced to be zero
	 * (inactive constraints). index[nsetp] through index[N-1]
	 * contain the indexes of the elements in x that are in set Z,
	 * the set of zero values; that is, the elements that are forced to be zero
	 * (active constraints).
	 */
	private int[] index;

	/**
	 * The number of elements in the set P; that is, the number of
	 * positive values (inactive constraints). An output of the solve()
	 * method.
	 */
	private int nsetp;


	// After solving, the orthogonal matrix of the A matrix and b vector
	// cloning originalA and originalB so that those values are kept for later use
	private double[] b;
	private double[][] a;

	// Working storage.
	private double[] w;
	private double[] zz;
	private double[] terms;

	// Maximum number of iterations.
	private int itmax;

	// Magic numbers.
	private static final double factor = 0.01;

	// Exported constructors.

	/**
	 * Construct a new nonnegative least squares problem of the given size.
	 * Fields M and N are set based on the observed values matrix. The array
	 * fields a, b, x, and index are
	 * allocated with the proper sizes but are not filled in.
	 *
	 * @exception  IllegalArgumentException
	 *     (unchecked exception) Thrown if M &lt;= 0 or N
	 *     &lt;= 0.
	 */
	public NonNegativeLeastSquares(){}

	/**
     * Loads model x and y sample data, overriding any previous sample.
     *
     * @param y the [n,1] array representing the y sample
     * @param x the [n,k] array representing the x sample
     * @throws MathIllegalArgumentException if the x and y array data are not
     *             compatible for the regression
     */
    public void newSampleData(double[] y, double[][] x) throws MathIllegalArgumentException {
    	this.measuredValues = y;
    	this.observedValues = x;
		// cloning y and x so that those values are kept for later use
    	try{
    		b = y.clone();
		}
		catch (NullPointerException e){
			DeconvolutionLogger.log.info("ERROR: Expression values are not read in correctly, check if input files are correct.");
			throw(e);
		}
		a = new double[x.length][];
		for(int z = 0; z < x.length; ++z)
			a[z] = x[z].clone();

		//Number of rows (input data points) in the least squares problem.
		int M = x.length;
		// Number of columns (output parameters) in the least squares problem.
		int N = x[0].length;
		if (M <= 0)
		{
			throw new IllegalArgumentException
			("NonNegativeLeastSquares(): M = " + M + " illegal");
		}
		if (N <= 0)
		{
			throw new IllegalArgumentException
			("NonNegativeLeastSquares(): N = " + N + " illegal");
		}

		this.M = M;
		this.N = N;
		//this.a = new double [M] [N];
		//this.b = new double [M];
		this.x = new double [N];
		this.index = new int [N];

		this.w = new double [N];
		this.zz = new double [M];
		this.terms = new double [2];
		this.itmax = 3*N;
		
		solve();
    }
	
	/**
	 * Solve this least squares minimization problem with nonnegativity
	 * constraints. The solve() method finds an approximate solution to
	 * the linear system of equations Ax = b, such that
	 * ||Ax&nbsp;-&nbsp;b||<SUP>2</SUP> is minimized, and such
	 * that x &gt;= 0. On input, the field a must be
	 * filled in with the matrix A and the field b must be
	 * filled in with the vector b for the problem to be solved. On
	 * output, the other fields are filled in with the solution as explained in
	 * the documentation for each field.
	 * 
	 * The MxN-element A matrix for the least squares
	 * problem. On input to the solve() method, originalA contains the
	 * matrix A. On output, originalA is cloned into a. a has been replaced with QA,
	 * where Q is an MxM-element orthogonal matrix
	 * generated during the solve() method's execution.
	 *
	 * The M-element b vector for the least squares problem. On
	 * input to the solve() method, originalB contains the vector
	 * b. On output, originalB is cloned into b. b has been replaced with Qb, where
	 * Q is an MxM-element orthogonal matrix generated
	 * during the solve() method's execution.
	 *
	 * @exception  TooManyIterationsException
	 *     (unchecked exception) Thrown if too many iterations occurred without
	 *     finding a minimum (more than 3N iterations).
	 */
	private void solve()
	{
		int i, iz, j, l, izmax, jz, jj, ip, ii;
		double sm, wmax, asave, unorm, ztest, up, alpha, t, cc, ss, temp;

		// Keep count of iterations.
		int iter = 0;

		// Initialise the arrays index and x.
		// index[0] through index[nsetp-1] = set P.
		// index[nsetp] through index[N-1] = set Z.
		for (i = 0; i < N; ++ i)
		{
			x[i] = 0.0;
			index[i] = i;
		}
		nsetp = 0;

		// Main loop begins here.
		mainloop: for (;;)
		{
			// Quit if all coefficients are already in the solution, or if M
			// columns of A have been triangularized.
			if (nsetp >= N || nsetp >= M) break mainloop;

			// Compute components of the dual (negative gradient) vector W.
			for (iz = nsetp; iz < N; ++ iz)
			{
				j = index[iz];
				sm = 0.0;
				for (l = nsetp; l < M; ++ l)
				{
					sm += a[l][j]*b[l];

				}
				w[j] = sm;
			}

			// Find a candidate j to be moved from set Z to set P.
			candidateloop: for (;;)
			{
				// Find largest positive W[j].
				wmax = 0.0;
				izmax = -1;
				for (iz = nsetp; iz < N; ++ iz)
				{
					j = index[iz];
					if (w[j] > wmax)
					{
						wmax = w[j];
						izmax = iz;
					}
				}

				// If wmax <= 0, terminate. This indicates satisfaction of the
				// Kuhn-Tucker conditions.
				if (wmax <= 0.0) break mainloop;
				iz = izmax;
				j = index[iz];

				// The sign of W[j] is okay for j to be moved to set P. Begin
				// the transformation and check new diagonal element to avoid
				// near linear independence.
				asave = a[nsetp][j];
				up = constructHouseholderTransform (nsetp, nsetp+1, a, j);
				unorm = 0.0;
				for (l = 0; l < nsetp; ++ l)
				{
					unorm += sqr (a[l][j]);
				}
				unorm = Math.sqrt (unorm);
				if (diff (unorm + Math.abs(a[nsetp][j])*factor, unorm) > 0.0)
				{
					// Column j is sufficiently independent. Copy B into ZZ,
					// update ZZ, and solve for ztest = proposed new value for
					// X[j].
					System.arraycopy (b, 0, zz, 0, M);
					applyHouseholderTransform (nsetp, nsetp+1, a, j, up, zz);
					ztest = zz[nsetp] / a[nsetp][j];

					// If ztest is positive, we've found our candidate.
					if (ztest > 0.0) break candidateloop;
				}

				// Reject j as a candidate to be moved from set Z to set P.
				// Restore a[nsetp][j], set w[j] = 0, and try again.
				a[nsetp][j] = asave;
				w[j] = 0.0;
			}

			// The index j = index[iz] has been selected to be moved from set Z
			// to set P. Update B, update indexes, apply Householder
			// transformations to columns in new set Z, zero subdiagonal
			// elements in column j, set w[j] = 0.
			System.arraycopy (zz, 0, b, 0, M);

			index[iz] = index[nsetp];
			index[nsetp] = j;
			++ nsetp;

			jj = -1;
			for (jz = nsetp; jz < N; ++ jz)
			{
				jj = index[jz];
				applyHouseholderTransform (nsetp-1, nsetp, a, j, up, a, jj);
			}

			for (l = nsetp; l < M; ++ l)
			{
				a[l][j] = 0.0;
			}

			w[j] = 0.0;

			// Solve the triangular system. Store the solution temporarily in
			// zz.
			for (l = 0; l < nsetp; ++ l)
			{
				ip = nsetp - l;
				if (l != 0)
				{
					for (ii = 0; ii < ip; ++ ii)
					{
						zz[ii] -= a[ii][jj] * zz[ip];
					}
				}
				-- ip;
				jj = index[ip];
				zz[ip] /= a[ip][jj];
			}

			// Secondary loop begins here.
			secondaryloop: for (;;)
			{
				// Increment iteration counter.
				++ iter;
				if (iter > itmax)
				{
					throw new RuntimeException
					("NonNegativeLeastSquares.solve(): Too many iterations");
				}

				// See if all new constrained coefficients are feasible. If not,
				// compute alpha.
				alpha = 2.0;
				for (ip = 0; ip < nsetp; ++ ip)
				{
					l = index[ip];
					if (zz[ip] <= 0.0)
					{
						t = -x[l] / (zz[ip] - x[l]);
						if (alpha > t)
						{
							alpha = t;
							jj = ip;
						}
					}
				}

				// If all new constrained coefficients are feasible then alpha
				// will still be 2. If so, exit from secondary loop to main
				// loop.
				if (Math.abs(alpha- 2.0) < 0.000000000000000000000000001) break secondaryloop;

				// Otherwise, use alpha (which will be between 0 and 1) to
				// interpolate between the old x and the new zz.
				for (ip = 0; ip < nsetp; ++ ip)
				{
					l = index[ip];
					x[l] += alpha * (zz[ip] - x[l]);
				}

				// Modify A and B and the index arrays to move coefficient i
				// from set P to set Z.
				i = index[jj];
				tertiaryloop: for (;;)
				{
					x[i] = 0.0;
					if (jj != nsetp-1)
					{
						++ jj;
						for (j = jj; j < nsetp; ++ j)
						{
							ii = index[j];
							index[j-1] = ii;
							a[j-1][ii] =
									computeGivensRotation
									(a[j-1][ii], a[j][ii], terms);
							a[j][ii] = 0.0;
							cc = terms[0];
							ss = terms[1];
							for (l = 0; l < N; ++ l)
							{
								if (l != ii)
								{
									// Apply Givens rotation to column l of A.
									temp = a[j-1][l];
									a[j-1][l] =  cc*temp + ss*a[j][l];
									a[j  ][l] = -ss*temp + cc*a[j][l];
								}
							}
							// Apply Givens rotation to B.
							temp = b[j-1];
							b[j-1] =  cc*temp + ss*b[j];
							b[j  ] = -ss*temp + cc*b[j];
						}
					}
					-- nsetp;
					index[nsetp] = i;

					// See if the remaining coefficients in set P are feasible.
					// They should be because of the way alpha was determined.
					// If any are infeasible it is due to roundoff error. Any
					// that are nonpositive will be set to 0 and moved from set
					// P to set Z.
					for (jj = 0; jj < nsetp; ++ jj)
					{
						i = index[jj];
						if (x[i] <= 0.0) continue tertiaryloop;
					}
					break tertiaryloop;
				}

				// Copy b into zz, then solve the tridiagonal system again and
				// continue the secondary loop.
				System.arraycopy (b, 0, zz, 0, M);
				for (l = 0; l < nsetp; ++ l)
				{
					ip = nsetp - l;
					if (l != 0)
					{
						for (ii = 0; ii < ip; ++ ii)
						{
							zz[ii] -= a[ii][jj] * zz[ip];
						}
					}
					-- ip;
					jj = index[ip];
					zz[ip] /= a[ip][jj];
				}
			}

			// Update x from zz.
			for (ip = 0; ip < nsetp; ++ ip)
			{
				i = index[ip];
				x[i] = zz[ip];
			}

			// All new coefficients are positive. Continue the main loop.
		}
	}

	// Hidden operations.

	/**
	 * Construct a Householder transformation. u is an
	 * MxN-element matrix used as an input and an output of this
	 * method.
	 *
	 * @param  ipivot
	 *     Index of the pivot element within the pivot vector.
	 * @param  i1
	 *     If i1 &lt; M, the transformation will be constructed
	 *     to zero elements indexed from i1 through M-1. If
	 *     i1 &gt;= M, an identity transformation will be
	 *     constructed.
	 * @param  u
	 *     An MxN-element matrix. On input, column
	 *     pivotcol of u contains the pivot vector. On output,
	 *     column pivotcol of u, along with the return value
	 *     (up), contains the Householder transformation.
	 * @param  pivotcol
	 *     Index of the column of u that contains the pivot vector.
	 *
	 * @return
	 *     The quantity up which is part of the Householder
	 *     transformation.
	 */
	private static double constructHouseholderTransform
	(int ipivot,
			int i1,
			double[][] u,
			int pivotcol)
	{
		int M = u.length;
		int j;
		double cl, clinv, sm, up;

		cl = Math.abs (u[ipivot][pivotcol]);

		// Construct the transformation.
		for (j = i1; j < M; ++ j)
		{
			cl = Math.max (Math.abs (u[j][pivotcol]), cl);
		}
		if (cl <= 0.0)
		{
			throw new IllegalArgumentException
			("NonNegativeLeastSquares.constructHouseholderTransform(): Illegal pivot vector");
		}
		clinv = 1.0 / cl;
		sm = sqr (u[ipivot][pivotcol] * clinv);
		for (j = i1; j < M; ++ j)
		{
			sm += sqr (u[j][pivotcol] * clinv);
		}
		cl = cl * Math.sqrt (sm);
		if (u[ipivot][pivotcol] > 0.0) cl = -cl;
		up = u[ipivot][pivotcol] - cl;
		u[ipivot][pivotcol] = cl;
		return up;
	}

	/**
	 * Apply a Householder transformation to one column of a matrix. u
	 * is an MxN-element matrix used as an input of this method.
	 * c is an MxN-element matrix used as an input and
	 * output of this method. ipivot, i1, u, and
	 * pivotcol must be the same as in a previous call of
	 * constructHouseholderTransform(), and up must be the
	 * value returned by that method call.
	 *
	 * @param  ipivot
	 *     Index of the pivot element within the pivot vector.
	 * @param  i1
	 *     If i1 &lt; M, the transformation will zero elements
	 *     indexed from i1 through M-1. If i1 &gt;=
	 *     M, the transformation is an identity transformation.
	 * @param  u
	 *     An MxN-element matrix. On input, column
	 *     pivotcol of u, along with up, contains the
	 *     Householder transformation. This must be the output of a previous
	 *     call of constructHouseholderTransform().
	 * @param  pivotcol
	 *     Index of the column of u that contains the Householder
	 *     transformation.
	 * @param  up
	 *     The rest of the Householder transformation. This must be the return
	 *     value of the same previous call of
	 *     constructHouseholderTransform().
	 * @param  c
	 *     An MxN-element matrix. On input, column
	 *     applycol of c contains the vector to which the
	 *     Householder transformation is to be applied. On output, column
	 *     applycol of c contains the transformed vector.
	 * @param  applycol
	 *     Index of the column of c to which the Householder
	 *     transformation is to be applied.
	 */
	private static void applyHouseholderTransform
	(int ipivot,
			int i1,
			double[][] u,
			int pivotcol,
			double up,
			double[][] c,
			int applycol)
	{
		int M = u.length;
		int i;
		double cl, b, sm;

		cl = Math.abs (u[ipivot][pivotcol]);
		if (cl <= 0.0)
		{
			throw new IllegalArgumentException
			("NonNegativeLeastSquares.applyHouseholderTransform(): Illegal pivot vector");
		}

		b = up * u[ipivot][pivotcol];
		// b must be nonpositive here. If b = 0, return.
		if (b == 0.0)
		{
			return;
		}
		else if (b > 0.0)
		{
			throw new IllegalArgumentException
			("NonNegativeLeastSquares.applyHouseholderTransform(): Illegal pivot element");
		}
		b = 1.0 / b;

		sm = c[ipivot][applycol] * up;
		for (i = i1; i < M; ++ i)
		{
			sm += c[i][applycol] * u[i][pivotcol];
		}
		if (sm != 0.0)
		{
			sm = sm * b;
			c[ipivot][applycol] += sm * up;
			for (i = i1; i < M; ++ i)
			{
				c[i][applycol] += sm * u[i][pivotcol];
			}
		}
	}

	/**
	 * Apply a Householder transformation to a vector. u is an
	 * MxN-element matrix used as an input of this method.
	 * c is an M-element array used as an input and output of
	 * this method. ipivot, i1, u, and
	 * pivotcol must be the same as in a previous call of
	 * constructHouseholderTransform(), and up must be the
	 * value returned by that method call.
	 *
	 * @param  ipivot
	 *     Index of the pivot element within the pivot vector.
	 * @param  i1
	 *     If i1 &lt; M, the transformation will zero elements
	 *     indexed from i1 through M-1. If i1 &gt;=
	 *     M, the transformation is an identity transformation.
	 * @param  u
	 *     An MxN-element matrix. On input, column
	 *     pivotcol of u, along with up, contains the
	 *     Householder transformation. This must be the output of a previous
	 *     call of constructHouseholderTransform().
	 * @param  pivotcol
	 *     Index of the column of u that contains the Householder
	 *     transformation.
	 * @param  up
	 *     The rest of the Householder transformation. This must be the return
	 *     value of the same previous call of
	 *     constructHouseholderTransform().
	 * @param  c
	 *     An M-element array. On input, c contains the vector
	 *     to which the Householder transformation is to be applied. On output,
	 *     c contains the transformed vector.
	 */
	private static void applyHouseholderTransform
	(int ipivot,
			int i1,
			double[][] u,
			int pivotcol,
			double up,
			double[] c)
	{
		int M = u.length;
		int i;
		double cl, b, sm;

		cl = Math.abs (u[ipivot][pivotcol]);
		if (cl <= 0.0)
		{
			throw new IllegalArgumentException
			("NonNegativeLeastSquares.applyHouseholderTransform(): Illegal pivot vector");
		}

		b = up * u[ipivot][pivotcol];
		// b must be nonpositive here. If b = 0, return.
		if (b == 0.0)
		{
			return;
		}
		else if (b > 0.0)
		{
			throw new IllegalArgumentException
			("NonNegativeLeastSquares.applyHouseholderTransform(): Illegal pivot element");
		}
		b = 1.0 / b;

		sm = c[ipivot] * up;
		for (i = i1; i < M; ++ i)
		{
			sm += c[i] * u[i][pivotcol];
		}
		if (sm != 0.0)
		{
			sm = sm * b;
			c[ipivot] += sm * up;
			for (i = i1; i < M; ++ i)
			{
				c[i] += sm * u[i][pivotcol];
			}
		}
	}

	/**
	 * Compute the sine and cosine terms of a Givens rotation matrix. The terms
	 * c and s are returned in terms[0] and
	 * terms[1], respectively, such that:
	 * <PRE>
	 *     [ c  s] * [a] = [sqrt(a^2+b^2)]
	 *     [-s  c]   [b]   [      0      ]
	 * </PRE>
	 *
	 * @param  a      Input argument.
	 * @param  b      Input argument.
	 * @param  terms  A 2-element array. On output, terms[0] contains
	 *                c and terms[1] contains s.
	 *
	 * @return  sqrt(a<SUP>2</SUP>+b<SUP>2</SUP>).
	 */
	private static double computeGivensRotation
	(double a,
			double b,
			double[] terms)
	{
		double xr, yr;

		if (Math.abs(a) > Math.abs(b))
		{
			xr = b/a;
			yr = Math.sqrt (1.0 + sqr (xr));
			terms[0] = sign (1.0/yr, a);
			terms[1] = terms[0]*xr;
			return Math.abs(a)*yr;
		}
		else if (b != 0.0)
		{
			xr = a/b;
			yr = Math.sqrt (1.0 + sqr (xr));
			terms[1] = sign (1.0/yr, b);
			terms[0] = terms[1]*xr;
			return Math.abs(b)*yr;
		}
		else
		{
			terms[0] = 0.0;
			terms[1] = 1.0;
			return 0.0;
		}
	}

	/**
	 * Determine if x differs from y, to machine precision.
	 *
	 * @return  0.0, if x is the same as y to machine precision; x-y (nonzero),
	 *          if x differs from y to machine precision.
	 */
	private static double diff
	(double x,
			double y)
	{
		return x - y;
	}

	/**
	 * Returns x^2.
	 */
	private static double sqr
	(double x)
	{
		return x*x;
	}

	/**
	 * Returns the number whose absolute value is x and whose sign is the same
	 * as that of y. x is assumed to be nonnegative.
	 */
	private static double sign
	(double x,
			double y)
	{
		return y >= 0.0 ? x : -x;
	}

	protected double calculateResidualSumOfSquares(){
		/**
		 * The squared Euclidean norm of the residual vector, ||Ax -
		 * b||<SUP>2</SUP>. An output of the solve() method.
		 */
		// Compute the squared Euclidean norm of the final residual vector.
		double normsqr = 0.0;
		for (int i = nsetp; i < M; ++ i)
		{	
			normsqr += sqr (b[i]);
		}

		return normsqr;
	}

	public double[] estimateRegressionParameters() {
		return this.x;
	}
	
	public double[] getPredictedExpressionValues(){
		if(predictedValues != null){
			return predictedValues;
		}
		predictedValues = new double[measuredValues.length];
		for(int i = 0; i < measuredValues.length; ++i ){
			double predictedValue = 0;
			for(int z = 0; z < x.length; ++z){
				predictedValue += x[z] * observedValues[i][z];
			}
			predictedValues[i] = predictedValue;
		}
		return predictedValues;
	}
	
	public double[] estimateResiduals() {
		double[] predictedValues = getPredictedExpressionValues();
		double[] residuals = new double[measuredValues.length];
		for(int i = 0; i < measuredValues.length; ++i ){
			residuals[i] = measuredValues[i] - predictedValues[i];
		}
		return(residuals);
	}

}
