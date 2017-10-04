package deconvolution;


import java.util.ArrayList;
import java.util.List;

import Jama.Matrix;
	/*
	 * Copyright 2008 Josh Vermaas, except he's nice and instead prefers
	 * this to be licensed under the LGPL. Since the license itself is longer
	 * than the code, if this truly worries you, you can look up the text at
	 * http://www.gnu.org/licenses/
	 * 
	 * Now onto the required "I didn't do it" part.
	 * This program is free software: you can redistribute it and/or modify
	    it under the terms of the GNU Lesser General Public License as published by
	    the Free Software Foundation, either version 3 of the License, or
	    (at your option) any later version.

	    This program is distributed in the hope that it will be useful,
	    but WITHOUT ANY WARRANTY; without even the implied warranty of
	    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
	    GNU General Public License for more details.

	    You should have received a copy of the GNU Lesser General Public License
	    along with this program.  If not, see <http://www.gnu.org/licenses/>.

	 */
	/**
	 * This is a perfectly functional Java NNLS solver that utilizes JAMA for all
	 * the Matrix calculations. As such, one will need Jama (currently at:
	 * http://math.nist.gov/javanumerics/jama/) to use this program.
	 * @author Josh Vermaas
	 *
	 * I got it from here: http://odell.radonc.med.ufl.edu/WOImageJ/woSource/NNLSSolver.html
	 *
	 */
	public class NonNegativeLeastSquaresJAMA {
		/**
		 * This is my own Java implementation of the NNLS algorithm
		 * as described in:
		 * Lawson and Hanson, "Solving Least Squares Problems", Prentice-Hall, 1974, Chapter 23, p. 161.
		 * It not only solves the least squares problem, but does so while also requiring
		 * that none of the answers be negative.
		 * @param arrayA The A in Ax=b
		 * @param arrayB The b in Ax=b
		 * @return The x in Ax=b
		 */
		public static Matrix solveNNLS(double[][] arrayA, double[][] arrayB)
		{
			Matrix A = new Matrix(arrayA);
			Matrix b = new Matrix(arrayB);
			List<Integer> p = new ArrayList<Integer>();
			List<Integer> z = new ArrayList<Integer>();
			int i = 0;
			int xm = A.getColumnDimension();
			int xn = 1;
			while (i < A.getColumnDimension())
				z.add(i++);
			Matrix x = new Matrix(xm,xn);
			/*
			 * You need a finite number of iterations. Without this condition, the finite precision nature
			 * of the math being done almost makes certain that the <1e-15 conditions won't ever hold up.
			 * However, after so many iterations, it should at least be close to the correct answer.
			 * For the intrepid coder, however, one could replace this again with an infinite while
			 * loop and make the <1e-15 conditions into something like c*norm(A) or c*norm(b).
			 */
			for(int iterations = 0; iterations < 300*A.getColumnDimension()*A.getRowDimension(); iterations++)
			{
				//System.out.println(z.size() + " " + p.size());
				Matrix w = A.transpose().times(b.minus(A.times(x)));
				//w.print(7, 5);
				if(z.size() == 0 || isAllNegative(w))
				{
					//System.out.println("Computation should break");
					//We are done with the computation. Break here!
					break;//Should break out of the outer while loop.
				}
				//Step 4
				int t = z.get(0);
				double max = w.get(t, 0);
				for (i = 1; i < z.size(); i++)
				{
					if (w.get(z.get(i), 0) > max)
					{
						t = z.get(i);
						max = w.get(z.get(i), 0);
					}
				}
				//Step 5
				p.add(t);
				z.remove((Integer)t);
				boolean allPositive = false;
				while(!allPositive)
				{
					//Step 6
					Matrix Ep = new Matrix(b.getRowDimension(),p.size());
					for (i = 0; i < p.size(); i++)
						for (int j = 0; j < Ep.getRowDimension(); j++)
							Ep.set(j, i, A.get(j, p.get(i)));
					Matrix Zprime = Ep.solve(b);
					Ep = null;
					Matrix Z = new Matrix(xm,xn);
					for (i = 0; i < p.size(); i++)
						Z.set(p.get(i), 0, Zprime.get(i, 0));
					//Step 7
					allPositive = true;
					for (i = 0; i < p.size(); i++)
						allPositive &= Z.get(p.get(i), 0) > 0;
					if (allPositive)
						x = Z;
					else
					{
						double alpha = Double.MAX_VALUE;
						for (i = 0; i < p.size(); i++)
						{
							int q = p.get(i);
							if (Z.get(q,0) <= 0)
							{
								double xq = x.get(q, 0);
								if (xq / (xq - Z.get(q,0)) < alpha)
									alpha = xq / (xq - Z.get(q,0));
							}
						}
						//Finished getting alpha. Onto step 10
						x = x.plus(Z.minus(x).times(alpha));
						for (i = p.size() - 1; i >= 0; i--)
							if (Math.abs(x.get(p.get(i),0)) < 1e-15)//Close enough to zero, no?
								z.add(p.remove(i));
					}
				}
			}
			return x;
		}
		private static boolean isAllNegative(Matrix w)
		{
			boolean result = true;
			int m = w.getRowDimension();
			for (int i = 0; i < m; i++)
				result &= w.get(i, 0) <= 1e-15;
			return result;
		}
	}


