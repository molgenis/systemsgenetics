package nl.systemsgenetics.genenetworkpathwayenrichment;

/**
 * Created by hwestra on 10/26/15.
 * Derived from Scott. A. Czepiel's implementation (http://czep.net/)
 */

import cern.colt.matrix.tdouble.DoubleMatrix2D;
import cern.colt.matrix.tdouble.algo.DenseDoubleAlgebra;
import cern.colt.matrix.tdouble.impl.DenseDoubleMatrix2D;
import cern.jet.stat.Gamma;
import umcg.genetica.io.text.TextFile;

import java.io.IOException;
import java.text.DecimalFormat;

public class LogisticRegressionOptimized {

    int max_iter = 100;
    double EPSILON = 1E-6;
    DoubleMatrix2D xtwx;
    private double[] g;
    private double[] numer;
    private DoubleMatrix2D pi;
    private DoubleMatrix2D H;


    public boolean debug = false;
    public boolean printIters = false;
    public String output;
    public boolean flipCoding;

    public LogisticRegressionOptimized() {

    }

    public LogisticRegressionOptimized(int maxiter) {
        max_iter = maxiter;
    }

    public LogisticRegressionOptimized(int maxiter, double epsilon) {
        max_iter = maxiter;
        EPSILON = epsilon;
    }

    public void setMax_iter(int max_iter) {
        this.max_iter = max_iter;
    }

    public void setEPSILON(double EPSILON) {
        this.EPSILON = EPSILON;
    }


    public LogisticRegressionResult binomial(double[] y, double[][] x) {
        if (x.length != y.length) {
            throw new IllegalArgumentException("Error in logistic regression: length of y does not match length of x");
        }
        DoubleMatrix2D xtmp = new DenseDoubleMatrix2D(x);
        return binomial(y, xtmp);
    }

    public LogisticRegressionResult binomial(double[] y, DoubleMatrix2D x) {
        if (x.rows() != y.length) {
            throw new IllegalArgumentException("Error in logistic regression: length of y does not match length of x");
        }


        DoubleMatrix2D ytmp = new DenseDoubleMatrix2D(y.length, 2);
        for (int i = 0; i < y.length; i++) {
            int index = (int) Math.abs(y[i]);
            ytmp.setQuick(i, index, 1);
        }

        return mlelr(x, ytmp);
    }

    public LogisticRegressionResult multinomial(double[][] y, double[][] x) {
        if (x.length != y.length) {
            throw new IllegalArgumentException("Error in logistic regression: length of y does not match length of x");
        }
        DoubleMatrix2D xtmp = new DenseDoubleMatrix2D(x);
        DoubleMatrix2D ytmp = new DenseDoubleMatrix2D(y);

        return multinomial(ytmp, xtmp);
    }

    public LogisticRegressionResult multinomial(double[][] y, DoubleMatrix2D xprime) {
        if (xprime.rows() != y.length) {
            throw new IllegalArgumentException("Error in logistic regression: length of y does not match length of x");
        }
        DoubleMatrix2D ytmp = new DenseDoubleMatrix2D(y);
        return multinomial(ytmp, xprime);
    }

    private DoubleMatrix2D formatBinomialData(DoubleMatrix2D y) {
        DoubleMatrix2D ytmp = new DenseDoubleMatrix2D(y.rows(), y.columns() + 1);
        for (int i = 0; i < y.rows(); i++) {
//			int
//			if (flipCoding) {
            int index = (int) Math.abs(y.get(i, 0) - 1);
            if (flipCoding) {
                index = (int) y.get(i, 0);
            }
//			}
            ytmp.set(i, index, 1);
        }
        return ytmp;
    }

    public LogisticRegressionResult multinomial(DoubleMatrix2D y, DoubleMatrix2D x) {
        if (x.rows() != y.rows()) {
            throw new IllegalArgumentException("Error in logistic regression: length of y does not match length of x");
        }

        DoubleMatrix2D ytmp = null;

//		System.out.println("ytmp: " + ytmp.rows() + "\t" + ytmp.columns());


        if (y.columns() == 1) {
            ytmp = formatBinomialData(y);
        } else {
            ytmp = new DenseDoubleMatrix2D(y.rows(), y.columns() + 1);
            int lastCol = ytmp.columns() - 1;
            for (int i = 0; i < ytmp.rows(); i++) {
                int nrnull = 0;
                for (int j = 0; j < y.columns(); j++) {
                    // y[ind][disease1] =0
                    // y[ind][disease2] =0 --> y[i][0] = 1

                    // y[ind][disease1] =1
                    // y[ind][disease2] =0 --> y[i][1] = 1

                    // y[ind][disease1] =0
                    // y[ind][disease2] =1 --> y[i][2] = 1
                    // y=1 --> y[i][1] = 1
                    double val = y.getQuick(i, j);
                    if (val == 0) {
                        nrnull++;
                    }

                    int index = y.columns() - j - 1;
                    // col0 -> col 1 //
                    // col1 -> col 0
                    ytmp.setQuick(i, index, val);
                }
                // set first column to 1 if 0 everywhere else (individual is in reference group)
                if (nrnull == y.columns()) {
                    ytmp.setQuick(i, lastCol, 1);
                }
            }
        }


        LogisticRegressionResult r = mlelr(x, ytmp);
        if (y.columns() == 1) {
            // TODO: find out what is wrong here...
            // flip beta's
            if (flipCoding) {
                double[][] betas = r.getBeta();
                for (int i = 0; i < betas[0].length; i++) {
                    betas[0][i] *= -1;
                }
            }
        }
        if (output != null || debug) {
            System.out.println("Input for the model will be written here: " + output);

            System.out.println(y.columns());
            System.out.println(ytmp.columns());

//			r = mlelr(x, ytmp);
//			for (int i = 0; i < r.getBeta().length; i++) {
//				for (int j = 0; j < r.getBeta()[i].length; j++) {
//					System.out.println(i + "\t" + j + "\tbeta: " + r.getBeta()[i][j] + "\tse: " + r.getStderrs()[i][j]);
//				}
//			}
            System.out.println("Deviance: " + r.getDeviance());


            int printrows = 10;
            if (y.rows() < 100) {
                printrows = y.rows();
            }
            try {

                TextFile outputf = null;
                if (output != null) {
                    outputf = new TextFile(output, TextFile.W);
                }
                for (int i = 0; i < y.rows(); i++) {
                    String ln = "" + i;
                    String lnout = "";
                    for (int j = 0; j < ytmp.columns(); j++) {
                        ln += "\t" + ytmp.getQuick(i, j);
                    }
                    ln += "\torig";
                    for (int j = 0; j < y.columns(); j++) {
                        ln += "\t" + y.getQuick(i, j);
                        if (j == 0) {
                            lnout += y.getQuick(i, j);
                        } else {
                            lnout += "\t" + y.getQuick(i, j);
                        }
                    }
                    ln += "\torigX:\t";
                    for (int j = 1; j < x.columns(); j++) {
                        ln += "\t" + x.getQuick(i, j);
//						if (j > 0) {
                        lnout += "\t" + x.getQuick(i, j);
//						}
                    }
                    if (outputf != null) {
                        outputf.writeln(lnout);
                    }
                    if (i < printrows) {
                        System.out.println(ln);
                    }
                }
                if (outputf != null) {
                    outputf.close();
                }

                return null;
            } catch (IOException e) {
                e.printStackTrace();
            }
        }


        return r;
    }


    public LogisticRegressionResult mlelr(
            DoubleMatrix2D X,
            DoubleMatrix2D Y) {

        double[] n = new double[Y.rows()];
        for (int i = 0; i < Y.rows(); i++) {
            n[i] = 1;
        }

        int i, j, k;

        int iter = 0;
        boolean converged = false;

        int N = X.rows();
        int K = X.columns();
        int J = Y.columns();

        int kJMinusOne = K * (J - 1);
        // start values for beta can be determined by linear regression of log(pij/piJ) on design matrix x
//		System.out.println("Nr betas: " + kJMinusOne);
        double[] beta = new double[kJMinusOne];

        double[] beta0 = new double[kJMinusOne];
//		double[] diff = new double[kJMinusOne];
//		boolean[] diffb = new boolean[kJMinusOne];
        //double[] beta_inf = new double[kJMinusOne];
        if (xtwx == null) {
            xtwx = new DenseDoubleMatrix2D(kJMinusOne, kJMinusOne);
        } else {
            if (xtwx.rows() != kJMinusOne) {
                xtwx = new DenseDoubleMatrix2D(kJMinusOne, kJMinusOne);
            } else {
                xtwx.assign(0);
            }
        }


        double[] loglike = new double[1];
        double[] deviance = new double[1];
        double loglike0 = 0;


        double deviance0 = 0;

        while (iter < max_iter && !converged) {

            for (k = 0; k < kJMinusOne; k++) {
                beta0[k] = beta[k];
            }

            // ! nr function needs error handling
            nr(X, Y, n, J, N, K, beta0, beta, xtwx, loglike, deviance);


            if (loglike[0] < loglike0 && iter > 0) {
                // backtracking code
                // run subiterations to halve the distance to prior iteration
                // until difference in log_like increases (which is when the model has converged)
                // remember: ml is about maximizing loglike
            }

//			// test for infinity of beta
//			for (k = 0; k < kJMinusOne; k++) {
//				if (beta_inf[k] != 0) {
//					beta[k] = beta_inf[k];
//				} else {
//					//  Math.sqrt(xtwx[k][k]) contains the variance of beta[k]
//					double betak = beta[k];
//					double absbeta = Math.abs(betak);
//					if (absbeta > (5d / xrange[k]) && Math.sqrt(xtwx[k][k]) >= (3 * absbeta)) {
//						beta_inf[k] = betak;
//					}
//				}
//			}

            // test for convergence
            converged = true;

//			double deltadeviance = Math.abs(deviance0 - deviance[0]);
//			if (deltadeviance / (0.1 + Math.abs(deviance[0])) > EPSILON) {
//				converged = false;
//			}


            for (k = 0; k < kJMinusOne; k++) {
                double beta0k = beta0[k];
                double delta = Math.abs(beta[k] - beta0k);
                if (Double.isNaN(delta) || Double.isInfinite(delta) || delta > EPSILON * Math.abs(beta0k)) {
                    converged = false;
//					diffb[k] = false;
                    break;
                }
            }

            if (Math.abs(deviance0 - deviance[0]) > EPSILON) {
                converged = false;
            }


            if (printIters) {


                System.out.println("Deviance: " + deviance[0] + " previous deviance: " + deviance0);
                for (k = 0; k < kJMinusOne; k++) {
                    double beta0k = beta0[k];
                    double delta = Math.abs(beta[k] - beta0k);
                    double eps = EPSILON * Math.abs(beta0k);
                    DecimalFormat f = new DecimalFormat("#.######");
                    String prefix = "";
                    if ((delta > eps)) {
                        prefix = (char) 27 + "[36m";
                    } else {
                        prefix = (char) 27 + "[37m";
                    }


                    System.out.println(prefix + "Col: " + k + "\tmodel converged:\t" + converged + "\titer:" + iter + "\tb0: " + f.format(beta0k) + "\tbk: " + f.format(beta[k]) + "\td: " + f.format(delta) + "\teps: " + f.format(eps) + "\tparam converged: " + !(delta > eps));
//					}


                }
            }

            if (iter == 0) {
                loglike0 = loglike[0];
            }
            deviance0 = deviance[0];

            iter++;
        }


        if (debug) {
            System.out.println();
            System.out.println("Output at final iteration...");
            System.out.println("Model converged?\t" + converged);
            for (k = 0; k < kJMinusOne; k++) {
                double beta0k = beta0[k];
                double delta = Math.abs(beta[k] - beta0k);
                double eps = EPSILON * Math.abs(beta0k);

                System.out.println("Col: " + k + "\tIter " + iter + "\tb0: " + beta0k + "\tbk: " + beta[k] + "\td: " + delta + "\teps: " + eps + "\tparam converged: " + !(delta > eps));
            }

//			System.out.println("Will try to rerun the MLE");
//			if (!converged) {
//
//				System.out.println("Rerunning MLE");
//
//				// try another round of MLE, but now reset the conflicting parameter's estimates...
//				for (k = 0; k < kJMinusOne; k++) {
//					double beta0k = beta0[k];
//					double delta = Math.abs(beta[k] - beta0k);
//					double eps = EPSILON * Math.abs(beta0k);
//
////					if (delta > eps) {
//					beta[k] = beta[k];
//					beta0[k] = beta0[k];
////					}
//				}
//
//				iter = 0;
//
//				while (iter < max_iter && !converged) {
//
//					for (k = 0; k < kJMinusOne; k++) {
//						beta0[k] = beta[k];
//					}
//
//					// ! nr function needs error handling
//					nr(X, Y, n, J, N, K, beta0, beta, xtwx, loglike, deviance);
//
//					if (loglike[0] < loglike0 && iter > 0) {
//						// backtracking code
//						// run subiterations to halve the distance to prior iteration
//						// until difference in log_like increases (which is when the model has converged)
//						// remember: ml is about maximizing loglike
//					}
//
////			// test for infinity of beta
////			for (k = 0; k < kJMinusOne; k++) {
////				if (beta_inf[k] != 0) {
////					beta[k] = beta_inf[k];
////				} else {
////					//  Math.sqrt(xtwx[k][k]) contains the variance of beta[k]
////					double betak = beta[k];
////					double absbeta = Math.abs(betak);
////					if (absbeta > (5d / xrange[k]) && Math.sqrt(xtwx[k][k]) >= (3 * absbeta)) {
////						beta_inf[k] = betak;
////					}
////				}
////			}
//
//					// test for convergence
//					converged = true;
//					for (k = 0; k < kJMinusOne; k++) {
//						double beta0k = beta0[k];
//						double delta = Math.abs(beta[k] - beta0k);
//						if (Double.isNaN(delta) || Double.isInfinite(delta) || delta > EPSILON * Math.abs(beta0k)) {
//							converged = false;
////					diffb[k] = false;
//							break;
//						}
//					}
//
//					if (iter == 0) {
//						loglike0 = loglike[0];
//					}
//
//					if (printIters) {
//						for (k = 0; k < kJMinusOne; k++) {
//							double beta0k = beta0[k];
//							double delta = Math.abs(beta[k] - beta0k);
//							double eps = EPSILON * Math.abs(beta0k);
//							DecimalFormat f = new DecimalFormat("#.######");
//							String prefix = "";
//							if ((delta > eps)) {
//								prefix = (char) 27 + "[35m";
//							} else {
//								prefix = (char) 27 + "[37m";
//							}
//							System.out.println(prefix + "Col: " + k + "\tmodel converged:\t" + converged + "\titer:" + iter + "\tb0: " + f.format(beta0k) + "\tbk: " + f.format(beta[k]) + "\td: " + f.format(delta) + "\teps: " + f.format(eps) + "\tparam converged: " + !(delta > eps));
////					}
//
//
//						}
//					}
//					iter++;
//				}
//			}
        }

        // chi2 tests for significance
//			double chi1 = 2 * (loglike[0] - loglike0);
//			int df1 = (K * (J - 1)) - J - 1;
//			double chiTest = 1 - ChiSquare.getP(df1, chi1);
//
//			double chi2 = deviance[0];
//			int df2 = (N * (J - 1)) - (K * (J - 1));
//			double chiTest2 = 1 - ChiSquare.getP(df2, chi2);

//			double[] sigprms = new double[kJMinusOne];
        double[] stderrs = new double[kJMinusOne];
//		double[] wald = new double[kJMinusOne];
        for (i = 0; i < kJMinusOne; i++) {
            double xtwxii = xtwx.get(i, i);
            if (xtwxii > 0) {
                stderrs[i] = Math.sqrt(xtwxii);
//				wald[i] = Math.pow(beta[i] / stderrs[i], 2);
//					sigprms[i] = 1 - ChiSquare.getP(1, wald[i]);
            } else {
//					sigprms[i] = -1;
            }
        }

        int nrDiseases = Y.columns() - 1;
        int nrPredictors = X.columns();
//			System.out.println(outputYdim + " x " + outputXdim);
        double[][] outputbeta = new double[nrDiseases][nrPredictors];
        double[][] outputse = new double[nrDiseases][nrPredictors];

        // output is inverted for some reason..
        int ctr = 0;

//		if (debug) {
//			System.out.println("Converged: " + converged);
//			System.out.println(beta.length);
//			System.out.println(X.columns());
//
//			System.out.println("model parameters:");
//			for (int q = 0; q < beta.length; q++) {
//				System.out.println(q + "\t" + beta[q] + "\t" + stderrs[q]);
//			}
//			System.out.println();
//
//		}
//			System.out.println(betasPerDisase);
        for (int disease = 0; disease < nrDiseases; disease++) {
            // order:
            // 2 disases -> disease0 --> 2 - 0 - 1 = 1
            //           -> disease1 --> 2 - 1 - 1 = 0
            for (int col = 0; col < nrPredictors; col++) {
//					System.out.println(disease + "\t" + col);
                outputbeta[disease][col] = beta[ctr];
                outputse[disease][col] = stderrs[ctr];
                ctr++;
            }
        }
        if (debug) {
            return new LogisticRegressionResult(outputbeta, outputse, deviance[0]);
        }

        if (converged) {
            return new LogisticRegressionResult(outputbeta, outputse, deviance[0]);
        } else {
            return null;
        }
    }


    private int nr(
            DoubleMatrix2D X,
            DoubleMatrix2D Y,
            double[] n,
            int J,
            int N,
            int K,
            double[] beta0,
            double[] beta1,
            DoubleMatrix2D xtwx,
            double[] loglike,
            double[] deviance
    ) {
        int kJMinusOne = (K * (J - 1));

        // if this is the first iteration, initialize (and save some GC time in the following iterations)
        if (pi == null) {
            pi = new DenseDoubleMatrix2D(N, J);
            H = new DenseDoubleMatrix2D(kJMinusOne, kJMinusOne);
            g = new double[kJMinusOne];
            numer = new double[J];
        } else {
            // check the size
            if (pi.rows() != N || pi.columns() != J) {
                pi = new DenseDoubleMatrix2D(N, J);
            } else {
                pi.assign(0);
            }
            if (H.rows() != kJMinusOne || H.columns() != kJMinusOne) {
                H = new DenseDoubleMatrix2D(kJMinusOne, kJMinusOne);
            } else {
                H.assign(0);
            }
            if (g.length != kJMinusOne) {
                g = new double[kJMinusOne];
            } else {
                for (int i = 0; i < g.length; i++) {
                    g[i] = 0;
                }
            }
            if (numer.length != J) {
                numer = new double[J];
            } else {
                for (int i = 0; i < numer.length; i++) {
                    numer[i] = 0;
                }
            }
        }


        double denom, q1, w1, w2, sum1;
        double devtmp;

        int i, j, k, jj, kk, kprime, jprime;

        loglike[0] = 0;
        deviance[0] = 0;

        int jMinusOne = J - 1;
        for (i = 0; i < N; i++) {
            denom = 1d;
            jj = 0;
            for (j = 0; j < jMinusOne; j++) {
                sum1 = 0;
                for (k = 0; k < K; k++) {
                    sum1 += X.get(i, k) * beta0[jj++];
                }
                numer[j] = Math.exp(sum1);
                denom += numer[j];

            }

            /* calculate predicted probabilities */
            for (j = 0; j < J - 1; j++) {
                pi.setQuick(i, j, numer[j] / denom);
            }

            /* omitted category */
            pi.setQuick(i, j, 1.0d / denom);

            /* increment log likelihood */
            loglike[0] += Gamma.logGamma(n[i] + 1);
            for (j = 0; j < J; j++) {
                double pij = pi.getQuick(i, j);
                double yij = Y.getQuick(i, j);
                loglike[0] = loglike[0] - Gamma.logGamma(yij + 1) + (yij * Math.log(pij));
//			}
//
//			/* increment deviance */
//			for (j = 0; j < J; j++) {
//				double yij = Y.get(i, j);
                if (yij > 0) {
                    devtmp = 2 * yij * Math.log(yij / (n[i] * pij));
                } else {
                    devtmp = 0;
                }
                deviance[0] += devtmp;
            }

            double ni = n[i];
            /* increment first and second derivatives */
            for (j = 0, jj = 0; j < J - 1; j++) {

                double yij = Y.getQuick(i, j);
                double pij = pi.getQuick(i, j);
                /* terms for first derivative, see Eq. 32 */
                q1 = yij - (ni * pij);

                /* terms for second derivative, see Eq. 37 */
                w1 = ni * pij * (1 - pij);

                for (k = 0; k < K; k++) {

                    double xik = X.getQuick(i, k);
                    /* first derivative term in Eq. 23 */
                    g[jj] += q1 * xik;

                    /* increment the current pop's contribution to the 2nd derivative */

                    /* jprime = j (see Eq. 37) */
                    kk = jj - 1;
                    for (kprime = k; kprime < K; kprime++) {
                        kk += 1;
                        double hjjkk = H.getQuick(jj, kk);
                        hjjkk += w1 * xik * X.getQuick(i, kprime);
                        H.setQuick(jj, kk, hjjkk);
                        H.setQuick(kk, jj, hjjkk);
                    }

                    /* jprime != j (see Eq. 37) */

                    for (jprime = j + 1; jprime < J - 1; jprime++) {
                        w2 = -ni * pij * pi.getQuick(i, jprime);
                        for (kprime = 0; kprime < K; kprime++) {
                            kk += 1;
                            double hjjkk = H.getQuick(jj, kk);
                            hjjkk += w2 * xik * X.getQuick(i, kprime);
                            H.setQuick(jj, kk, hjjkk);
                            H.setQuick(kk, jj, hjjkk);
                        }
                    }
                    jj++;
                }
            }

        }

        /* compute xtwx * beta0 + x(y-mu) (see Eq. 40) */
        for (i = 0; i < kJMinusOne; i++) {
            sum1 = 0;
            for (j = 0; j < kJMinusOne; j++) {
                sum1 += H.getQuick(i, j) * beta0[j];
            }
            g[i] += sum1;
        }

        /* invert xtwx */
        if (cholesky(H)) return -1;
        if (backsub(H)) return -2;
        if (trimult(H, xtwx)) return -3;

        /* solve for new betas */
        for (i = 0; i < kJMinusOne; i++) {
            sum1 = 0;
            for (j = 0; j < kJMinusOne; j++) {
                sum1 += xtwx.getQuick(i, j) * g[j];
            }
            beta1[i] = sum1;
        }

        return 0;
    }

    private DenseDoubleAlgebra dda = new DenseDoubleAlgebra();

    private boolean trimult(DoubleMatrix2D in, DoubleMatrix2D out) {
        int i, j, k, m;
        double sum;
        int order = in.rows();

        for (i = 0; i < order; i++) {
            for (j = 0; j < order; j++) {
                sum = 0;
                if (i > j) {
                    m = i;
                } else {
                    m = j;
                }
                for (k = m; k < order; k++) {
                    sum += in.getQuick(i, k) * in.getQuick(j, k);
                }
                out.setQuick(i, j, sum);
            }
        }
        return false;
    }

    private boolean backsub(DoubleMatrix2D x) {
        int i, j, k;
        double sum;
        int order = x.rows();

        if (x.get(0, 0) == 0) {
            return true; // throw new ArithmeticException("Problem with backsubstitution: x[0][0] == 0d");
        }

        x.setQuick(0, 0, (1d / x.getQuick(0, 0)));
        for (i = 1; i < order; i++) {
            if (x.getQuick(i, i) == 0) {
                return true; // throw new ArithmeticException("Problem with backsubstitution: x[i][i] == 0d");
            }
            x.setQuick(i, i, (1d / x.getQuick(i, i)));
            for (j = 0; j < i; j++) {
                sum = 0;
                for (k = j; k < i; k++) {
                    sum += x.getQuick(j, k) * x.getQuick(k, i);
                }
                x.setQuick(j, i, (-sum * x.getQuick(i, i)));
            }
        }
        return false;
    }

    private boolean cholesky(DoubleMatrix2D x) {
        int i, j, k;
        double sum;
        int order = x.rows();

        for (i = 0; i < order; i++) {
            sum = 0;
            double xii = x.getQuick(i, i);
            for (j = 0; j < i; j++) {
                double xji = x.getQuick(j, i);
                sum += xji * xji;
            }
            if (sum >= xii) {
                return true; // throw new ArithmeticException("Problem with Cholesky operation: sum>x[i][i]");
            }
            x.setQuick(i, i, Math.sqrt(xii - sum));
            for (j = i + 1; j < order; j++) {
                sum = 0;
                for (k = 0; k < i; k++) {
                    sum += x.getQuick(k, i) * x.getQuick(k, j);
                }
                x.setQuick(i, j, (x.getQuick(i, j) - sum) / x.getQuick(i, i));
            }
        }
        return false;
    }


}
