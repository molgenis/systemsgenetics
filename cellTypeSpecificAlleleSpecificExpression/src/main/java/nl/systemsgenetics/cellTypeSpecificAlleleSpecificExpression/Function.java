/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package nl.systemsgenetics.cellTypeSpecificAlleleSpecificExpression;

import org.apache.commons.math3.analysis.MultivariateFunction;
import org.ejml.simple.SimpleMatrix;
/**
 *
 * @author adriaan
 */
public interface Function extends MultivariateFunction
{
	
	public int getDim();

	/**
	 * Vector version of the function evaluation
	 * Evaluate the function at x
	 * @param x
	 * @return
	 */
        
	public double value(double [] x);
	
	/**
	 * Matrix version of the function evaluation
	 * @param xx - n x d array where each row represents an input to the function
	 * @return
	 */
	public double [] value(SimpleMatrix xx);
	
}

