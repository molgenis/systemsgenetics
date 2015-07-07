/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package nl.systemsgenetics.cell_type_specific_ase;



import org.apache.commons.math3.analysis.MultivariateFunction;
import org.apache.commons.math3.optim.InitialGuess;
import org.apache.commons.math3.optim.MaxEval;
import org.apache.commons.math3.optim.PointValuePair;
import org.apache.commons.math3.optim.nonlinear.scalar.GoalType;
import org.apache.commons.math3.optim.nonlinear.scalar.ObjectiveFunction;
import org.apache.commons.math3.optim.nonlinear.scalar.noderiv.MultiDirectionalSimplex;
import org.apache.commons.math3.optim.nonlinear.scalar.noderiv.SimplexOptimizer;
import org.apache.commons.math3.optim.univariate.SearchInterval;
import org.junit.Assert;
import org.junit.Test;
/**
 *
 * @author adriaan
 */
public class OptimizerTest 
{
        
	
    public OptimizerTest() {
        MultivariateFunction branin = new Branin();
        SimplexOptimizer simplexOptimizer = new SimplexOptimizer(1e-6, 1e-6);

        PointValuePair solution = simplexOptimizer.optimize(new ObjectiveFunction(branin), new MaxEval(1000), GoalType.MINIMIZE, 
                        new SearchInterval(0.0, 1.0), new InitialGuess(new double[]{0.5, 0.5}), new MultiDirectionalSimplex(2));

        double [] x = solution.getFirst();
        System.out.println("Minimum:" + solution.getValue() + " found at, (" + (x[0]*15 - 5) + ", " + 15*x[1] + ")");
        System.out.println("Difference to true value: " + (solution.getValue() - 0.397887));

//		
//		if (!NumericalUtils.(solution.getValue(), 0.397887, 1e-6))
//		{
//			Assert.fail("Optimization did not work!");
//		}

		
    }

}
