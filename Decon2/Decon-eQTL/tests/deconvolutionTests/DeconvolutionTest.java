package deconvolutionTests;

import static org.junit.Assert.*;

import java.io.File;
import java.util.ArrayList;
import java.util.Arrays;

import org.apache.commons.io.FileUtils;
import org.apache.commons.io.LineIterator;
import org.junit.Before;
import org.junit.Test;
import Decon_eQTL.Deconvolution;

public class DeconvolutionTest {
	private double [] y;
	private double [][] x;
	private double [][] xb;
	private double sumOfSquaresA;
	private double sumOfSquaresNoInterceptA;
	private double sumOfSquaresNoInterceptB;
	private double sumOfSquaresB;
	private int degreesOfFreedomA;
	private int degreesOfFreedomB;

    @Before
    public void init() {
    	// set the variables that will be used for multiple tests
		y = new double[] {-0.48812477, 0.33458213, -0.52754476, -0.79863471, -0.68544309, -0.12970239, 0.02355622, -0.31890850, 0.34725819,  0.08108851};
		x = new double[][] {{1,0}, {0,0}, {1,0}, {2,1}, {0,1}, {0,0}, {1,0}, {0,0}, {1,0}, {0,0}};
		xb = new double[][] {{1,0,0}, {0,0,0}, {1,0,0}, {2,1,2}, {0,1,0}, {0,0,0}, {1,0,0}, {0,0,0}, {1,0,0}, {0,0,0}};
		//InteractionModel noInteraction = new InteractionModel();
		//noInteraction.SetExpressionValues(y);
		//noInteraction.SetObservedValues(x);
		//InteractionModel interaction = new InteractionModel();
		//interaction.SetExpressionValues(y);
		//interaction.SetObservedValues(xb);
		try {
		//	sumOfSquaresA = Deconvolution.calculateSumOfSquaresNNLS(noInteraction, false);
		//sumOfSquaresB = Deconvolution.calculateSumOfSquaresOLS(interaction, false);
		//sumOfSquaresNoInterceptA = Deconvolution.calculateSumOfSquaresOLS(noInteraction, true);
		//sumOfSquaresNoInterceptB = Deconvolution.calculateSumOfSquaresOLS(interaction, true);
		} catch (Exception e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		} 
		degreesOfFreedomA = y.length - (x[0].length + 1);
		degreesOfFreedomB = y.length - (xb[0].length + 1);
    }
    
	@Test
	public void calculateSumOfSquaresTest() {
		// test if correct sum of squares is given for simple model (x), model with interaction (xb) and both with and without intercept
		assertEquals(7, degreesOfFreedomA, 0);
		assertEquals(6, degreesOfFreedomB, 0);
		assertEquals(0.7854, sumOfSquaresNoInterceptA, 0.0001);
		assertEquals(0.7798, sumOfSquaresA, 0.0001);
		assertEquals(0.7707, sumOfSquaresNoInterceptB, 0.0001);
		assertEquals(0.7705, sumOfSquaresB, 0.0001);
	}
	
	@Test
	public void anovaTest() {
		/**pvalue same as in R with the following code:
		 *   test_trait <- c( -0.48812477 , 0.33458213, -0.52754476, -0.79863471, -0.68544309, -0.12970239,  0.02355622, -0.31890850,0.34725819 , 0.08108851)
		 *   geno_A <- c(1, 0, 1, 2, 0, 0, 1, 0, 1, 0)
		 *   geno_B <- c(0, 0, 0, 1, 1, 0, 0, 0, 0, 0)
		 *   fit <- lm(test_trait ~ geno_A+geno_B -1 )
		 *   fit2 <- lm(test_trait ~ geno_A + geno_B + geno_A:geno_B-1)
		 *   anova(fit,fit2)
		 **/
		// test if correct sum of squares is given for simple model (x), model with interaction (xb) and both with and without intercept
		//assertEquals("Anova of lm fits with intercept", 0.782, Deconvolution.anova(sumOfSquaresA, sumOfSquaresB,degreesOfFreedomA, degreesOfFreedomB, false), 0.001);
		// test if correct sum of squares is given for simple model (x), model with interaction (xb) and both with and without intercept
		//assertEquals("Anova of lm fits without intercept", 0.711, Deconvolution.anova(sumOfSquaresNoInterceptA, sumOfSquaresNoInterceptB,degreesOfFreedomA, degreesOfFreedomB, true), 0.001);
	}
	
	@Test
	public void mainTest() throws Exception {
		File counts = new File("tests/resources/counts.txt");
		File expTable = new File("tests/resources/expTable_Corrected_addMean_snpname.txt");
		File dsgTable = new File("tests/resources/dsgTable_testing_snpname.txt");
		String[] args = {"-o","tests/resources/decovolutionTestResult.txt","-c",counts.getAbsolutePath(),"-e",
						 expTable.getAbsolutePath(), "-g", dsgTable.getAbsolutePath()};
		Deconvolution.main(args);
		LineIterator deconResults = FileUtils.lineIterator(new File("tests/resources/decovolutionTestResult.txt"), "UTF-8");
		LineIterator deconExpected = FileUtils.lineIterator(new File("tests/resources/deconExpected.txt"), "UTF-8");
		//test if header is same
		assertEquals("File header the same",deconExpected.next(),deconResults.next());
		while (deconResults.hasNext() && deconExpected.hasNext()){
			ArrayList<String> deconResultsStringVector = new ArrayList<String>(Arrays.asList(deconResults.next().split("\t")));
			ArrayList<String> deconExpectedStringVector = new ArrayList<String>(Arrays.asList(deconExpected.next().split("\t")));
			assertEquals("QTL name the same", deconExpectedStringVector.remove(0), deconResultsStringVector.remove(0));
			//double[] deconResultsVector = Deconvolution.StringVectorToDoubleVector(deconResultsStringVector);
			//double[] deconExpectedVector = Deconvolution.StringVectorToDoubleVector(deconExpectedStringVector);
			//for (int i = 0; i < deconResultsVector.length; i++){
			//	assertEquals("R p-value vs Java p-value", deconExpectedVector[i], deconResultsVector[i], 0.0002);
			//}
		}
	}
	
	@Test
	public void permutationTestTest() throws Exception {
		LineIterator expressionIterator = FileUtils.lineIterator(new File("tests/resources/expTable_Corrected_addMean_head.txt"), "UTF-8");
		LineIterator genotypeIterator = FileUtils.lineIterator(new File("tests/resources/dsgTable_testing_head.txt"), "UTF-8");
		LineIterator cellcountIterator = FileUtils.lineIterator(new File("tests/resources/counts.txt"), "UTF-8");
		ArrayList<String> celltypes = new ArrayList<String>(Arrays.asList(cellcountIterator.next().split("\t")));
		celltypes.removeAll(Arrays.asList("", null));
		cellcountIterator.close();
		expressionIterator.next();
		genotypeIterator.next();
		//PermutationResult permutation_result = Deconvolution.permutationTest("tests/resources/expTable_Corrected_addMean_head.txt", "tests/resources/dsgTable_testing_head.txt", cellcountTable, 10, 1, true);
		/* TODO: write tests on permutation_result */
		fail("Not yet implemented");
	}
}
