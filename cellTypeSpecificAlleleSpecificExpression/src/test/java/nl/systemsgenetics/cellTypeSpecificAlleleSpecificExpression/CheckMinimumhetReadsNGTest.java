/*
 * Copyright (C) 2016 adriaan
 *
 * This program is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program.  If not, see <http://www.gnu.org/licenses/>.
 */
package nl.systemsgenetics.cellTypeSpecificAlleleSpecificExpression;

import java.util.ArrayList;
import junit.framework.Assert;
import static nl.systemsgenetics.cellTypeSpecificAlleleSpecificExpression.LinearModelTestNGTest.asRef;
import static org.testng.Assert.*;
import org.testng.annotations.AfterClass;
import org.testng.annotations.AfterMethod;
import org.testng.annotations.BeforeClass;
import org.testng.annotations.BeforeMethod;
import org.testng.annotations.Test;

/**
 *
 * @author adriaan
 */
public class CheckMinimumhetReadsNGTest {
    
    static int[] asRef      = { 6, 7, 0, 0, 0, 0, 6, 29, 9, 4, 5, 3, 0, 9, 0, 0, 9, 3, 7, 5, 4, 0, 3 , 10, 3, 0,11, 7, 7, 0, 8, 7,11, 8, 0, 0, 13,10,10, 3,10, 8,16, 8, 10, 0, 9, 0, 2, 2, 0, 0, 4, 2, 8 ,10, 0 };
    static int[] asAlt      = { 1, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 2, 0, 0, 1, 2, 0, 0, 0, 0, 1, 0, 0, 0, 0, 1, 2, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 };
    static double[] dispersion = {0.3164, 0.3527, 0.3917, 0.3236, 0.3525, 0.3574, 0.3621, 0.2974, 0.3274, 0.3309, 0.2965, 0.3058, 0.3466, 0.3197, 0.3353, 0.3514, 0.3171, 0.3329, 0.3806, 0.3566, 0.3235, 0.3275, 0.3245, 0.2899, 0.3007, 0.3119, 0.3405, 0.3188, 0.3496, 0.3584, 0.3325, 0.3242, 0.3496, 0.3207, 0.3014, 0.3168, 0.3108, 0.2883, 0.3319, 0.3087, 0.2937, 0.3301, 0.3241, 0.3241, 0.3069, 0.2958, 0.3195, 0.3279, 0.3362, 0.3464, 0.3817, 0.3281, 0.3441, 0.3849, 0.3232, 0.3206, 0.2932 };
    static double[] cellProp   = {0.5650, 0.5240, 0.5670, 0.5270, 0.5760, 0.4400, 0.5140, 0.3610, 0.5470, 0.4330, 0.4480, 0.5950, 0.5130, 0.5090, 0.4900, 0.5220, 0.4820, 0.7200, 0.6050, 0.4920, 0.4820, 0.4580, 0.6930, 0.6250, 0.5490, 0.5150, 0.3950, 0.6000, 0.5560, 0.5480, 0.4620, 0.6060, 0.5390, 0.4240, 0.3900, 0.6320, 0.6060, 0.3780, 0.5930, 0.4830, 0.4490, 0.5090, 0.3850, 0.4530, 0.4730, 0.6390, 0.3680, 0.5620, 0.6730, 0.5280, 0.5060, 0.6310, 0.5620, 0.7970, 0.6430, 0.6360, 0.4540};
    //turn this into individual snp data.

    private ArrayList<IndividualSnpData> snpData = new ArrayList<IndividualSnpData>();
    private double asRatio;

    
    
    public CheckMinimumhetReadsNGTest() {
    }
    
    
    // TODO add test methods here.
    // The methods must be annotated with annotation @Test. For example:
    //
    // @Test
    // public void hello() {}

    
    @BeforeClass
    public void setUpClass() throws Exception {
        int refSum = 0;
        int altSum = 0;
        for(int i =0;i<57;i++){
            
            String snpLine = "0\t1\trsTest\tA\tG\t";
            snpLine += Integer.toString(asRef[i]) + "\t";
            snpLine += Integer.toString(asAlt[i]) + "\t";
            snpLine += "0\t";
            snpLine += "[A, G]";

            String thisName = "ind" + Integer.toString(i);
            IndividualSnpData thisSnp = new IndividualSnpData(thisName, snpLine);
            thisSnp.setCellTypeProp(cellProp[i]);
            thisSnp.setDispersion(dispersion[i]);
            
            snpData.add(thisSnp);
        }

    }

    
    @AfterClass
    public static void tearDownClass() throws Exception {
    }

    @Test
    public void BinomialTestSanity(){
        
        GlobalVariables.minHetReads = 0.1;
        double minhets = GlobalVariables.minHetReads;
        BinomialTest checkBinomial = new BinomialTest(snpData);
        
        Assert.assertEquals(checkBinomial.getChromosome(),"0");
        Assert.assertEquals(checkBinomial.binomRatio,0.8571428571428571);

    }
    
    
    @Test
    public void BetaBinomialTestSanity() throws Exception{
        
        GlobalVariables.minHetReads = 0.1;
        double minhets = GlobalVariables.minHetReads;
        BetaBinomialTest checkBetaBinomial = new BetaBinomialTest(snpData);
        
        Assert.assertEquals(checkBetaBinomial.getChromosome(),"0");
        //please note, this is NOT the same as the binom ratio of the binomial test, as this is 
        // calculated as the alphaParam / (alphaParam + betaParam)
        // in the BetaBinomialTest.java
        Assert.assertEquals(checkBetaBinomial.binomRatio,0.8181818181818182);

    }
    
    
    
    @BeforeMethod
    public void setUpMethod() throws Exception {
    }

    @AfterMethod
    public void tearDownMethod() throws Exception {
    }
}
