/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package nl.systemsgenetics.cellTypeSpecificAlleleSpecificExpression;

import static java.lang.Math.abs;
import junit.framework.Assert;
import org.testng.annotations.Test;

/**
 *
 * @author adriaan
 */
public class LikelihoodFunctionsNGTest {
    
    public LikelihoodFunctionsNGTest() {
    }
    // TODO add test methods here.
    // The methods must be annotated with annotation @Test. For example:
    //
    // @Test
    // public void hello() {}
    
    @Test
    public void testBinomialLikelihood1(){
        
        // Currently just making sure that we have around the same value everywhere
        // I want to change this later in some way that makes this correct, but 
        // currently just using the values that have output themselves.
        Double testValue1 = LikelihoodFunctions.BinomLogLik(0.5, 1, 1);
        Assert.assertEquals(testValue1.compareTo(0.6931471805599453), 0);
        
        Double testValue2 = LikelihoodFunctions.BinomLogLik(1, 100, 0);
        Assert.assertEquals(testValue2, -0.0);
        
        Double testValue3 = LikelihoodFunctions.BinomLogLik(1, 100, 1);
        Assert.assertEquals(testValue3, Double.POSITIVE_INFINITY);
        
        Double testValue4 = LikelihoodFunctions.BinomLogLik(0, 1, 1);
        Assert.assertEquals(testValue4, Double.POSITIVE_INFINITY);

        //Some stuff that was created before
        
        
    }
    
    @Test
    public void chiSq1DfToPvalTest(){
        //Check based on wikipedia, and online P value tables:
        
        Double pval = LikelihoodFunctions.determinePvalFrom1DFchiSq(0.004);
        Assert.assertEquals(abs(pval - 0.95) < 0.01, true);
        
        
        pval = LikelihoodFunctions.determinePvalFrom1DFchiSq(0.15);
        Assert.assertEquals(abs(pval - 0.7) < 0.01, true);
        
        
        pval = LikelihoodFunctions.determinePvalFrom1DFchiSq(1.642);
        Assert.assertEquals(abs(pval - 0.2) < 0.0001, true);
        
        
        pval = LikelihoodFunctions.determinePvalFrom1DFchiSq(2.706);
        Assert.assertEquals(abs(pval - 0.1) < 0.0001, true);

        pval = LikelihoodFunctions.determinePvalFrom1DFchiSq(3.841);
        Assert.assertEquals(abs(pval - 0.05) < 0.0001, true);
        
        
        pval = LikelihoodFunctions.determinePvalFrom1DFchiSq(6.635);
        Assert.assertEquals(abs(pval - 0.01) < 0.00001, true);
        
        
        pval = LikelihoodFunctions.determinePvalFrom1DFchiSq(10.828);
        Assert.assertEquals(abs(pval - 0.001) < 0.00001, true);
        
    }
    
    @Test
    public void BetaBinomTest(){
        //Reference values are taken from tests that I did and validated 
        //using the WASP python version.
        
        
        //Make sure that the error and the hetprobability is correct:
        GlobalVariables.hetProb  = 0.980198;
        GlobalVariables.seqError = 0.005;
        
        //These values were checked in WASP implementation.
        
        Assert.assertEquals(abs( LikelihoodFunctions.BetaBinomLogLik(0.132202148437500,1.50728797912598,1.47147178649902,3,1) - 2.7693565626714465) < 0.0000001, true);
        Assert.assertEquals(abs( LikelihoodFunctions.BetaBinomLogLik(0.142944335937500,2.51745270192623,0.934578329324722,0,0) - -0.01960849080585195) < 0.0000001, true);
        Assert.assertEquals(abs( LikelihoodFunctions.BetaBinomLogLik(0.149169921875000,3.21029663085938,2.11465454101563,23,15) - 25.82528716887173) < 0.0000001, true);
        Assert.assertEquals(abs( LikelihoodFunctions.BetaBinomLogLik(0.131591796875000,1.65016174316406,1.67405700683594,12,8) - 13.951339387441848) < 0.0000001, true);
        Assert.assertEquals(abs( LikelihoodFunctions.BetaBinomLogLik(0.149169921875000,2.46604236960411,0.689421594142914,2,1) - 2.0716073059203124) < 0.0000001, true);
        Assert.assertEquals(abs( LikelihoodFunctions.BetaBinomLogLik(0.132202148437500,0.476916499435902,2.34163384512067,3,1) - 5.3287565521228535) < 0.0000001, true);
        Assert.assertEquals(abs( LikelihoodFunctions.BetaBinomLogLik(0.144653320312500,3.52430664002895,3.28067292273045,1,0) - 0.6397056472882476) < 0.0000001, true);
        Assert.assertEquals(abs( LikelihoodFunctions.BetaBinomLogLik(0.141601562500000,0.992187663912773,2.39106918871403,0,0) - -0.01960849080585195) < 0.0000001, true);
        Assert.assertEquals(abs( LikelihoodFunctions.BetaBinomLogLik(0.140258789062500,0.624863505363464,2.00421983003616,1,9) - 3.861381935653935) < 0.0000001, true);
        Assert.assertEquals(abs( LikelihoodFunctions.BetaBinomLogLik(0.142822265625000,2.61502656340599,4.68805944919586,0,2) - 0.8483676462787743) < 0.0000001, true);
        Assert.assertEquals(abs( LikelihoodFunctions.BetaBinomLogLik(0.131591796875000,3.27853605151176,2.28978887200356,7,5) - 8.267991523632132) < 0.0000001, true);
        Assert.assertEquals(abs( LikelihoodFunctions.BetaBinomLogLik(0.144653320312500,4.03051721304655,1.64812137186527,5,2) - 4.2781478286599395) < 0.0000001, true);
        Assert.assertEquals(abs( LikelihoodFunctions.BetaBinomLogLik(0.144653320312500,4.03051721304655,1.64812137186527,5,2) - 4.2781478286599395) < 0.0000001, true);
        Assert.assertEquals(abs( LikelihoodFunctions.BetaBinomLogLik(0.142944335937500,0.606403350830078,2.74511146545410,0,5) - 0.9249243659749231) < 0.0000001, true);
        Assert.assertEquals(abs( LikelihoodFunctions.BetaBinomLogLik(0.131591796875000,3.51683580875397,0.525007486343384,1,1) - 2.2158433296833007) < 0.0000001, true);
        Assert.assertEquals(abs( LikelihoodFunctions.BetaBinomLogLik(0.131591796875000,2.36124551296234,1.02028298377991,2,0) - 0.6909150126333738) < 0.0000001, true);
        Assert.assertEquals(abs( LikelihoodFunctions.BetaBinomLogLik(0.142944335937500,0.619698881171644,1.93650450790301,0,0) - -0.01960849080585195) < 0.0000001, true);
        Assert.assertEquals(abs( LikelihoodFunctions.BetaBinomLogLik(0.170288085937500,1.88455958201702e-09,4.53212775289600,0,8) - 7.786623935533109E-4) < 0.0000001, true);
        Assert.assertEquals(abs( LikelihoodFunctions.BetaBinomLogLik(0.170288085937500,3.18143445253372,2.11250561475754,1,0) - 0.49616756206323065) < 0.0000001, true);
        Assert.assertEquals(abs( LikelihoodFunctions.BetaBinomLogLik(0.142944335937500,1.66619110107422,1.69025421142578,0,0) - -0.01960849080585195) < 0.0000001, true);
        Assert.assertEquals(abs( LikelihoodFunctions.BetaBinomLogLik(0.132202148437500,1.81266784667969,1.42890930175781,5,4) - 6.277860299969595) < 0.00000001, true);
        Assert.assertEquals(abs( LikelihoodFunctions.BetaBinomLogLik(0.141601562500000,1.81266784667969,1.42890930175781,6,3) - 6.012788658310896) < 0.00000001, true);
        Assert.assertEquals(abs( LikelihoodFunctions.BetaBinomLogLik(0.140258789062500,4.76576232910156,3.95823669433594,80,64) - 99.63044991758545) < 0.00000001, true);
        Assert.assertEquals(abs( LikelihoodFunctions.BetaBinomLogLik(0.142822265625000,4.76576232910156,3.95823669433594,64,53) - 81.2212298838738) < 0.00000001, true);
        Assert.assertEquals(abs( LikelihoodFunctions.BetaBinomLogLik(0.144653320312500,4.76576232910156,3.95823669433594,108,89) - 136.4819944170418) < 0.00000001, true);
        Assert.assertEquals(abs( LikelihoodFunctions.BetaBinomLogLik(0.149169921875000,4.76576232910156,3.95823669433594,273,216) - 336.9018295650372) < 0.00000001, true);
        Assert.assertEquals(abs( LikelihoodFunctions.BetaBinomLogLik(0.132202148437500,4.76576232910156,3.95823669433594,116,103) - 152.25781881114017) < 0.00000001, true);
        Assert.assertEquals(abs( LikelihoodFunctions.BetaBinomLogLik(0.140258789062500,2.21351242065430,2.26744461059570,74,90) - 113.7839331418034) < 0.00000001, true);
        Assert.assertEquals(abs( LikelihoodFunctions.BetaBinomLogLik(0.149169921875000,2.21351242065430,2.26744461059570,132,113) - 170.19797394218892) < 0.00000001, true);
        Assert.assertEquals(abs( LikelihoodFunctions.BetaBinomLogLik(0.141601562500000,2.21351242065430,2.26744461059570,31,32) - 44.09746187566005) < 0.00000001, true);
        Assert.assertEquals(abs( LikelihoodFunctions.BetaBinomLogLik(0.170288085937500,3.65502193570137,3.12404105067253,4,1) - 3.2365768631369005) < 0.00000001, true);
        Assert.assertEquals(abs( LikelihoodFunctions.BetaBinomLogLik(0.131591796875000,3.65502193570137,3.12404105067253,2,5) - 5.0913275014018895) < 0.00000001, true);
    
    }

    


}
        
        
