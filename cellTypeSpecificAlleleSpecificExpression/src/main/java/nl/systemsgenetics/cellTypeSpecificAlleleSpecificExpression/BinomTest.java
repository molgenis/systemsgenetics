/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */

package nl.systemsgenetics.cellTypeSpecificAlleleSpecificExpression;

import static java.lang.Math.log;
import java.util.ArrayList;
import org.apache.commons.math3.distribution.BinomialDistribution;
import org.apache.commons.math3.distribution.ChiSquaredDistribution;

/**
 *
 * @author Adriaan van der Graaf 2015
 * This is a class that creates and stores Binomial test statistics,
 * It is created in BinomialTest.java class where it is filled.
 */
public class BinomTest {
    
    private double pVal;
    private double chiSq;
    private double binomRatio;
    private double nullLogLik;
    private double altLogLik;

    
    public BinomTest(ArrayList<Integer> asRef, ArrayList<Integer> asAlt){
        int totalRef;
        int totalAlt;
        
        totalRef = 0;
        totalAlt = 0;
        
        for(int i=0; i < asRef.size(); i++){ 

            //get the totals for the files.
            totalRef += asRef.get(i);
            totalAlt += asAlt.get(i);
        }
        
        int[] asRefArray = new int[asRef.size()];
        int[] asAltArray = new int[asAlt.size()];
        
        for(int j =0; j< asRef.size(); j++){
            
            asRefArray[j] = asRef.get(j);
            asAltArray[j] = asAlt.get(j);
        
        }
        
        //binomial ratio is not determined using Maximum likelihoof algorithm, 
        //just by the actual proportion, can be proven for the binomial distribution
        if(totalRef !=0 ){
            binomRatio = (totalRef * 1.0)  / ( (totalRef + totalAlt) * 1.0 );
        } else{
            binomRatio = 0.0;
        }
        //Null is a ratio set to 0.5 as below.
        nullLogLik = likelihoodFunctions.BinomLogLik(0.5, asRefArray, asAltArray);
        //Alt is a ratio based on the binomRatio.
        altLogLik =  likelihoodFunctions.BinomLogLik(binomRatio, asRefArray, asAltArray);
        //chi squared statistic is determined based on both null and alt loglikelihoods.
        chiSq = likelihoodFunctions.ChiSqFromLogLik(nullLogLik, altLogLik);
        
        
        //determine P value based on distribution

        pVal = likelihoodFunctions.determinePvalFrom1DFchiSq(chiSq);
    
    }
    


    /**
     * @return the pVal
     */
    public double getpVal() {
        return pVal;
    }

    /**
     * @return the chiSq
     */
    public double getChiSq() {
        return chiSq;
    }

    /**
     * @return the binomRatio
     */
    public double getBinomRatio() {
        return binomRatio;
    }

    /**
     * @return the nullLogLik
     */
    public double getNullLogLik() {
        return nullLogLik;
    }

    /**
     * @return the altLogLik
     */
    public double getAltLogLik() {
        return altLogLik;
    }

    
}
