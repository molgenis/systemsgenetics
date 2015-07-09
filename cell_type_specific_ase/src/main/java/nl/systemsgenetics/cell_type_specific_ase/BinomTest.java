/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package nl.systemsgenetics.cell_type_specific_ase;

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
    
    private final double pVal;
    private final double chiSq;
    private final double binomRatio;
    private final double nullLogLik;
    private final double altLogLik;

    
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
        
        //binomial ratio is not determined using Maximum likelihoof algorithm, 
        //just by the actual proportion, can be proven for the binomial distribution
        if(totalRef !=0 ){
            binomRatio = (totalRef * 1.0)  / ( (totalRef + totalAlt) * 1.0 );
        } else{
            binomRatio = 0.0;
        }
        //Null is a ratio set to 0.5 as below.
        nullLogLik = BinomLogLik(0.5, asRef, asAlt);
        //Alt is a ratio based on the binomRatio.
        altLogLik =  BinomLogLik(binomRatio, asRef, asAlt);
        //chi squared statistic is determined based on both null and alt loglikelihoods.
        chiSq = -2.0 * (nullLogLik - altLogLik);
        
        //determine P value based on distribution
        ChiSquaredDistribution distribution = new ChiSquaredDistribution(1);
        pVal = 1 - distribution.cumulativeProbability(chiSq);
    
    }
    
    private double BinomLogLik(double d, ArrayList<Integer> asRef, ArrayList<Integer> asAlt) {
        
        double logLik = 0; 
        
        for(int i = 0; i < asRef.size(); i++){
            
            int total = asRef.get(i) + asAlt.get(i);
            //determine likelihood here.
            BinomialDistribution binomDist = new BinomialDistribution(total, d);
            logLik += log(binomDist.probability(asRef.get(i)));
        
        }
        
        return logLik;
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
