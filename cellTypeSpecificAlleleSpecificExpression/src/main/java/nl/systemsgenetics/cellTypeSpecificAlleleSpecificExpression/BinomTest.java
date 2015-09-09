/*
 * Copyright (C) 2015 Adriaan van der Graaf
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
        nullLogLik = LikelihoodFunctions.BinomLogLik(0.5, asRefArray, asAltArray);
        //Alt is a ratio based on the binomRatio.
        altLogLik =  LikelihoodFunctions.BinomLogLik(binomRatio, asRefArray, asAltArray);
        //chi squared statistic is determined based on both null and alt loglikelihoods.
        chiSq = LikelihoodFunctions.ChiSqFromLogLik(nullLogLik, altLogLik);
        
        
        //determine P value based on distribution

        pVal = LikelihoodFunctions.determinePvalFrom1DFchiSq(chiSq);
    
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
