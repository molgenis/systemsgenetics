/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package umcg.genetica.util;

import static org.testng.Assert.*;
import org.testng.annotations.BeforeMethod;
import org.testng.annotations.Test;

/**
 *
 * @author MarcJan
 */
public class RankArrayNGTest {
    
    public RankArrayNGTest() {
    }

    @BeforeMethod
    public void setUpMethod() throws Exception {
    }

    /**
     * Test of rank method, of class RankArray.
     */
    @Test
    public void testRank_doubleArr_boolean() {
        double[] ranks = {4.0d, 2.0d, 0.0d, 7.0d, 5.0d, 8.0d, 9.0d, 1.0d, 6.0d, 3.0d};
        double[] ranks1 = {4.0d, 5.0d, 2.0d, 0.0d, 8.0d, 6.0d, 9.0d, 10.0d, 1.0d, 7.0d, 3.0d};
        
        double[] ranksT = {5.0d, 2.0d, 0.0d, 7.0d, 5.0d, 8.0d, 9.0d, 1.0d, 5.0d, 3.0d};
        double[] ranksT1 = {5.5d, 5.5d, 2.0d, 0.0d, 8.0d, 5.5d, 9.0d, 10.0d, 1.0d, 5.5d, 3.0d};
        
        double[] arrayOfDouble = {5.0d, 3.0d, 0.0d, 7.0d, 5.0d, 8.0d, 9.0d, 1.0d, 5.0d, 4.0d};
        double[] arrayOfDouble1 = {5.0d, 5.0d, 3.0d, 0.0d, 7.0d, 5.0d, 8.0d, 9.0d, 1.0d, 5.0d, 4.0d};
        
        RankArray rda = new RankArray();
        
        double[] valsRanked = rda.rank(arrayOfDouble, false);
        
        boolean correct = true;
        for(int i=0; i < ranks.length; ++i){
//            System.out.println(ranks[i]+"\t"+valsRanked[i]);
            if(ranks[i]!=valsRanked[i]){
                correct = false;
                break;
            }
        }
        assertEquals(correct, true);
        
        double[] valsRanked2 = rda.rank(arrayOfDouble1, false);
        correct = true;
        for(int i=0; i < ranks1.length; ++i){
//            System.out.println(ranks1[i]+"\t"+valsRanked2[i]);
            if(ranks1[i]!=valsRanked2[i]){
                correct = false;
                break;
            }
        }
        assertEquals(correct, true);
        
        double[] valsRankedT = rda.rank(arrayOfDouble, true);
        
        correct = true;
        for(int i=0; i < ranks.length; ++i){
//            System.out.println(ranksT[i]+"\t"+valsRankedT[i]);
            if(ranksT[i]!=valsRankedT[i]){
                correct = false;
//                break;
            }
        }
        assertEquals(correct, true);
        
        double[] valsRankedT2 = rda.rank(arrayOfDouble1, true);
        correct = true;
        for(int i=0; i < ranks1.length; ++i){
//            System.out.println(ranksT1[i]+"\t"+valsRankedT2[i]);
            if(ranksT1[i]!=valsRankedT2[i]){
                correct = false;
                break;
            }
        }
        assertEquals(correct, true);
        
    }

    /**
     * Test of rank method, of class RankArray.
     */
    @Test
    public void testRank_floatArr_boolean() {
    }
}