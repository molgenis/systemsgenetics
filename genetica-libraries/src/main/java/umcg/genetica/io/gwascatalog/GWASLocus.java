/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package umcg.genetica.io.gwascatalog;

import java.util.HashSet;

/**
 *
 * @author harmjan
 */
public class GWASLocus {
    int id; 
    int start;
    int end;
    byte chr;
    HashSet<GWASSNP> snps = new HashSet<GWASSNP>();
    
	    
}
