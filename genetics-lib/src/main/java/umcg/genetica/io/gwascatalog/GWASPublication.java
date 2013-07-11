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
public class GWASPublication {
    HashSet<GWASTrait> traits = new HashSet<GWASTrait>();
    HashSet<GWASSNP> snps = new HashSet<GWASSNP>();
    HashSet<GWASLocus> loci = new HashSet<GWASLocus>();
    
    int id;
    String name;
    
}
