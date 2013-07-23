/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package umcg.genetica.io.gwascatalog;

import java.util.HashMap;
import java.util.HashSet;
import umcg.genetica.containers.Pair;

/**
 *
 * @author harmjan
 */
public class GWASPublication {

    HashSet<GWASTrait> traits = new HashSet<GWASTrait>();
    HashSet<GWASSNP> snps = new HashSet<GWASSNP>();
    HashSet<GWASLocus> loci = new HashSet<GWASLocus>();
    HashMap<Pair<GWASTrait, GWASSNP>, Double> assoc = new HashMap<Pair<GWASTrait, GWASSNP>, Double>();
    int id;
    String name;

    public void setPValueAssociatedWithTrait(GWASSNP gwasSNPObj, GWASTrait traitObj, Double pval) {
        assoc.put(new Pair<GWASTrait, GWASSNP>(traitObj, gwasSNPObj), pval);
    }

    public Double getPValueAssociatedWithTrait(GWASSNP gwasSNPObj, GWASTrait traitObj) {
        return assoc.get(new Pair<GWASTrait, GWASSNP>(traitObj, gwasSNPObj));
    }
}
