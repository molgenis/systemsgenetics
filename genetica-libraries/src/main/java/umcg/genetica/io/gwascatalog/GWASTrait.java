/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package umcg.genetica.io.gwascatalog;

import java.util.ArrayList;
import java.util.HashSet;
import java.util.List;

/**
 *
 * @author harmjan
 */
public class GWASTrait {

    int id;
    String name;
    String cleanName;
    HashSet<GWASLocus> loci = new HashSet<GWASLocus>();
    HashSet<GWASPublication> publishedIn = new HashSet<GWASPublication>();
    HashSet<GWASSNP> snps = new HashSet<GWASSNP>();
    HashSet<GWASSNP> strongestSNPAssociation = new HashSet<GWASSNP>();
    GWASSNP[] snpArray = null;
    HashSet<String> reportedGenes = new HashSet<String>();
    HashSet<String> mappedGenes = new HashSet<String>();

    public String getName() {
        return name;
    }

    public String getCleanName() {
        return cleanName;
    }

    public GWASSNP[] getSNPs() {
        if (snpArray == null) {
            snpArray = snps.toArray(new GWASSNP[0]);
        }
        return snpArray;
    }

    public GWASSNP[] getSNPs(double pThreshold) {
        if (snpArray == null) {
            snpArray = new GWASSNP[snps.size()];
            snps.toArray(snpArray);
        }
        List<GWASSNP> limited = new ArrayList<GWASSNP>();
        for (GWASSNP snp : snpArray) {
            Double p = snp.getPValueAssociatedWithTrait(this);
            if (p != null && p <= pThreshold) {
                limited.add(snp);
            }
        }
        return limited.toArray(new GWASSNP[0]);
    }

    public HashSet<String> getReportedGenes() {
        return reportedGenes;
    }

    public void setReportedGenes(HashSet<String> reportedGenes) {
        this.reportedGenes = reportedGenes;
    }

    public void appendReportedGenes(HashSet<String> reportedGenes) {
        this.reportedGenes.addAll(reportedGenes);
    }

    public HashSet<String> getMappedGenes() {
        return mappedGenes;
    }

    public void setMappedGenes(HashSet<String> mappedGenes) {
        this.mappedGenes = mappedGenes;
    }

    public void appendMappedGenes(HashSet<String> mappedGenes) {
        this.mappedGenes.addAll(reportedGenes);
    }

    @Override
    public String toString() {
        return name;
    }

    public void addTopSNP(GWASSNP gwasTopSNPObj) {
        this.strongestSNPAssociation.add(gwasTopSNPObj);
        this.snps.add(gwasTopSNPObj);
    }
    
    public GWASSNP[] getTopAssociations(){
        return strongestSNPAssociation.toArray(new GWASSNP[0]);
    }
    
    public boolean isTopAssociation(GWASSNP snp){
        return strongestSNPAssociation.contains(snp);
    }
}
