package nl.systemsgenetics.geneticriskscorecalculator;

import gnu.trove.map.hash.TObjectDoubleHashMap;
import java.util.ArrayList;
import java.util.List;
import java.util.regex.Pattern;
import org.molgenis.genotype.Alleles;
import org.molgenis.genotype.RandomAccessGenotypeData;
import org.molgenis.genotype.Alleles;
import org.molgenis.genotype.GenotypeDataException;
import org.molgenis.genotype.RandomAccessGenotypeData;
import org.molgenis.genotype.multipart.IncompatibleMultiPartGenotypeDataException;
import org.molgenis.genotype.variant.GeneticVariant;
/**
 *
 * @author Patrick Deelen
 */
public class SimpleGeneticRiskScoreCalculator implements GeneticRiskScoreCalculator {

    public String phenotype;

    public List<Integer> chr = new ArrayList<Integer>();
    public List<Integer> pos = new ArrayList<Integer>();

    public List<String> rsid = new ArrayList<String>();
    public List<String> riskallele = new ArrayList<String>();

    @Override
    public TObjectDoubleHashMap<String> calculateRiskScores(RandomAccessGenotypeData genotypeData) {
        int index;
        TObjectDoubleHashMap<String> scores = new TObjectDoubleHashMap<String>(genotypeData.getSampleNames().length); // size of the hashmap

        List<Integer> score = new ArrayList<Integer>();
        for (index = 0; index < genotypeData.getSampleNames().length; index++) {
            score.add(0);
        }

        int nrGCAT = 0;
        int nrnonGCAT = 0;

        for (index = 0; index < chr.size(); index++) {

            Integer currentChr = chr.get(index);
            int currentPos = Integer.valueOf(pos.get(index)); // necessary?

            System.out.print("SNP" + index + " " + rsid.get(index) + "  " + chr.get(index) + "  " + pos.get(index) + "  " + riskallele.get(index) + "  ");

            try {
                //snpVariantByPos = inputGenotypes.getSnpVariantByPos("1", 161012760);
                GeneticVariant snpVariantByPos = genotypeData.getSnpVariantByPos(currentChr.toString(), currentPos);

                System.out.print("MINORALLELE  " + snpVariantByPos.getMinorAllele() + "  ALLALLELES  " + snpVariantByPos.getVariantAlleles() + "  Is GC AT snp: " + snpVariantByPos.getVariantAlleles().isAtOrGcSnp() + "  MAF  "
                        + snpVariantByPos.getMinorAlleleFrequency() + " ");

                if (snpVariantByPos.getVariantAlleles().isAtOrGcSnp()) {
                    nrGCAT++;
                } else {
                    nrnonGCAT++;
                }

                List<Alleles> alleles = snpVariantByPos.getSampleVariants();
                float[] dosages = snpVariantByPos.getSampleDosages();
                List<String> myids = snpVariantByPos.getAllIds();

                System.out.print(myids.get(0) + " ");
                int index2 = 0;
                for (Alleles allele : alleles) {
                    if (index2 < 20) {
                        System.out.print(allele.toString() + " ");  // only alleles
                    }
                    if (snpVariantByPos.getVariantAlleles().isAtOrGcSnp() == false && (riskallele.get(index).equals(allele.getAllelesAsString().get(0)) || riskallele.get(index).equals(allele.getAllelesAsString().get(1)))) {
                        score.set(index2, score.get(index2) + 1);
                    }
                    index2++;
                }
                System.out.println();

            } catch (NullPointerException ex) {
                System.out.println("null pointer exception");

            }

        }

        for (index = 0; index < genotypeData.getSampleNames().length; index++) {
            scores.put(genotypeData.getSampleNames()[index], score.get(index));
        }

        return (scores);
    }

    @Override
    public TObjectDoubleHashMap<String> calculateRiskScores(RandomAccessGenotypeData genotypeData, PhenotypeData phenotypeData) {
        throw new UnsupportedOperationException("Not supported yet.");
    }

    @Override
    public String getPhenotype() {
        return (this.phenotype);
    }

    public void setPhenotype(String phenotype) {
        this.phenotype = phenotype;
    }

    public void addChr(Integer chr) {
        this.chr.add(chr);
    }

    public void addPos(Integer pos) {
        this.pos.add(pos);
    }

    public void addRsid(String rsid) {
        this.rsid.add(rsid);
    }

    public void addRiskallele(String riskallele) {
        this.riskallele.add(riskallele);
    }

}
