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
 * @author Rudi Alberts
 */
public class SimpleGeneticRiskScoreCalculator implements GeneticRiskScoreCalculator {

    private String phenotype;

    private List<Integer> chr = new ArrayList<Integer>();
    private List<Integer> pos = new ArrayList<Integer>();

    private List<String> rsid = new ArrayList<String>();
    private List<String> riskallele = new ArrayList<String>();
    private List<String> otherallele = new ArrayList<String>();

    private List<Double> pvalue = new ArrayList<Double>();
    private List<Double> ororbeta = new ArrayList<Double>();

    // sum and scale betas, 1 or 2 per individual = no. of risk alleles per ind
    @Override
    public TObjectDoubleHashMap<String> calculateRiskScores(RandomAccessGenotypeData genotypeData, String onlyCount, String harmonizedData, double inclusionThreshold) {
        int index;
        int useDosages = 1;
        TObjectDoubleHashMap<String> scores = new TObjectDoubleHashMap<String>(genotypeData.getSampleNames().length); // size of the hashmap

        List<Double> score = new ArrayList<Double>();
        for (index = 0; index < genotypeData.getSampleNames().length; index++) {
            score.add(index, 0.0);
        }

        int nrGCAT = 0;
        int nrnonGCAT = 0;

        int intOnlyCount = 0;
        if (onlyCount == null) {
            System.out.println("onlyCount is null: we add the betas ");
        } else {
            intOnlyCount = 1;
            System.out.println("onlyCount is not null: we only count risk alleles ");
        }

        int intHarmonizedData = 0;
        if (harmonizedData == null) {
            System.out.println("harmonizedData is null: we don't assume harmonized data");
        } else {
            intHarmonizedData = 1;
            System.out.println("harmonizedData is not null: we assume harmonized data");
        }
        System.out.println("risk snp inclusion threshold is:  " + inclusionThreshold);
        
        String complementriskallele = "K";
        
        
        double maxRisk = 0;
        System.out.println("-----------------------------------------------");
        System.out.println("Processing: " + phenotype);
        System.out.println("-----------------------------------------------");
        System.out.println("Size of chr " + chr.size());
        System.out.println("Size of pos " + pos.size());
        System.out.println("Size of ororbeta " + ororbeta.size());

        // if the length of ororbeta is different from the length of chr, then for some snps the effect size is not given
        // leave the whole phenotype out for the moment.
        // loop over risk snps
        for (index = 0; index < chr.size(); index++) {

            Integer currentChr = chr.get(index);
            int currentPos = pos.get(index);
            double currentOrorbeta = ororbeta.get(index);

            System.out.print("SNP" + index + " " + rsid.get(index) + "  " + chr.get(index) + "  " + pos.get(index) + "  " + riskallele.get(index) + " " + otherallele.get(index) + " effectsize " + ororbeta.get(index) + " Pvalue " + pvalue.get(index) + "   ");

            if (pvalue.get(index) < inclusionThreshold) {

                try {
                //snpVariantByPos = inputGenotypes.getSnpVariantByPos("1", 161012760);

                // get alleles etc for all individuals for current snp
                    GeneticVariant snpVariantByPos = genotypeData.getSnpVariantByPos(currentChr.toString(), currentPos);

                    System.out.print("         Chr: " + snpVariantByPos.getSequenceName() + " Pos " + snpVariantByPos.getStartPos()
                            + " MINORALLELE  " + snpVariantByPos.getMinorAllele()
                            + "  ALLALLELES  " + snpVariantByPos.getVariantAlleles()
                            + "  FIRST: " + snpVariantByPos.getVariantAlleles().get(0)
                            + "  SECOND: " + snpVariantByPos.getVariantAlleles().get(1)
                            + "  Is GC AT snp: " + snpVariantByPos.getVariantAlleles().isAtOrGcSnp()
                            + "  MAF  " + snpVariantByPos.getMinorAlleleFrequency()
                            + " ");

                    if (snpVariantByPos.getVariantAlleles().isAtOrGcSnp()) {
                        nrGCAT++;
                    } else {
                        nrnonGCAT++;
                    }

                    List<Alleles> alleles = snpVariantByPos.getSampleVariants();
                    float[] dosages = snpVariantByPos.getSampleDosages();
                    List<String> myids = snpVariantByPos.getAllIds();

                    String originalriskallele = riskallele.get(index);
                    String usedriskallele = riskallele.get(index);
                    String secondallele = otherallele.get(index);

                    if (ororbeta.get(index) < 0) {
                        usedriskallele = otherallele.get(index);

                    }

                    System.out.print(myids.get(0) + " ORIGRISKALLELE " + originalriskallele + " ");
                    System.out.print(" USEDRISKALLELE " + usedriskallele + " ");

                    //System.out.print("alleles size: " + alleles.size());
                    //System.out.print("dosages size: " + dosages.length);

                    // do not use the SNP if intHarmonized == 0 AND it is a GC/AT snp
                    // that means, do this if intHarmonzed == 1 OR GC/AT is false

                    if (intHarmonizedData==1 | snpVariantByPos.getVariantAlleles().isAtOrGcSnp()==false) {
                        int index2;

                        for (index2 = 0; index2 < alleles.size(); index2++) {
                            if (index2 < 20) {
                                System.out.print(" -- allele " + alleles.get(index2).getAllelesAsString() + "__" + alleles.get(index2).getAllelesAsString().get(0) + "__" + alleles.get(index2).getAllelesAsString().get(1));
                                System.out.print(" dosage: " + dosages[index2]);
                            }

                            if (intOnlyCount == 0) {

                                // used risk allel is gelijk aan FIRST allel, neem dosage
                                if (usedriskallele.equals(snpVariantByPos.getVariantAlleles().get(0).toString())) {
                                    score.set(index2, score.get(index2) + dosages[index2] * Math.abs(ororbeta.get(index)));
                                    if (index2 < 20) {
                                        System.out.print(" added (first allele, orig strand) " + dosages[index2] * Math.abs(ororbeta.get(index)) + "**" + alleles.get(index2).getAllelesAsString().get(0) + "**");
                                    }
                                }

                                // needed if risk allel is reported as second allele, we are still looking at the same strand
                                // used risk allel is gelijk aan SECOND allel, neem 2-dosage
                                if (usedriskallele.equals(snpVariantByPos.getVariantAlleles().get(1).toString())) {
                                    score.set(index2, score.get(index2) + (2 - dosages[index2]) * Math.abs(ororbeta.get(index)));
                                    if (index2 < 20) {
                                        System.out.print(" added (second allele, orig strand) " + (2 - dosages[index2]) * Math.abs(ororbeta.get(index)) + "**" + alleles.get(index2).getAllelesAsString().get(0) + "**");
                                    }
                                }

                                // now we will look at the other strand 
                                // to do so we keep the genotype the same and take the reverse of the risk allele
                                // only do if intHarmonizedData==0
                                if (intHarmonizedData == 0) {
                                    if (usedriskallele.equals("C")) {
                                        complementriskallele = "G";
                                    }
                                    if (usedriskallele.equals("G")) {
                                        complementriskallele = "C";
                                    }
                                    if (usedriskallele.equals("A")) {
                                        complementriskallele = "T";
                                    }
                                    if (usedriskallele.equals("T")) {
                                        complementriskallele = "A";
                                    }

                                    // used risk allel is gelijk aan FIRST allel, neem dosage
                                    if (complementriskallele.equals(snpVariantByPos.getVariantAlleles().get(0).toString())) {
                                        score.set(index2, score.get(index2) + dosages[index2] * Math.abs(ororbeta.get(index)));
                                        if (index2 < 20) {
                                            System.out.print(" added (first allele, other strand) " + dosages[index2] * Math.abs(ororbeta.get(index)) + "**" + alleles.get(index2).getAllelesAsString().get(0) + "**");
                                        }
                                    }

                                    // needed if risk allel is reported as second allele, we are still looking at the same strand
                                    // used risk allel is gelijk aan SECOND allel, neem 2-dosage
                                    if (complementriskallele.equals(snpVariantByPos.getVariantAlleles().get(1).toString())) {
                                        score.set(index2, score.get(index2) + (2 - dosages[index2]) * Math.abs(ororbeta.get(index)));
                                        if (index2 < 20) {
                                            System.out.print(" added (second allele,other strand) " + (2 - dosages[index2]) * Math.abs(ororbeta.get(index)) + "**" + alleles.get(index2).getAllelesAsString().get(0) + "**");
                                        }
                                    }

                                }

                            } else { // here only count risk alleles

                                // used risk allel is gelijk aan FIRST allel, neem dosage
                                if (usedriskallele.equals(snpVariantByPos.getVariantAlleles().get(0).toString())) {
                                    score.set(index2, score.get(index2) + dosages[index2]);
                                    if (index2 < 20) {
                                        System.out.print(" added (first allele, orig strand) " + dosages[index2] + "**" + alleles.get(index2).getAllelesAsString().get(0) + "**");
                                    }
                                }

                                // used risk allel is gelijk aan SECOND allel, neem 2-dosage
                                if (usedriskallele.equals(snpVariantByPos.getVariantAlleles().get(1).toString())) {
                                    score.set(index2, score.get(index2) + (2 - dosages[index2]));
                                    if (index2 < 20) {
                                        System.out.print(" added (second allele, orig strand) " + (2 - dosages[index2]) + "**" + alleles.get(index2).getAllelesAsString().get(0) + "**");
                                    }
                                }

                                // only do if intHarmonizedData==0
                                if (intHarmonizedData == 0) {

                                    if (usedriskallele.equals("C")) {
                                        complementriskallele = "G";
                                    }
                                    if (usedriskallele.equals("G")) {
                                        complementriskallele = "C";
                                    }
                                    if (usedriskallele.equals("A")) {
                                        complementriskallele = "T";
                                    }
                                    if (usedriskallele.equals("T")) {
                                        complementriskallele = "A";
                                    }

                                    // used risk allel is gelijk aan FIRST allel, neem dosage
                                    if (complementriskallele.equals(snpVariantByPos.getVariantAlleles().get(0).toString())) {
                                        score.set(index2, score.get(index2) + dosages[index2]);
                                        if (index2 < 20) {
                                            System.out.print(" added (first allele, other strand) " + dosages[index2] + "**" + alleles.get(index2).getAllelesAsString().get(0) + "**");
                                        }
                                    }

                                    // used risk allel is gelijk aan SECOND allel, neem 2-dosage
                                    if (complementriskallele.equals(snpVariantByPos.getVariantAlleles().get(1).toString())) {
                                        score.set(index2, score.get(index2) + (2 - dosages[index2]));
                                        if (index2 < 20) {
                                            System.out.print(" added (second allele, other strand) " + (2 - dosages[index2]) + "**" + alleles.get(index2).getAllelesAsString().get(0) + "**");
                                        }
                                    }

                                }
                            }
                        }
                    } else {
                        System.out.print("      SNP excluded because non-harmonized data and SNP is a GC/AT snp");                       
                    }
                // if genotype data is present, add two times beta to max risk

                //maxRisk = maxRisk + 2 * currentOrorbeta;
                } catch (NullPointerException ex) {
                    System.out.print("       No genotype data for this genomic position");

                }
            } else {
                // threshold not met
                System.out.print("      SNP excluded because Pvalue not low enough");

            }

            System.out.println();
        }

        //System.out.println("MAX RISK: " + maxRisk);
        System.out.println("IN CALCULATOR FIRST SCORE: " + score.get(0));
        System.out.println("IN CALCULATOR SECND SCORE: " + score.get(1));

        // with scaling
        //for (index = 0; index < genotypeData.getSampleNames().length; index++) {
        //    scores.put(genotypeData.getSampleNames()[index], score.get(index) / maxRisk );  // here score is scaled between 0 and 1 with 1 corresponding to maxrisk
        //}
        // without scaling
        for (index = 0; index < genotypeData.getSampleNames().length; index++) {
            scores.put(genotypeData.getSampleNames()[index], score.get(index));  // here score is scaled between 0 and 1 with 1 corresponding to maxrisk
        }

        return (scores);
    }

    @Override
    public TObjectDoubleHashMap<String> calculateRiskScores(RandomAccessGenotypeData genotypeData, PhenotypeData phenotypeData
    ) {
        throw new UnsupportedOperationException("Not supported yet.");
    }

    @Override
    public String getPhenotype() {
        return this.phenotype;
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

    public void addOtherallele(String riskallele) {
        this.otherallele.add(riskallele);
    }

    public void addOrorbeta(Double ororbeta) {
        this.ororbeta.add(ororbeta);
    }

    public void addPvalue(Double pvalue) {
        this.pvalue.add(pvalue);
    }
}
