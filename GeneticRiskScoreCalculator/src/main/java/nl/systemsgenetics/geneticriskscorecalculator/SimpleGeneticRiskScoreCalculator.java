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
    public TObjectDoubleHashMap<String> calculateRiskScores(RandomAccessGenotypeData genotypeData, String onlyCount, double inclusionThreshold) {
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
        System.out.println("inclusion threshold is:  " + inclusionThreshold);

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
                // reverse risk allele if not found in data (on other strand)
                    //if (!snpVariantByPos.getVariantAlleles().isAtOrGcSnp()) {
                    //    if (!originalriskallele.equals(snpVariantByPos.getVariantAlleles().get(0).toString())) {
                    //        if (!originalriskallele.equals(snpVariantByPos.getVariantAlleles().get(1).toString())) {

                //            if (originalriskallele.equals("C")) {
                    //                usedriskallele = "G";
                    //            }
                    //            if (originalriskallele.equals("G")) {
                    //                usedriskallele = "C";
                    //            }
                    //            if (originalriskallele.equals("A")) {
                    //                usedriskallele = "T";
                    //            }
                    //            if (originalriskallele.equals("T")) {
                    //                usedriskallele = "A";
                    //            }
                    //        }
                    //    }
                    //}
                    System.out.print(myids.get(0) + " ORIGRISKALLELE " + originalriskallele + " ");
                    System.out.print(" USEDRISKALLELE " + usedriskallele + " ");

                //System.out.print("alleles size: " + alleles.size());
                    //System.out.print("dosages size: " + dosages.length);
                    if (useDosages == 0) {
                        // probably not needed anymore. for hard genotypes dosages can also be used (0,1,2)
                        // if needed, update to usedriskallele instead of riskallele
                        int index2 = 0;
                        //loop over all alleles (individuals) for this snp
                            
                        for (Alleles allele : alleles) {
                            if (index2 < 20) {
                                System.out.print(allele.toString() + " ");  // only alleles
                            }

                            // add beta to score max. two times per individual
                            if (snpVariantByPos.getVariantAlleles().isAtOrGcSnp() == false && (riskallele.get(index).equals(allele.getAllelesAsString().get(0)))) {
                                score.set(index2, score.get(index2) + currentOrorbeta);
                            }
                            if (snpVariantByPos.getVariantAlleles().isAtOrGcSnp() == false && (riskallele.get(index).equals(allele.getAllelesAsString().get(1)))) {
                                score.set(index2, score.get(index2) + currentOrorbeta);
                            }
                            index2++;
                        }
                    }

                    // float[] dosages = snpVariantByPos.getSampleDosages();
                    if (useDosages == 1) {
                        int index2 = 0;

                        for (index2 = 0; index2 < alleles.size(); index2++) {
                            if (index2 < 20) {
                                System.out.print(" -- allele " + alleles.get(index2).getAllelesAsString() + "__" + alleles.get(index2).getAllelesAsString().get(0) + "__" + alleles.get(index2).getAllelesAsString().get(1));
                                System.out.print(" dosage: " + dosages[index2]);
                            }
//                        if (snpVariantByPos.getVariantAlleles().isAtOrGcSnp() == false && (riskallele.get(index).equals(allele.getAllelesAsString().get(0)) || riskallele.get(index).equals(allele.getAllelesAsString().get(1)))) {
//                            score.set(index2, score.get(index2) + 1);
//                        }
                            // if usedallele is equal to first allele, add dosage to score
                            // if usedallele is equal to seond allele, add 2-dosage to score

                        // this is without GC AT snps and counting the dosage
                        //if (snpVariantByPos.getVariantAlleles().isAtOrGcSnp() == false && usedriskallele.equals(alleles.get(index2).getAllelesAsString().get(0)) ) {
                            //    score.set(index2, score.get(index2) + dosages[index2]);
                            //}
                            //if (snpVariantByPos.getVariantAlleles().isAtOrGcSnp() == false && usedriskallele.equals(alleles.get(index2).getAllelesAsString().get(1)) ) {
                            //    score.set(index2, score.get(index2) + 2 - dosages[index2]);
                            //}               
                        // this is for all snps, and if beta is negative, we take the other allele as risk allele and use abs(beta), and counting the beta's
                            // dosage geeft al nul of 1 of 2 voor het risk allel
                            // dus 2 chromosomen worden al gecheckt.
                            // als je het otherallele als risk allel neemt omdat beta negatief moet je 2-dosage nemen
                        // het eerste allel is niet altijd het minor allel
                            // dosage wordt altijd gegeven voor het eerste allel
                            // als je de dosage van "C" wilt weten moet je checken of "C" het eerste allel is, dan dosage nemen
                            // of of "C" het tweede allel is, dan 2-dosage nemen
                            if (intOnlyCount == 0) {

                                // used risk allel is gelijk aan FIRST allel, neem dosage
                                if (usedriskallele.equals(snpVariantByPos.getVariantAlleles().get(0).toString())) {
                                    score.set(index2, score.get(index2) + dosages[index2] * Math.abs(ororbeta.get(index)));
                                    if (index2 < 20) {
                                        System.out.print(" added " + dosages[index2] * Math.abs(ororbeta.get(index)) + "**" + alleles.get(index2).getAllelesAsString().get(0) + "**");
                                    }
                                }
                                // used risk allel is gelijk aan SECOND allel, neem 2-dosage
                                if (usedriskallele.equals(snpVariantByPos.getVariantAlleles().get(1).toString())) {
                                    score.set(index2, score.get(index2) + (2 - dosages[index2]) * Math.abs(ororbeta.get(index)));
                                    if (index2 < 20) {
                                        System.out.print(" added " + (2 - dosages[index2]) * Math.abs(ororbeta.get(index)) + "**" + alleles.get(index2).getAllelesAsString().get(0) + "**");
                                    }
                                }

                            } else { // here only count risk alleles

                                // used risk allel is gelijk aan FIRST allel, neem dosage
                                if (usedriskallele.equals(snpVariantByPos.getVariantAlleles().get(0).toString())) {
                                    score.set(index2, score.get(index2) + dosages[index2]);
                                    if (index2 < 20) {
                                        System.out.print(" added " + dosages[index2] + "**" + alleles.get(index2).getAllelesAsString().get(0) + "**");
                                    }
                                }
                                // used risk allel is gelijk aan SECOND allel, neem 2-dosage
                                if (usedriskallele.equals(snpVariantByPos.getVariantAlleles().get(1).toString())) {
                                    score.set(index2, score.get(index2) + (2 - dosages[index2]));
                                    if (index2 < 20) {
                                        System.out.print(" added " + (2 - dosages[index2]) + "**" + alleles.get(index2).getAllelesAsString().get(0) + "**");
                                    }
                                }

                            }

                        }
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

        System.out.println("MAX RISK: " + maxRisk);
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
