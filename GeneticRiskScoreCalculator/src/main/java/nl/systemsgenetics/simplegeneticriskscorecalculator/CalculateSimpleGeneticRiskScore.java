/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package nl.systemsgenetics.simplegeneticriskscorecalculator;

import gnu.trove.map.hash.THashMap;
import java.io.File;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.HashSet;
import java.util.Map.Entry;
import java.util.logging.Level;
import java.util.logging.Logger;
import org.molgenis.genotype.RandomAccessGenotypeData;
import static org.molgenis.genotype.util.LdCalculator.calculateRsquare;
import org.molgenis.genotype.util.LdCalculatorException;
import org.molgenis.genotype.variant.GeneticVariant;
import umcg.genetica.console.ProgressBar;
import umcg.genetica.io.text.TextFile;
import umcg.genetica.math.matrix2.DoubleMatrixDataset;

/**
 *
 * @author MarcJan
 */
class CalculateSimpleGeneticRiskScore {

    private static final String[] chrOrder = {"1", "2", "3", "4", "5", "6", "7", "8", "9", "10", "11", "12", "13", "14", "15", "16", "17", "18", "19", "20", "21", "22"};

    static DoubleMatrixDataset<String, String> calculate(RandomAccessGenotypeData genotypeData, THashMap<String, THashMap<String, THashMap<String, ArrayList<RiskEntry>>>> risks, File outputFolder, double rSquare, double windowSize, boolean debugMode, double[] pValueThreshold, boolean sumRisk) {
        ArrayList<String> keys = new ArrayList<String>();
        for (Entry<String, THashMap<String, THashMap<String, ArrayList<RiskEntry>>>> riskScorePheno : risks.entrySet()) {
            for (Entry<String, THashMap<String, ArrayList<RiskEntry>>> riskScorePheno2 : riskScorePheno.getValue().entrySet()) {
                keys.add(riskScorePheno.getKey() + riskScorePheno2.getKey());
            }
        }

        DoubleMatrixDataset<String, String> scores = new DoubleMatrixDataset<String, String>(keys, Arrays.asList(genotypeData.getSampleNames()));

        ProgressBar p = new ProgressBar(risks.size() * chrOrder.length);

        for (int counter = 0; counter < chrOrder.length; counter++) {
            for (Entry<String, THashMap<String, THashMap<String, ArrayList<RiskEntry>>>> riskScorePheno : risks.entrySet()) {
                HashSet<String> excludeList = new HashSet<String>();
                for (double pVal : pValueThreshold){
                    String key = "_P" + pVal;
                    THashMap<String, ArrayList<RiskEntry>> riskScorePheno2 = riskScorePheno.getValue().get(key);
                    String NameOfEntry = riskScorePheno.getKey() + key;
                    int rowNr = scores.getHashRows().get(NameOfEntry);
                    try {
                        TextFile out = null;
                        if (debugMode) {
                            System.out.println(NameOfEntry);

                            out = new TextFile(outputFolder + File.separator + NameOfEntry + "Chr" + chrOrder[counter] + ".log", TextFile.W);

                            out.write("SNPs used for GRS calculation:\n");
                        }
                        int nrSNPs = 0;

//                        System.out.println("Processing chromosome:\t" + chrOrder[counter]);
                        if (riskScorePheno2.containsKey(chrOrder[counter])) {

                            ArrayList<RiskEntry> valueE2 = riskScorePheno2.get(chrOrder[counter]);

                            int nrSNPsThisChr = valueE2.size();
                            boolean[] excludeSNPs = new boolean[nrSNPsThisChr];

                            //Get the original entries back, so we are sure we dont need to do to many look ups.
                            if (excludeList.size() > 0) {
                                for (int snp = 0; snp < nrSNPsThisChr; snp++) {
                                    if (excludeList.contains(valueE2.get(snp).getRsName())) {
                                        excludeSNPs[snp] = true;
                                    }
                                }
                            }

                            //Actual scoring.
                            for (int snp = 0; snp < nrSNPsThisChr; snp++) {
                                if (!excludeSNPs[snp]) {
                                    RiskEntry riskE = valueE2.get(snp);
                                    //System.out.println(snpID + "\t" + c + "\t" + chrPos + "\t" + object.doubleValue);
                                    GeneticVariant var1 = genotypeData.getSnpVariantByPos(riskE.getChr(), riskE.getPos());
                                    //Check if at least 75% of the sampels have information for the SNP otherwise it is removed by default.
                                    if (var1.getCallRate() < 0.75) {
                                        excludeSNPs[snp] = true;
                                        excludeList.add(riskE.getRsName());
                                        continue;
                                    }
                                    
                                    if (debugMode) {
                                        out.write(riskE.InfoToString() + "\n");
                                    }
                                    
                                    double or = riskE.getOr();
                                    boolean riskCodedAsTwo;
                                    if(sumRisk && or <0){
                                        or = or*-1;
                                        riskCodedAsTwo = false;
                                        if (!(riskE.getAllele()==(var1.getRefAllele().getAlleleAsSnp()) || riskE.getAllele()==(var1.getRefAllele().getComplement().getAlleleAsSnp()))) {
                                            riskCodedAsTwo = true;
                                        }
                                    } else {
                                        riskCodedAsTwo = true;
                                        if (!(riskE.getAllele()==(var1.getRefAllele().getAlleleAsSnp()) || riskE.getAllele()==(var1.getRefAllele().getComplement().getAlleleAsSnp()))) {
                                            riskCodedAsTwo = false;
                                        }
                                    }
//                                    StringBuilder Genos = new StringBuilder();
//                                    StringBuilder Genos1 = new StringBuilder();
//                                    StringBuilder Genos2 = new StringBuilder();
//                                    StringBuilder Genos3 = new StringBuilder();
                                    for (int sample = 0; sample < var1.getSampleCalledDosages().length; sample++) {
                                        if (var1.getSampleCalledDosages()[sample] != -1) {
                                            if(riskCodedAsTwo){
                                                scores.getMatrix().setQuick(rowNr, sample, (scores.getMatrix().getQuick(rowNr, sample) + (or * var1.getSampleCalledDosages()[sample])));
//                                                Genos2.append(or * var1.getSampleCalledDosages()[sample]).append(",");
                                            } else {
                                                scores.getMatrix().setQuick(rowNr, sample, (scores.getMatrix().getQuick(rowNr, sample) + (or * Math.abs(var1.getSampleCalledDosages()[sample] - 2))));
//                                                Genos2.append(or * Math.abs(var1.getSampleCalledDosages()[sample]-2)).append(",");
                                            }
//                                            Genos.append(var1.getSampleCalledDosages()[sample]).append(",");
//                                            Genos1.append(var1.getSampleVariants().get(sample).toString());Genos1.append(",");
//                                            Genos3.append(scores.getMatrix().getQuick(rowNr, sample)).append(",");
                                        }
                                    }
                                    //I have now removed the conversion from 0 to 2 and vice versa from the calculations. Why is this needed?
//                                            System.out.println("");
//                                            System.out.println("SNP: "+riskE.getRsName());
//                                            System.out.println("Allele in data: "+var1.getRefAllele().toString());
//                                            System.out.println("Allele: "+riskE.getAllele());
//                                            System.out.println("Diction: "+direction);
//                                            System.out.println("Or: "+or);
//                                            System.out.println("Genotypes: "+Genos1.toString());
//                                            System.out.println("Dosages: "+Genos.toString());
//                                            System.out.println("Scores: "+Genos2.toString());
//                                            System.out.println("ScoresT: "+Genos3.toString());
//                                            System.out.println("");

                                    nrSNPs++;

                                    for (int t = snp + 1; t < nrSNPsThisChr; t++) {
                                        if (!excludeSNPs[t]) {
                                            RiskEntry riskE2 = valueE2.get(t);
                                            if (Math.abs(riskE2.getPos() - riskE.getPos()) <= windowSize) {
                                                GeneticVariant var2 = genotypeData.getSnpVariantByPos(riskE2.getChr(), riskE2.getPos());
                                                if (var2.getCallRate() < 0.75) {
                                                    excludeSNPs[t] = true;
                                                    excludeList.add(riskE2.getRsName());
                                                    continue;
                                                }
                                                try {
                                                    if (calculateRsquare(var1, var2) >= rSquare) {
                                                        excludeSNPs[t] = true;
                                                        excludeList.add(riskE2.getRsName());
                                                    }
                                                } catch (LdCalculatorException ex) {
                                                    Logger.getLogger(CalculateSimpleGeneticRiskScore.class.getName()).log(Level.SEVERE, null, ex);
                                                }

                                            }
                                        }
                                    }
                                }
                            }
                        }
                        if (debugMode) {
                            out.write("Total SNPs used: " + nrSNPs);
                            out.close();
                        }
                    } catch (IOException ex) {
                        Logger.getLogger(CalculateSimpleGeneticRiskScore.class.getName()).log(Level.SEVERE, null, ex);
                    }
                }
                p.iterate();
            }
        }
        p.close();

        return scores;
    }

    static DoubleMatrixDataset<String, String> calculateTwoStages(RandomAccessGenotypeData genotypeData, THashMap<String, THashMap<String, THashMap<String, ArrayList<RiskEntry>>>> risks, File outputFolder, double rSquare, double[] windowSize, boolean debugMode, double[] pValueThreshold, boolean sumRisk) {
        ArrayList<String> keys = new ArrayList<String>();
        for (Entry<String, THashMap<String, THashMap<String, ArrayList<RiskEntry>>>> riskScorePheno : risks.entrySet()) {
            for (Entry<String, THashMap<String, ArrayList<RiskEntry>>> riskScorePheno2 : riskScorePheno.getValue().entrySet()) {
                keys.add(riskScorePheno.getKey() + riskScorePheno2.getKey());
            }
        }

        DoubleMatrixDataset<String, String> scores = new DoubleMatrixDataset<String, String>(keys, Arrays.asList(genotypeData.getSampleNames()));
        ProgressBar p = new ProgressBar(risks.size() * chrOrder.length);

        for (int counter = 0; counter < chrOrder.length; counter++) {
            for (Entry<String, THashMap<String, THashMap<String, ArrayList<RiskEntry>>>> riskScorePheno : risks.entrySet()) {
                HashSet<String> excludeList = new HashSet<String>();

                for (double pVal : pValueThreshold){
                    String key = "_P" + pVal;
                    THashMap<String, ArrayList<RiskEntry>> riskScorePheno2 = riskScorePheno.getValue().get(key);
                    try {
                        String NameOfEntry = riskScorePheno.getKey() + key;
                        int rowNr = scores.getHashRows().get(NameOfEntry);

                        TextFile out = null;
                        if (debugMode) {
                            System.out.println(NameOfEntry);

                            out = new TextFile(outputFolder + File.separator + NameOfEntry + "Chr" + chrOrder[counter] + ".log", TextFile.W);

                            out.write("SNPs used for GRS calculation:\n");
                        }

                        int nrSNPs = 0;

//                        System.out.println("Processing chromosome:\t" + chrOrder[counter]);
                        if (riskScorePheno2.containsKey(chrOrder[counter])) {
                            ArrayList<RiskEntry> valueE2 = riskScorePheno2.get(chrOrder[counter]);

                            int nrSNPsThisChr = valueE2.size();
                            boolean[] excludeSNPs = new boolean[nrSNPsThisChr];

                            //Get the original entries back, so we are sure we dont need to do to many look ups.
                            if (excludeList.size() > 0) {
                                for (int snp = 0; snp < nrSNPsThisChr; snp++) {
                                    if (excludeList.contains(valueE2.get(snp).getRsName())) {
                                        excludeSNPs[snp] = true;
                                    }
                                }
                            }
                            //Loop 1, pre-filtering.
                            for (int snp = 0; snp < nrSNPsThisChr; snp++) {
                                if (!excludeSNPs[snp]) {
                                    RiskEntry riskE = valueE2.get(snp);
                                    //System.out.println(snpID + "\t" + c + "\t" + chrPos + "\t" + object.doubleValue);
                                    GeneticVariant var1 = genotypeData.getSnpVariantByPos(riskE.getChr(), riskE.getPos());
                                    //Check if at least 75% of the sampels have information for the SNP otherwise it is removed by default.
                                    if (var1.getCallRate() < 0.75) {
                                        excludeSNPs[snp] = true;
                                        excludeList.add(riskE.getRsName());
                                        continue;
                                    }

                                    for (int t = snp + 1; t < nrSNPsThisChr; t++) {
                                        if (!excludeSNPs[t]) {
                                            RiskEntry riskE2 = valueE2.get(t);
                                            if (Math.abs(riskE2.getPos() - riskE.getPos()) <= windowSize[0]) {
                                                GeneticVariant var2 = genotypeData.getSnpVariantByPos(riskE2.getChr(), riskE2.getPos());
                                                if (var2.getCallRate() < 0.75) {
                                                    excludeSNPs[t] = true;
                                                    excludeList.add(riskE2.getRsName());
                                                    continue;
                                                }
                                                try {
                                                    if (calculateRsquare(var1, var2) >= rSquare) {
                                                        excludeSNPs[t] = true;
                                                        excludeList.add(riskE2.getRsName());
                                                    }
                                                } catch (LdCalculatorException ex) {
                                                    Logger.getLogger(CalculateSimpleGeneticRiskScore.class.getName()).log(Level.SEVERE, null, ex);
                                                }
                                            }
                                        }
                                    }
                                }
                            }

                            //Loop 2, Actual scoring.
                            for (int snp = 0; snp < nrSNPsThisChr; snp++) {
                                if (!excludeSNPs[snp]) {
                                    RiskEntry riskE = valueE2.get(snp);
                                    //System.out.println(snpID + "\t" + c + "\t" + chrPos + "\t" + object.doubleValue);
                                    GeneticVariant var1 = genotypeData.getSnpVariantByPos(riskE.getChr(), riskE.getPos());
                                    if (debugMode) {
                                        out.write(riskE.InfoToString() + "\n");
                                    }
                                    
                                    double or = riskE.getOr();
                                    boolean riskCodedAsTwo;
                                    if(sumRisk && or <0){
                                        or = or*-1;
                                        riskCodedAsTwo = false;
                                        if (!(riskE.getAllele()==(var1.getRefAllele().getAlleleAsSnp()) || riskE.getAllele()==(var1.getRefAllele().getComplement().getAlleleAsSnp()))) {
                                            riskCodedAsTwo = true;
                                        }
                                    } else {
                                        riskCodedAsTwo = true;
                                        if (!(riskE.getAllele()==(var1.getRefAllele().getAlleleAsSnp()) || riskE.getAllele()==(var1.getRefAllele().getComplement().getAlleleAsSnp()))) {
                                            riskCodedAsTwo = false;
                                        }
                                    }

//                                    StringBuilder Genos = new StringBuilder();
//                                    StringBuilder Genos1 = new StringBuilder();
//                                    StringBuilder Genos2 = new StringBuilder();
//                                    StringBuilder Genos3 = new StringBuilder();
                                    for (int sample = 0; sample < var1.getSampleCalledDosages().length; sample++) {
                                        if (var1.getSampleCalledDosages()[sample] != -1) {
                                            if(riskCodedAsTwo){
                                                scores.getMatrix().setQuick(rowNr, sample, (scores.getMatrix().getQuick(rowNr, sample) + (or * var1.getSampleCalledDosages()[sample])));
//                                                Genos2.append(or * var1.getSampleCalledDosages()[sample]).append(",");
                                            } else {
                                                scores.getMatrix().setQuick(rowNr, sample, (scores.getMatrix().getQuick(rowNr, sample) + (or * Math.abs(var1.getSampleCalledDosages()[sample] - 2))));
//                                                Genos2.append(or * Math.abs(var1.getSampleCalledDosages()[sample]-2)).append(",");
                                            }
//                                            Genos.append(var1.getSampleCalledDosages()[sample]).append(",");
//                                            Genos1.append(var1.getSampleVariants().get(sample).toString());Genos1.append(",");
//                                            Genos3.append(scores.getMatrix().getQuick(rowNr, sample)).append(",");
                                        }
                                    }
                                    //I have now removed the conversion from 0 to 2 and vice versa from the calculations. Why is this needed?
//                                    System.out.println("");
//                                    System.out.println("SNP: "+riskE.getRsName());
//                                    System.out.println("Allele in data: "+var1.getRefAllele().toString());
//                                    System.out.println("Allele: "+riskE.getAllele());
//                                    System.out.println("Risk coded as two: "+riskCodedAsTwo);
//                                    System.out.println("Or: "+or);
//                                    System.out.println("Dosages: "+Genos.toString());
//                                    System.out.println("Genotypes: "+Genos1.toString());
//                                    System.out.println("Scores: "+Genos2.toString());
//                                    System.out.println("ScoresT: "+Genos3.toString());
//                                    System.out.println("");

                                    nrSNPs++;

                                    for (int t = snp + 1; t < nrSNPsThisChr; t++) {
                                        if (!excludeSNPs[t]) {
                                            RiskEntry riskE2 = valueE2.get(t);
                                            if (Math.abs(riskE2.getPos() - riskE.getPos()) <= windowSize[1]) {
                                                GeneticVariant var2 = genotypeData.getSnpVariantByPos(riskE2.getChr(), riskE2.getPos());
                                                try {
                                                    if (calculateRsquare(var1, var2) >= rSquare) {
                                                        excludeSNPs[t] = true;
                                                        excludeList.add(riskE2.getRsName());
                                                    }
                                                } catch (LdCalculatorException ex) {
                                                    Logger.getLogger(CalculateSimpleGeneticRiskScore.class.getName()).log(Level.SEVERE, null, ex);
                                                }
                                            }
                                        }
                                    }
                                }
                            }
                        }
                        if (debugMode) {
                            out.write("Total SNPs used: " + nrSNPs);
                            out.close();
                        }
                    } catch (IOException ex) {
                        Logger.getLogger(CalculateSimpleGeneticRiskScore.class.getName()).log(Level.SEVERE, null, ex);
                    }
                }
                p.iterate();
            }
        }
        p.close();
        return scores;
    }

}
