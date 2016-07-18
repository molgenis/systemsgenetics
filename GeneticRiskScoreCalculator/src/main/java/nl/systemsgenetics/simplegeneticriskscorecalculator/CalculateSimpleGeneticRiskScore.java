/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package nl.systemsgenetics.simplegeneticriskscorecalculator;

import java.io.File;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.HashMap;
import java.util.HashSet;
import java.util.LinkedHashMap;
import java.util.Map.Entry;
import java.util.logging.Level;
import java.util.logging.Logger;
import org.molgenis.genotype.RandomAccessGenotypeData;
import org.molgenis.genotype.util.Ld;
import org.molgenis.genotype.util.LdCalculatorException;
import org.molgenis.genotype.variant.GeneticVariant;
import umcg.genetica.io.text.TextFile;
import umcg.genetica.math.matrix2.DoubleMatrixDataset;

/**
 *
 * @author MarcJan
 */
class CalculateSimpleGeneticRiskScore {

    private static final String[] chrOrder = {"1", "2", "3", "4", "5", "6", "7", "8", "9", "10", "11", "12", "13", "14", "15", "16", "17", "18", "19", "20", "21", "22"};

    static DoubleMatrixDataset<String, String> calculate(RandomAccessGenotypeData genotypeData, HashMap<String, LinkedHashMap<String, HashMap<String, ArrayList<RiskEntry>>>> risks, File outputFolder, double rSquare, double windowSize) {
        ArrayList<String> keys = new ArrayList<String>();
        for (Entry<String, LinkedHashMap<String, HashMap<String, ArrayList<RiskEntry>>>> riskScorePheno : risks.entrySet()) {
            for (Entry<String, HashMap<String, ArrayList<RiskEntry>>> riskScorePheno2 : riskScorePheno.getValue().entrySet()) {
                keys.add(riskScorePheno.getKey() + riskScorePheno2.getKey());
            }
        }

        DoubleMatrixDataset<String, String> scores = new DoubleMatrixDataset<String, String>(keys, Arrays.asList(genotypeData.getSampleNames()));

        for (Entry<String, LinkedHashMap<String, HashMap<String, ArrayList<RiskEntry>>>> riskScorePheno : risks.entrySet()) {
            HashMap<String, HashSet<String>> chrExcludeList = new HashMap<String, HashSet<String>>();
            try {
                for (Entry<String, HashMap<String, ArrayList<RiskEntry>>> riskScorePheno2 : riskScorePheno.getValue().entrySet()) {
                    String NameOfEntry = riskScorePheno.getKey() + riskScorePheno2.getKey();
                    int rowNr = scores.getHashRows().get(NameOfEntry);
                    System.out.println(NameOfEntry);
                    TextFile out = new TextFile(outputFolder + File.separator + NameOfEntry + ".log", TextFile.W);

                    out.write("SNPs used for GRS calculation:\n");
                    int nrSNPs = 0;
                    for (int counter = 0; counter < chrOrder.length; counter++) {
//                        System.out.println("Processing chromosome:\t" + chrOrder[counter]);

                        if (riskScorePheno2.getValue().containsKey(chrOrder[counter])) {
                            if(!chrExcludeList.containsKey(chrOrder[counter])){
                                chrExcludeList.put(chrOrder[counter], new HashSet<String>());
                            }
                            HashSet<String> excludeList = chrExcludeList.get(chrOrder[counter]);
                            ArrayList<RiskEntry> valueE2 = riskScorePheno2.getValue().get(chrOrder[counter]);

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

                                    out.write(riskE.InfoToString() + "\n");

                                    double or = riskE.getOr();
                                    int direction = -1;
                                    if (!(riskE.getAllele().equals(var1.getRefAllele().toString()) || riskE.getAllele().equals(var1.getRefAllele().getComplement().toString()))) {
                                        direction = 1;
                                    }

                                    //        StringBuilder Genos = new StringBuilder();
                                    //        StringBuilder Genos2 = new StringBuilder();
                                    //        StringBuilder Genos3 = new StringBuilder();
                                    for (int sample = 0; sample < var1.getSampleCalledDosages().length; sample++) {
                                        if (var1.getSampleCalledDosages()[sample] != -1) {
                                            scores.getMatrix().setQuick(rowNr, sample, (scores.getMatrix().getQuick(rowNr, sample) + (direction * or * Math.abs(var1.getSampleCalledDosages()[sample] - 2))));
                                            //        Genos.append(Math.abs(var1.getSampleCalledDosages()[sample]-2)).append(",");
                                            //        Genos2.append(direction * or * Math.abs(var1.getSampleCalledDosages()[sample]-2)).append(",");
                                            //        Genos3.append(scores.getMatrix().getQuick(rowNr, sample)).append(",");
                                        }
                                    }
                                    //        System.out.println("");
                                    //        System.out.println("SNP: "+riskE.getRsName());
                                    //        System.out.println("Allele in data: "+var1.getRefAllele().toString());
                                    //        System.out.println("Allele: "+riskE.getAllele());
                                    //        System.out.println("Diction: "+direction);
                                    //        System.out.println("Or: "+or);
                                    //        System.out.println("Dosages: "+Genos.toString());
                                    //        System.out.println("Scores: "+Genos2.toString());
                                    //        System.out.println("ScoresT: "+Genos3.toString());
                                    //        System.out.println("");

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
                                                    Ld ld = var1.calculateLd(var2);
                                                    if (ld.getR2() >= rSquare) {
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

                    }
                    out.write("Total SNPs used: " + nrSNPs);
                    out.close();
                }

            } catch (IOException ex) {
                Logger.getLogger(CalculateSimpleGeneticRiskScore.class.getName()).log(Level.SEVERE, null, ex);
            }
        }

        return scores;
    }

    static DoubleMatrixDataset<String, String> calculateTwoStages(RandomAccessGenotypeData genotypeData, HashMap<String, LinkedHashMap<String, HashMap<String, ArrayList<RiskEntry>>>> risks, File outputFolder, double rSquare, double[] windowSize) {
        ArrayList<String> keys = new ArrayList<String>();
        for (Entry<String, LinkedHashMap<String, HashMap<String, ArrayList<RiskEntry>>>> riskScorePheno : risks.entrySet()) {
            for (Entry<String, HashMap<String, ArrayList<RiskEntry>>> riskScorePheno2 : riskScorePheno.getValue().entrySet()) {
                keys.add(riskScorePheno.getKey() + riskScorePheno2.getKey());
            }
        }

        DoubleMatrixDataset<String, String> scores = new DoubleMatrixDataset<String, String>(keys, Arrays.asList(genotypeData.getSampleNames()));

        for (Entry<String, LinkedHashMap<String, HashMap<String, ArrayList<RiskEntry>>>> riskScorePheno : risks.entrySet()) {
            HashMap<String, HashSet<String>> chrExcludeList = new HashMap<String, HashSet<String>>();
            for (Entry<String, HashMap<String, ArrayList<RiskEntry>>> riskScorePheno2 : riskScorePheno.getValue().entrySet()) {
                try {
                    String NameOfEntry = riskScorePheno.getKey() + riskScorePheno2.getKey();
                    int rowNr = scores.getHashRows().get(NameOfEntry);
                    System.out.println(NameOfEntry);
                    TextFile out = new TextFile(outputFolder + File.separator + NameOfEntry + ".log", TextFile.W);

                    out.write("SNPs used for GRS calculation:\n");
                    int nrSNPs = 0;
                    for (int counter = 0; counter < chrOrder.length; counter++) {
//                        System.out.println("Processing chromosome:\t" + chrOrder[counter]);
                        if (riskScorePheno2.getValue().containsKey(chrOrder[counter])) {
                            if(!chrExcludeList.containsKey(chrOrder[counter])){
                                chrExcludeList.put(chrOrder[counter], new HashSet<String>());
                            }

                            HashSet<String> excludeList = chrExcludeList.get(chrOrder[counter]);
                            ArrayList<RiskEntry> valueE2 = riskScorePheno2.getValue().get(chrOrder[counter]);

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
                                                    Ld ld = var1.calculateLd(var2);
                                                    if (ld.getR2() >= rSquare) {
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
                                    
                                    out.write(riskE.InfoToString() + "\n");

                                    double or = riskE.getOr();
                                    int direction = -1;
                                    if (!(riskE.getAllele().equals(var1.getRefAllele().toString()) || riskE.getAllele().equals(var1.getRefAllele().getComplement().toString()))) {
                                        direction = 1;
                                    }

                                    //        StringBuilder Genos = new StringBuilder();
                                    //        StringBuilder Genos2 = new StringBuilder();
                                    //        StringBuilder Genos3 = new StringBuilder();
                                    for (int sample = 0; sample < var1.getSampleCalledDosages().length; sample++) {
                                        if (var1.getSampleCalledDosages()[sample] != -1) {
                                            scores.getMatrix().setQuick(rowNr, sample, (scores.getMatrix().getQuick(rowNr, sample) + (direction * or * Math.abs(var1.getSampleCalledDosages()[sample] - 2))));
                                            //        Genos.append(Math.abs(var1.getSampleCalledDosages()[sample]-2)).append(",");
                                            //        Genos2.append(direction * or * Math.abs(var1.getSampleCalledDosages()[sample]-2)).append(",");
                                            //        Genos3.append(scores.getMatrix().getQuick(rowNr, sample)).append(",");
                                        }
                                    }
                                    //        System.out.println("");
                                    //        System.out.println("SNP: "+riskE.getRsName());
                                    //        System.out.println("Allele in data: "+var1.getRefAllele().toString());
                                    //        System.out.println("Allele: "+riskE.getAllele());
                                    //        System.out.println("Diction: "+direction);
                                    //        System.out.println("Or: "+or);
                                    //        System.out.println("Dosages: "+Genos.toString());
                                    //        System.out.println("Scores: "+Genos2.toString());
                                    //        System.out.println("ScoresT: "+Genos3.toString());
                                    //        System.out.println("");

                                    nrSNPs++;

                                    for (int t = snp + 1; t < nrSNPsThisChr; t++) {
                                        if (!excludeSNPs[t]) {
                                            RiskEntry riskE2 = valueE2.get(t);
                                            if (Math.abs(riskE2.getPos() - riskE.getPos()) <= windowSize[1]) {
                                                GeneticVariant var2 = genotypeData.getSnpVariantByPos(riskE2.getChr(), riskE2.getPos());
                                                try {
                                                    Ld ld = var1.calculateLd(var2);
                                                    if (ld.getR2() >= rSquare) {
                                                        excludeSNPs[t] = true;
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

                    }
                    out.write("Total SNPs used: " + nrSNPs);
                    out.close();
                } catch (IOException ex) {
                    Logger.getLogger(CalculateSimpleGeneticRiskScore.class.getName()).log(Level.SEVERE, null, ex);
                }

            }
        }

        return scores;
    }

}
