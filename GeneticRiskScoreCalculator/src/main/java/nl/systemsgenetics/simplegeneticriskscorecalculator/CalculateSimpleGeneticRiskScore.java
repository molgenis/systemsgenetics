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
import java.util.Map;
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

    static DoubleMatrixDataset<String, String> calculate(RandomAccessGenotypeData genotypeData, HashMap<String, HashMap<String, ArrayList<RiskEntry>>> risks, File outputFolder, double rSquare, double windowSize) {
        DoubleMatrixDataset<String, String> scores = new DoubleMatrixDataset<String, String>(new ArrayList<String>(risks.keySet()), Arrays.asList(genotypeData.getSampleNames()));
        
        for(Map.Entry<String, HashMap<String, ArrayList<RiskEntry>>> riskScorePheno : risks.entrySet()){
            int rowNr = scores.getHashRows().get(riskScorePheno.getKey());
//            System.out.println(riskScorePheno.getKey() +"\t"+rowNr);
            
            try {
                TextFile out = new TextFile(outputFolder+File.separator+riskScorePheno.getKey()+".log", TextFile.W);
                
                out.write("SNPs used for GRS calculation:\n");
                double[] risksPerSample = new double[genotypeData.getSampleNames().length];
                int nrSNPs = 0;
                for(Map.Entry<String, ArrayList<RiskEntry>> e2 : riskScorePheno.getValue().entrySet()){

                    System.out.println("Processing chromosome:\t" + e2.getKey());
                    
                    int nrSNPsThisChr = e2.getValue().size();
                    boolean[] excludeSNPs = new boolean[nrSNPsThisChr];
                    for (int snp = 0; snp < nrSNPsThisChr; snp++) {
                        RiskEntry riskE = e2.getValue().get(snp);
                        if (!excludeSNPs[snp]) {
                            //System.out.println(snpID + "\t" + c + "\t" + chrPos + "\t" + object.doubleValue);
                            GeneticVariant var1 = genotypeData.getSnpVariantByPos(riskE.getChr(), riskE.getPos());
                            //Check if at least 75% of the sampels have information for the SNP otherwise it is removed by default.
                            if (var1.getCallRate() < 0.75) {
                                excludeSNPs[snp] = true;
                                continue;
                            }

                            out.write(riskE.InfoToString() + "\n");

                            double or = riskE.getOr();
                            int direction = 1;
                            if (riskE.getAllele().equals(var1.getMinorAllele().toString()) || riskE.getAllele().equals(var1.getMinorAllele().getComplement().toString())) {
                                direction = -1;
                            }

                            for (int sample = 0; sample < var1.getSampleCalledDosages().length; sample++) {
                                if (var1.getSampleCalledDosages()[sample] != -1) {
                                    scores.getMatrix().setQuick(rowNr, sample, (scores.getMatrix().getQuick(rowNr, sample)+(direction * or * var1.getSampleCalledDosages()[sample])));
                                }
                            }

                            nrSNPs++;

                            //Test here all snp entries around the test snp, selected same number of snps as window size in bases.
                            //So we ware almost sure that we test al relevant entries, not all positions have SNPs so this should be save. Actual distance is calculated in the if statement below.

                            for (int t = snp+1; t < nrSNPsThisChr; t++) {
                                RiskEntry riskE2 = e2.getValue().get(t);
                                if(Math.abs(riskE2.getPos()-riskE.getPos())<=windowSize){
                                    GeneticVariant var2 = genotypeData.getSnpVariantByPos(riskE2.getChr(), riskE2.getPos());
                                    if(var2.getCallRate()<0.75){
                                        excludeSNPs[t] = true;
                                        continue;
                                    }
                                    try {
                                        Ld ld = var1.calculateLd(var2);
                                        if(ld.getR2()>=rSquare){
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
                out.write("Total SNPs used: "+nrSNPs);
                out.close();
            } catch (IOException ex) {
                Logger.getLogger(CalculateSimpleGeneticRiskScore.class.getName()).log(Level.SEVERE, null, ex);
            }
        }
        
        
        
        return scores;
    }
    
}
