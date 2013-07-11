/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package umcg.genetica.io.trityper;

import java.io.IOException;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.Map.Entry;
import java.util.Set;
import java.util.HashSet;
import umcg.genetica.io.Gpio;
import umcg.genetica.io.text.TextFile;
import umcg.genetica.math.stats.Log2Transform;
import umcg.genetica.math.stats.QuantileNormalization;

/**
 *
 * @author harmjan
 */
public final class TriTyperGeneticalGenomicsDataset implements Comparable<TriTyperGeneticalGenomicsDataset> {

    private TriTyperGenotypeData genotypeData;
    private TriTyperExpressionData expressionData;
    private HashMap<String, String> genotypeToExpressionCouplings;
    private TriTyperGeneticalGenomicsDatasetSettings settings;
    private short[] expressionToGenotypeIdArray;
    private short totalGGSamples;
    private boolean expressionDataLoadedCorrectly = true;
    private short[] genotypeToExpressionIdArray;

    public TriTyperGeneticalGenomicsDataset(TriTyperGeneticalGenomicsDatasetSettings settings) throws IOException, Exception {
        this.settings = settings;

        settings.genotypeLocation = Gpio.formatAsDirectory(settings.genotypeLocation);

        if (settings.expressionLocation == null) {
            settings.expressionLocation = settings.genotypeLocation + "ExpressionData.txt";
        }

        // load the genotype metadata
        genotypeData = new TriTyperGenotypeData();
        genotypeData.load(settings.genotypeLocation);
        HashSet<String> includedExpressionIndividuals = new HashSet<String>();
        Boolean[] isIncluded = genotypeData.getIsIncluded();

        // preload the sample coupling file
        loadCouplings();

        // determine which expression samples to include
        Set<Entry<String, String>> entries = genotypeToExpressionCouplings.entrySet();
        for (Entry<String, String> entry : entries) {
            String genotypeIndividual = entry.getKey();
            Integer genotypeIndividualId = genotypeData.getIndividualId(genotypeIndividual);

            if (genotypeIndividualId != null && isIncluded[genotypeIndividualId]) {
                includedExpressionIndividuals.add(entry.getValue());
            }
        }

        // load the expression data
        expressionData = new TriTyperExpressionData();
        expressionData.confineToProbes(settings.tsProbesConfine);
        expressionData.setConfineToProbesThatMapToChromosome(settings.confineProbesToProbesMappingToAnyChromosome);
        expressionData.setConfineToProbesThatMapToChromosome(settings.confineProbesToProbesThatMapToChromosome);
        expressionData.setIncludeIndividuals(includedExpressionIndividuals);
        expressionDataLoadedCorrectly = expressionData.load(settings.expressionLocation, settings.probeannotation, settings.expressionplatform, (settings.cisAnalysis && settings.transAnalysis));

        pruneGenotypeToExpressionCouplings();


        if (settings.quantilenormalize) {
            QuantileNormalization.quantilenormalize(expressionData.getMatrix());
        }

        if (settings.logtransform) {
            Log2Transform.log2transform(expressionData.getMatrix());
        }


    }

    /**
     * @return the genotypeData
     */
    public TriTyperGenotypeData getGenotypeData() {
        return genotypeData;
    }

    /**
     * @param genotypeData the genotypeData to set
     */
    public void setGenotypeData(TriTyperGenotypeData genotypeData) {
        this.genotypeData = genotypeData;
    }

    /**
     * @return the expressionData
     */
    public TriTyperExpressionData getExpressionData() {
        return expressionData;
    }

    /**
     * @param expressionData the expressionData to set
     */
    public void setExpressionData(TriTyperExpressionData expressionData) {
        this.expressionData = expressionData;
    }

    /**
     * @return the genotypeToExpressionCouplings
     */
    public HashMap<String, String> getGenotypeToExpressionCouplings() {
        return genotypeToExpressionCouplings;
    }

    /**
     * @param genotypeToExpressionCouplings the genotypeToExpressionCouplings to
     * set
     */
    public void setGenotypeToExpressionCouplings(HashMap<String, String> genotypeToExpressionCouplings) {
        this.genotypeToExpressionCouplings = genotypeToExpressionCouplings;
    }

    /**
     * @return the settings
     */
    public TriTyperGeneticalGenomicsDatasetSettings getSettings() {
        return settings;
    }

    /**
     * @param settings the settings to set
     */
    public void setSettings(TriTyperGeneticalGenomicsDatasetSettings settings) {
        this.settings = settings;
    }

    public int getTotalGGSamples() {
        return totalGGSamples;
    }

    private void loadCouplings() throws IOException {
        genotypeToExpressionCouplings = new HashMap<String, String>();
        String genotypeToExpressionCoupling = settings.genotypeToExpressionCoupling;
        if (genotypeToExpressionCoupling != null && genotypeToExpressionCoupling.trim().length() > 0) {
            if (!Gpio.exists(genotypeToExpressionCoupling)) {
                throw new IOException("Error: genotype to expression coupling file: " + genotypeToExpressionCoupling + " does not exist.");
            }
            TextFile in = new TextFile(genotypeToExpressionCoupling, TextFile.R);
            String[] elems = in.readLineElemsReturnReference(TextFile.tab);

            while (elems != null) {
                if (elems.length > 1) {
                    String key = new String(elems[0].getBytes("UTF-8"));
                    String value = new String(elems[1].getBytes("UTF-8"));
                    if (genotypeToExpressionCouplings.get(key) != null) {
                        System.out.println("ERROR: your genotype to expression coupling file contains duplicate entries for individual: " + key);
                        System.exit(0);
                    } else {
                        genotypeToExpressionCouplings.put(key, value);
                    }

                }
                elems = in.readLineElemsReturnReference(TextFile.tab);
            }

            in.close();
        } else {
            Boolean[] isIncluded = genotypeData.getIsIncluded();
            int i = 0;
            String[] individuals = genotypeData.getIndividuals();
            for (String ind : individuals) {
                if (isIncluded[i] != null && isIncluded[i]) {
                    if (genotypeToExpressionCouplings.get(ind) != null) {
                        System.out.println("ERROR: your genotype data contains duplicate individuals: " + ind);
                        System.exit(0);
                    } else {
                        genotypeToExpressionCouplings.put(ind, ind);
                    }
                }
                i++;
            }
        }
    }

    public int[] getExpressionToGenotypeIdArray() {
        int[] intExpToGArr = new int[expressionToGenotypeIdArray.length];
        for(int i=0; i<intExpToGArr.length; i++){
            intExpToGArr[i] = expressionToGenotypeIdArray[i];
        }
        return intExpToGArr;
    }
    
    public short[] getExpressionToGenotypeIdArrayShort() {
        return expressionToGenotypeIdArray;
    }

    @Override
    public int compareTo(TriTyperGeneticalGenomicsDataset o) {
        int numIndsOther = o.getGenotypeData().getIndividuals().length;
        return genotypeData.getIndividuals().length - numIndsOther;
    }

    public boolean equals(TriTyperGeneticalGenomicsDataset o) {
        int numIndsOther = o.getGenotypeData().getIndividuals().length;
        if (genotypeData.getIndividuals().length == numIndsOther) {
            return true;
        } else {
            return false;
        }
    }

    /**
     * Permutes the mapping between genotype and gene expression samples
     */
    public void permuteSampleLables() {
        ArrayList<Short> alIndWGA = new ArrayList<Short>();
        int numSamples = expressionToGenotypeIdArray.length;
        for (int i = 0; i < numSamples; i++) {
            if (expressionToGenotypeIdArray[i] != -1) {
                alIndWGA.add(expressionToGenotypeIdArray[i]);
            }
        }

        short[] indWGANew = new short[numSamples];
        
        genotypeToExpressionIdArray = new short[genotypeData.getIndividuals().length];
        for (int i = 0; i < numSamples; i++) {
            if (expressionToGenotypeIdArray[i] == -1) {
                indWGANew[i] = -1;
            } else {
                short genotypeId = alIndWGA.remove((int) (Math.random() * (double) alIndWGA.size()));
                indWGANew[i] = genotypeId;
                
                genotypeToExpressionIdArray[genotypeId] = (short) i;
            }
        }
        expressionToGenotypeIdArray = indWGANew;
    }

    public void resetGenotypeToExpressionCouplings() throws IOException {
        loadCouplings();
    }

    public void pruneGenotypeToExpressionCouplings() {
        // now check whether each genotype is actually linked to an expression individual...
        String[] individuals = genotypeData.getIndividuals();

        Boolean[] isReallyIncluded = new Boolean[individuals.length];
        HashMap<String, String> realGenotypeToExpressionCouplings = new HashMap<String, String>();
        totalGGSamples = 0;
        for (int i = 0; i < isReallyIncluded.length; i++) {
            String genotypeInd = individuals[i];
            if (!genotypeToExpressionCouplings.containsKey(genotypeInd)) {
                isReallyIncluded[i] = false;
            } else {
                String coupledExpressionSample = genotypeToExpressionCouplings.get(genotypeInd);
                if (coupledExpressionSample != null) {
                    Integer expressionSampleId = expressionData.getIndividualId(coupledExpressionSample);
                    if (expressionSampleId == null) {
                        isReallyIncluded[i] = false;
                    } else {
                        isReallyIncluded[i] = true;
                        realGenotypeToExpressionCouplings.put(genotypeInd, coupledExpressionSample);
                        totalGGSamples++;
                    }
                }
            }
        }

        // exclude genotypes for which no expression data is available.
        genotypeData.setIsIncluded(isReallyIncluded);
        genotypeToExpressionCouplings = realGenotypeToExpressionCouplings;

        // couple expression IDs to genotype IDs for quick reference
        Set<Entry<String, String>> entries = realGenotypeToExpressionCouplings.entrySet();
        expressionToGenotypeIdArray = new short[totalGGSamples];
        HashSet<Integer> visitedNumbers = new HashSet<Integer>();
        for (Entry<String, String> entry : entries) {
            Integer expressionIndId = expressionData.getIndividualId(entry.getValue());
            Integer genotypeIndId = genotypeData.getIndividualId(entry.getKey());
            if (expressionIndId != null && genotypeIndId != null) {
                if (visitedNumbers.contains(expressionIndId)) {
                    System.out.println("ERROR: your dataset contains duplicate samples!");
                } else {
                    expressionToGenotypeIdArray[expressionIndId] = genotypeIndId.shortValue();
                    visitedNumbers.add(expressionIndId);
                }
            }
        }
    }

    public HashMap<Integer, Integer> getGenotypeToExpressionIdHash() {
        HashMap<Integer, Integer> gte = new HashMap<Integer, Integer>();
        int i = 0;
        for (int entry : expressionToGenotypeIdArray) {
            Integer genotypeIndId = entry;
            Integer expressionIndId = i;
            if (expressionIndId != null && genotypeIndId != null) {
                gte.put(genotypeIndId, expressionIndId);
            }
            i++;
        }
        return gte;
    }

    public HashMap<Integer, Integer> getExpressionToGenotypeIdHash() {
        HashMap<Integer, Integer> etg = new HashMap<Integer, Integer>();
        int i = 0;
        for (int entry : expressionToGenotypeIdArray) {
            Integer genotypeIndId = entry;
            Integer expressionIndId = i;
            if (expressionIndId != null && genotypeIndId != null) {
                etg.put(expressionIndId, genotypeIndId);
            }
            i++;
        }
        return etg;
    }

    /**
     * @return the expressionDataLoadedCorrectly
     */
    public boolean isExpressionDataLoadedCorrectly() {
        return expressionDataLoadedCorrectly;
    }

    /**
     * @param expressionDataLoadedCorrectly the expressionDataLoadedCorrectly to
     * set
     */
    public void setExpressionDataLoadedCorrectly(boolean expressionDataLoadedCorrectly) {
        this.expressionDataLoadedCorrectly = expressionDataLoadedCorrectly;
    }
}
