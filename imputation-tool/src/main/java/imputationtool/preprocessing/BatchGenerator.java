/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package imputationtool.preprocessing;

import java.io.IOException;
import java.util.ArrayList;
import java.util.HashMap;
import umcg.genetica.io.text.TextFile;
import umcg.genetica.io.trityper.TriTyperGenotypeData;

/**
 *
 * @author harmjan
 */
public class BatchGenerator {

    public void generateBatchesFromPedAndMap(String input, String output, int batchSize) {
    }

    public void generateBatchesFromTriTyper(String input, String output, int batchSize) throws IOException {
        TriTyperGenotypeData dataset = new TriTyperGenotypeData();
        dataset.load(input);

        ArrayList<String> individuals = new ArrayList<String>();
        HashMap<String, Boolean> individualsAffected = new HashMap<String, Boolean>();

        for (int i = 0; i < dataset.getIndividuals().length; i++) {
            if (dataset.getIsIncluded()[i]) {
                String ind = dataset.getIndividuals()[i];
                individuals.add(ind);
                individualsAffected.put(ind, dataset.getIsCase()[i]);
            }
        }
        generateBatches(individuals, individualsAffected, batchSize, output);
    }

    private void generateBatches(ArrayList<String> individuals, HashMap<String, Boolean> individualsAffected, int requestedBatchSize, String output) throws IOException {
        int numSamples = individuals.size();
        int numAffected = 0;

        ArrayList<String> affectedIndividuals = new ArrayList<String>();
        ArrayList<String> unaffectedIndividuals = new ArrayList<String>();

        for (int i = 0; i < numSamples; i++) {
            if (individualsAffected.get(individuals.get(i))) {
                numAffected++;
                affectedIndividuals.add(individuals.get(i));
            } else {
                unaffectedIndividuals.add(individuals.get(i));
            }
        }

        if (requestedBatchSize > numSamples) {
            requestedBatchSize = numSamples;
        }

        int numBatches = (int) Math.floor((double) numSamples / requestedBatchSize);

        int remainder = numSamples % requestedBatchSize;
        int[] finalBatchSizes = new int[numBatches];

        for (int i = 0; i < numBatches; i++) {
            finalBatchSizes[i] = requestedBatchSize;
        }

        int numAffectedPerBatch = 0;
        int numUnaffectedPerBatch = 0;

        if (numAffected > 0) {
            numAffectedPerBatch = (int) Math.floor((double) numAffected / numBatches);
            numUnaffectedPerBatch = (int) Math.floor((double) (numSamples - numAffected) / numBatches);
        }

        while (remainder > 0) {
            for (int b = 0; b < numBatches; b++) {
                if (remainder > 0) {
                    finalBatchSizes[b]++;
                    remainder--;
                }
            }
        }

        String[] batchNames = getBatches(numBatches);


        TextFile out = new TextFile(output + "/batches.txt", TextFile.W);
        int sum = 0;
        for (int b = 0; b < numBatches; b++) {
            // get batch descriptor
            String batchname = batchNames[b];

            int numSamplesInBatch = finalBatchSizes[b];
            int currentSamplesInBatch = 0;
            int currentSamplesUnaffecteInBatch = 0;
            int currentSamplesAffectedInBatch = 0;

            while (currentSamplesInBatch < numSamplesInBatch) {
                boolean chooseAffected = false;
                String ind = "";
                if (numAffected > 0) {
                    if (Math.random() > 0.5 && affectedIndividuals.size() > 0 && currentSamplesAffectedInBatch < numAffectedPerBatch) {
                        ind = affectedIndividuals.remove((int) Math.floor(Math.random() * affectedIndividuals.size()));
                        currentSamplesAffectedInBatch++;
                        currentSamplesInBatch++;

                    } else if (unaffectedIndividuals.size() > 0 && currentSamplesUnaffecteInBatch < numUnaffectedPerBatch) {
                        ind = unaffectedIndividuals.remove((int) Math.floor(Math.random() * unaffectedIndividuals.size()));
                        currentSamplesUnaffecteInBatch++;
                        currentSamplesInBatch++;
                    } else if (affectedIndividuals.size() > 0) {
                        ind = affectedIndividuals.remove((int) Math.floor(Math.random() * affectedIndividuals.size()));
                        currentSamplesAffectedInBatch++;
                        currentSamplesInBatch++;
                    } else if (unaffectedIndividuals.size() > 0) {
                        ind = unaffectedIndividuals.remove((int) Math.floor(Math.random() * unaffectedIndividuals.size()));
                        currentSamplesUnaffecteInBatch++;
                        currentSamplesInBatch++;
                    }

                } else {
                    ind = unaffectedIndividuals.remove((int) Math.floor(Math.random() * unaffectedIndividuals.size()));
                    currentSamplesInBatch++;
                }
                out.write(batchname + "\t" + ind + "\n");
            }
            sum += currentSamplesInBatch;
            System.out.println("Samples in batch\t" + batchname + " - " + numSamplesInBatch + "\tAffected:\t" + currentSamplesAffectedInBatch + "\tUnaffected\t" + (currentSamplesInBatch - currentSamplesAffectedInBatch));

        }
        System.out.println("Total samples: " + sum);
        out.close();

    }

    private String[] getBatches(int numBatches) {

        String[] batches = new String[numBatches];
        String[] alphabet = {"a", "b", "c", "d", "e", "f", "g", "h", "i", "j", "k", "l", "m", "n", "o", "p", "q", "r", "s", "t", "u", "v", "w", "x", "y", "z"};

        String firstletter = "a";
        int alphacounter = 0;
        int betacounter = 0;
        for (int i = 0; i < numBatches; i++) {
            if (i % 26 == 0) {
                firstletter = alphabet[alphacounter];
                alphacounter++;
                betacounter = 0;
            }

            batches[i] = firstletter + alphabet[betacounter];
            betacounter++;
        }
        return batches;
    }
}
