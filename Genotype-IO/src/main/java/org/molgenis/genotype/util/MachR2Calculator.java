package org.molgenis.genotype.util;

import org.molgenis.genotype.GenotypeDataException;

/**
 *
 * @author Patrick Deelen
 */
public class MachR2Calculator {

    /**
     * Calculate the MACH r2 measure
     *
     * For formula see: doi:10.1038/nrg2796 S3
     *
     * @param dosages
     * @return
     */
    public static double calculateMachR2(float[][] probs) {

        //Here we perform the conversion to dosages manually to make sure all probs are used regardless of the calling threshold
        float[] dosages = ProbabilitiesConvertor.convertProbabilitiesToDosage(probs, 0);

        int nonMissingCount = 0;
        double dosageSum = 0;
        double dosageSqrSum = 0;

        for (float dosage : dosages) {
            if (dosage > 2) {
                throw new GenotypeDataException("Error in calculating MACH r2, found dosage larger than 2: " + dosage);
            }
            if (dosage >= 0) {
                ++nonMissingCount;
                dosageSum += dosage;
                dosageSqrSum += Math.pow(dosage, 2);
            }
            //else missing and ignore
        }

        double estimatedAlleleFrequency = dosageSum / (2 * nonMissingCount);

        if (estimatedAlleleFrequency <= 0 || estimatedAlleleFrequency >= 1) {
            return 1;
        }
        
        double machR2 = ((dosageSqrSum / nonMissingCount) - Math.pow((dosageSum / nonMissingCount), 2)) / (2 * estimatedAlleleFrequency * (1 - estimatedAlleleFrequency));
        return machR2 > 1.0 ? 1.0 : machR2;
    }

}
