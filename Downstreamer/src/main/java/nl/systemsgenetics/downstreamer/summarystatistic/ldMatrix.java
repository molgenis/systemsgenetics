package nl.systemsgenetics.downstreamer.summarystatistic;

import org.apache.log4j.Logger;
import org.molgenis.genotype.RandomAccessGenotypeData;
import org.molgenis.genotype.util.LdCalculatorException;
import org.molgenis.genotype.variant.GeneticVariant;
import org.molgenis.genotype.variantFilter.VariantIdIncludeFilter;

import java.util.HashMap;
import java.util.Map;
import java.util.Set;

/**
 * The type Ld matrix.
 */
public class ldMatrix {

    private static Logger LOGGER = Logger.getLogger(ldMatrix.class);
    private double[][] r2matrix;
    private String[] dimnames;
    private Map<String, Integer> index;


    public ldMatrix(ldMatrix other) {
        this.r2matrix = other.getR2matrix();
        this.dimnames = other.getDimnames();
        this.index = other.getIndex();
    }


    /**
     * Instantiates a new Ld matrix. Pairwise R2 values are stored in the upper triangle of the matrix.
     * Current variant limit is 1000 to avoid accidental computation of crazy LD matrices.
     *
     * @param genotypes the genotypes
     * @param variants  the variants
     * @throws LdCalculatorException the ld calculator exception
     */
    public ldMatrix(RandomAccessGenotypeData genotypes, Set<String> variants) throws LdCalculatorException {
        this.calculateR2Matrix(genotypes, variants);
    }

    /**
     * Instantiates a new Ld matrix. Pairwise R2 values are stored in the upper triangle of the matrix.
     * Current variant limit is 1000 to avoid accidental computation of crazy LD matrices.
     *
     * @param genotypes the genotypes
     * @param locus    the locus
     * @throws LdCalculatorException the ld calculator exception
     */
    public ldMatrix(RandomAccessGenotypeData genotypes, Locus locus) throws LdCalculatorException {
        this.calculateR2Matrix(genotypes, locus.getRecords().keySet());
    }


    /**
     * Instantiates a new Ld matrix. Pairwise R2 values are stored in the upper triangle of the matrix.
     * Current variant limit is X to avoid accidental computation of crazy LD matrices.
     *
     * @param genotypes the genotypes
     * @param variants  the variants to calculate pairwise LD between
     * @throws LdCalculatorException the ld calculator exception
     */
    private void calculateR2Matrix(RandomAccessGenotypeData genotypes, Set<String> variants) throws LdCalculatorException {

        // TODO: either make max variant size variable, or down-sample the amount of variants if it exceeds the limit.
        // TODO: Refactor to DoubleMatrix2d or a single double array, should be easy since already using maps for indices.

        // Some sanity checking to avoid idiotic running times when many variants are provided
        if (variants.size() > 3000) {
            // throw new LdCalculatorException("More than 3000 variants provided. Currently this is not supported due to running time concerns");
            LOGGER.warn("Over 3000 variants provided, might take a while to run");
        }

        LOGGER.info(variants.size() + " variants provided for pairwise LD calculation");

        // Initiate the ld matrix and dimmames
        r2matrix = new double[variants.size()][variants.size()];
        dimnames =  variants.toArray(new String[variants.size()]);
        index = new HashMap<>();

        int nVar = 0;
        for (String variant: dimnames) {
            index.put(variant, nVar);
            nVar ++;
        }

        // This filter is to ensure no excess pairwise comparisons are done
        Map<String, GeneticVariant> geneticVariants = genotypes.getVariantIdMap(new VariantIdIncludeFilter(variants));

        // Counter to keep track of the duplicated calculations
        // To not calculate LD between the same variants set this to 1
        int k = 1;

        // Loop over all rows
        for (int i=0; i < dimnames.length; i++) {
            GeneticVariant var1 = geneticVariants.get(dimnames[i]);
            // Loop over all columns
            for (int j=0; j < dimnames.length; j++) {
                r2matrix[i][j] = -9;
                // calculate LD only if current column is higher or equal to k
                // This is done so the same value is not calculated twice
                if (j >= k ) {
                    GeneticVariant var2 = geneticVariants.get(dimnames[j]);
                    r2matrix[i][j] = var1.calculateLd(var2).getR2();
                }
            }
            k++;
        }
    }

    /**
     * Gets r 2 for a variant pair. r2 measures are stored in the upper triangle of the matrix. If the assesed
     * pair is in the lower triangle, the combination is flipped so the r2 is returned.
     *
     * @param variant1 the variant 1
     * @param variant2 the variant 2
     * @return the r 2
     */
    public double getR2(String variant1, String variant2) {
        int idx1 = index.get(variant1);
        int idx2 = index.get(variant2);
        double out = r2matrix[idx1][idx2];

        if (out == -9) {
            return r2matrix[idx2][idx1];
        } else {
            return out;
        }
    }


    public static Logger getLOGGER() {
        return LOGGER;
    }

    public double[][] getR2matrix() {
        return r2matrix;
    }

    public String[] getDimnames() {
        return dimnames;
    }

    public Map<String, Integer> getIndex() {
        return index;
    }
}
