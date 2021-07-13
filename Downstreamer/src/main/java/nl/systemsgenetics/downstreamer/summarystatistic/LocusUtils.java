package nl.systemsgenetics.downstreamer.summarystatistic;


import nl.systemsgenetics.downstreamer.summarystatistic.filters.PvalueFilterSmaller;
import nl.systemsgenetics.downstreamer.summarystatistic.filters.RegionFilter;
import nl.systemsgenetics.downstreamer.summarystatistic.filters.SummaryStatisticsRecordFilter;
import nl.systemsgenetics.downstreamer.summarystatistic.filters.VariantIdFilter;
import org.apache.commons.cli.MissingOptionException;
import org.apache.log4j.Logger;
import org.molgenis.genotype.RandomAccessGenotypeData;
import org.molgenis.genotype.util.Ld;
import org.molgenis.genotype.util.LdCalculatorException;
import org.molgenis.genotype.variant.GeneticVariant;
import org.molgenis.genotype.variantFilter.VariantIdIncludeFilter;

import java.util.*;

/**
 * The type Region utils.
 */
public class LocusUtils {

    private static final Logger LOGGER = Logger.getLogger(LocusUtils.class);


    /**
     * Select top snps based on a given window size and p-value threshold
     *
     * @param args    the args
     * @param records the records
     * @return the map
     * @throws MissingOptionException the missing option exception
     */
    public static Map<String, SummaryStatisticRecord> selectTopHits(String[] args, Map<String, SummaryStatisticRecord> records) throws MissingOptionException {

        LOGGER.info("Selecting the top effects in a certain locus");

        int nTotal = records.size();
        double upper = Double.parseDouble(args[0]);
        int window = Integer.parseInt(args[1]);

        LOGGER.info("P threshold: " + upper);
        LOGGER.info("Window size: " + window);

        Map<String, SummaryStatisticRecord> topHits = new HashMap<>();
        // Define the regions
        List<Locus> loci = makeLoci(records, new PvalueFilterSmaller(upper), new PvalueFilterSmaller(2), window);

        for (Locus locus : loci) {
            SummaryStatisticRecord topHit = locus.getMinPvalRecord();
            topHits.put(topHit.getPrimaryVariantId(), topHit);
        }

        LOGGER.info(topHits.size() + " top effects out of " + nTotal);
        return topHits;
    }

    /**
     * Clump locus objects that have a ld matrix available. Clumps all variants in the locus regardless of p-value.
     * If you want to set a p-value threshold please create the locus objects trough "makeRegions"
     *
     * @param ldThresh  the ld threshold
     * @param inputLoci the top regions
     * @return the list
     * @throws LdCalculatorException the ld calculator exception if LD matrix is not available
     */
    public static List<Locus> selectIndepTopHits(double ldThresh, List<Locus> inputLoci, Set<String> availibleVariants) throws LdCalculatorException {

        List<Locus> inputRegions2 = new ArrayList<>(inputLoci);
        List<Locus> topLoci = new ArrayList<>();

        // Perform clumping in the previously defined regions
        for (Locus topLocus : inputRegions2) {

            // Select all the available variants for that sample
            List<SummaryStatisticsRecordFilter> filters = new ArrayList<>();
            filters.add(new VariantIdFilter(availibleVariants));
            topLocus.filterRecords(filters);

            // If there are no variants in the locus remaining, exit
            if (topLocus.getRecords().size() == 0) {
                continue;
            }

            // Initialize the top variant as the first independent one
            SummaryStatisticRecord indep = topLocus.getMinPvalRecord();
            topLocus.addIndepVariant(indep.getPrimaryVariantId());

            // If there is only one top hit in the locus skip the LD calculation
            if (topLocus.getRecords().size() == 1) {
                topLoci.add(topLocus);
                continue;
            }

            if (topLocus.getLdMatrix() == null) {
                throw new LdCalculatorException("LD matrix is missing from locus object. Cannot perform clumping. Please instanciate locus object with LD matrix");
            }

            // Loop over all variants in the locus passing the lower bound
            for (SummaryStatisticRecord record : topLocus.getRecords().values()) {

                // Check the LD with each defined independent variant in the locus
                boolean isIndependentVariant = true;

                for (String indepVar : topLocus.getIndepVariants()) {

                    // Check the R2 between the 2 variants
                    if (ldThresh < topLocus.getLdMatrix().getR2(indepVar, record.getPrimaryVariantId())) {
                        isIndependentVariant = false;
                        break;
                    }
                }

                // If the LD is lower then the threshold add a new independent effect to the locus
                if (isIndependentVariant) {
                    topLocus.addIndepVariant(record.getPrimaryVariantId());
                }
            }

            topLoci.add(topLocus);
        }

        return topLoci;
    }

    /**
     * Select top effects independent of each other based on LD (clumping)
     *
     * @param args      the args
     * @param reference the reference
     * @param records   the records
     * @return the list
     * @throws LdCalculatorException the ld calculator exception
     */
    public static List<Locus> selectIndepTopHits(String[] args, RandomAccessGenotypeData reference, Map<String, SummaryStatisticRecord> records) throws LdCalculatorException {

        // Define parameters
        // TODO: Fix this atm rater ugly to parse these here
        double upper = Double.parseDouble(args[0]);
        double lower = Double.parseDouble(args[1]);
        int window = Integer.parseInt(args[2]);
        double ldThresh = Double.parseDouble(args[3]);

        //LOGGER.info("Clumping using following parameters:");
        //LOGGER.info("\tUpper P threshold:\t" + upper);
        //LOGGER.info("\tLower P threshold:\t" + lower);
        //LOGGER.info("\tWindow size:\t\t" + window);
        //LOGGER.info("\tLD R2 threshold:\t" + ldThresh);
        //LOGGER.warn("Current implementation does not support a missing top variant, if top variant is missing, the locus is omitted");

        // Define the regions to clump in based on upper and lower bounds
        List<Locus> topLoci = makeLoci(records, new PvalueFilterSmaller(upper), new PvalueFilterSmaller(lower), window);

        // Perform clumping in the defined regions
        for (Locus topLocus : topLoci) {

            // Initialize the top variant as the first independent one
            SummaryStatisticRecord indep = topLocus.getMinPvalRecord();
            topLocus.addIndepVariant(indep.getPrimaryVariantId());

            // If there is only one top hit in the locus skip the LD calculation
            if (topLocus.getRecords().size() == 1) {
                continue;
            }

            // Select the variant in the reference data
            GeneticVariant iVar = reference.getSnpVariantByPos(indep.getSequenceName().replace("chr", ""), (int) indep.getPosition());
            GeneticVariant oVar;

            // Skip if the reference is missing
            if (iVar == null) {
                LOGGER.warn("Top variant missing from the reference, skipping locus");
                continue;
            }

            // Loop over all variants in the locus passing the lower bound
            for (SummaryStatisticRecord record : topLocus.getRecords().values()) {

                // Check the LD with each defined independent variant in the locus
                boolean isIndependentVariant = true;

                for (String indepVar : topLocus.getIndepVariants()) {
                    indep = topLocus.getRecords().get(indepVar);
                    iVar = reference.getSnpVariantByPos(indep.getSequenceName().replace("chr", ""), (int) indep.getPosition());
                    oVar = reference.getSnpVariantByPos(record.getSequenceName().replace("chr", ""), (int) record.getPosition());

                    // If either the independent variant or the other is missing, skip
                    if (iVar == null | oVar == null) {
                        continue;
                    }

                    // Check the R2 between the 2 variants
                    Ld curLd = iVar.calculateLd(oVar);
                    if (ldThresh < curLd.getR2()) {
                        isIndependentVariant = false;
                        break;
                    }
                }

                // If the LD is lower then the threshold add a new independent effect to the locus
                if (isIndependentVariant) {
                    topLocus.addIndepVariant(record.getPrimaryVariantId());
                }
            }
        }

        return topLoci;
    }

    /**
     * Select all the variants around a certain snp based on a window
     *
     * @param args    the args
     * @param topHits the top hits
     * @param records the records
     * @return the map
     */
    public static Map<String, Locus> selectTopLoci(String[] args, Map<String, SummaryStatisticRecord> topHits, Map<String, SummaryStatisticRecord> records) {

        LOGGER.info("Selecting the loci around the top effects");

        int window = Integer.parseInt(args[1]);
        Map<String, Locus> out = new HashMap<>();
        SummaryStatisticsRecordFilter filter;

        for (SummaryStatisticRecord topHit : topHits.values()) {

            // Determine the regions to extract
            Locus curLocus = new Locus(topHit, window);
            filter = new RegionFilter(curLocus);

            // Apply the specified filter to the set
            for (SummaryStatisticRecord curRecord : records.values()) {
                if (filter.passesFilter(curRecord)) {
                    curLocus.addRecord(curRecord);
                }
            }
            out.put(topHit.getPrimaryVariantId(), curLocus);

            LOGGER.info("Selected locus around: " + topHit.getPrimaryVariantId());
        }

        return out;
    }

    /**
     * Generate regions based on a window and the p-value. Variants are first sorted on p-value. In case of ties beta's are used if present.
     * Then the ranked list is iteratively looped over to identify any variants within a window of the current top SNPs.
     * If they are within a window they are added to the locus. If they are not within an existing locus a new locus is made
     * with that SNP as the top effect if it meets the lowerFiler. The resulting locus object can then be used for clumping
     * based on LD
     * <p>
     * Algo is greedy, adds overlapping snps between regions only to the first locus it encounters, this is always the
     * locus with the most significant top snp
     *
     * @param records     the records
     * @param upperFilter the upper filter
     * @param lowerFilter the lower filter
     * @param window      the window
     * @return the list
     */
    public static List<Locus> makeLoci(Map<String, SummaryStatisticRecord> records, PvalueFilterSmaller upperFilter, PvalueFilterSmaller lowerFilter, int window) {
        List<Locus> loci = new ArrayList<>();

        // first sort the statistics by pvalue
        List<SummaryStatisticRecord> statisticsList = new ArrayList<>(records.values());
        Collections.sort(statisticsList);
        //LOGGER.warn("Sorting pvalues. Does not properly support ties & 0 p-values");

        for (SummaryStatisticRecord curRecord : statisticsList) {

            // If the record meets at least the lower p-value threshold check if it is the locus
            if (lowerFilter.passesFilter(curRecord)) {
                // Check if it is located in an existing locus
                boolean isInExistingRegion = false;
                for (Locus locus : loci) {
                    if (locus.isInRegion(curRecord)) {
                        locus.addRecord(curRecord);
                        isInExistingRegion = true;

                        // Only add a snp to the first locus in case of overlaps
                        break;
                    }
                }
                // If it is not in an existing locus and it meets the upper threshold for significance create a new locus
                if (!isInExistingRegion) {
                    if (upperFilter.passesFilter(curRecord)) {
                        loci.add(new Locus(curRecord, window));
                    }
                }
            }
        }

        return loci;
    }

    /**
     * Generate a list of loci with pre-defined index snps. No checking for overlaps or Pvalue is done, so its assumed
     * the index snps are correct.
     * @param records
     * @param indexSnps
     * @param window
     * @return
     */

    public static List<Locus> makeLociWithGivenIndexSnps(Map<String, SummaryStatisticRecord> records, Set<String> indexSnps, int window) {
        List<Locus> loci = new ArrayList<>();

        for (String curSnp : indexSnps) {
            SummaryStatisticRecord  curRecord = records.get(curSnp);
            loci.add(new Locus(curRecord, window));
        }

        return loci;
    }

    /**
     * Calculate LD r2 measures for all variants in a given locus of the specified variant record.
     *
     * @param genotypeData the genotype data
     * @param record       the record
     * @param window       the window
     * @param ldThreshold  the ld threshold
     * @return the map
     * @throws LdCalculatorException the ld calculator exception
     */
    public static Map<String, Double> selectTaggingVariants(RandomAccessGenotypeData genotypeData, SummaryStatisticRecord record, int window, double ldThreshold) throws LdCalculatorException {

        Iterable<GeneticVariant> locus = genotypeData.getVariantsByRange(record.getSequenceName().replace("chr", ""), (int) (record.getPosition() - window), (int) (record.getPosition() + window));
        GeneticVariant topHit = genotypeData.getVariantIdMap(new VariantIdIncludeFilter(record.getPrimaryVariantId())).get(record.getPrimaryVariantId());

        if (topHit == null) {
            throw new IllegalArgumentException("Variant not present in reference genotype, cannot perform tagging");
        }

        Map<String, Double> out = new HashMap<>();

        for (GeneticVariant refVar : locus) {

            Ld ld = topHit.calculateLd(refVar);

            if (ld.getR2() >= ldThreshold) {
                out.put(refVar.getPrimaryVariantId(), ld.getR2());
            }
        }

        return out;
    }

    /**
     * Checks if A overlaps B.
     *
     * @param a range a
     * @param b range b
     * @return
     */
    public static boolean partialGenomicRangeOverlap(OverlappableGenomicRange a, OverlappableGenomicRange b) {
        return partialGenomicRangeOverlapWindow(a, b, 0);
    }


    /**
     * Checks if A overlaps B.
     *
     * @param a The smaller range
     * @param b The larger range
     * @return
     */
    public static boolean partialGenomicRangeOverlapWindow(OverlappableGenomicRange a, OverlappableGenomicRange b, int window) {

        if (a.getSequenceName().toLowerCase().replace("chr", "").equals(b.getSequenceName().toLowerCase().replace("chr", ""))) {

            final int windowStart = Math.min(a.getStart(), a.getEnd()) - window;
            final int windowEnd = Math.max(a.getStart(), a.getEnd()) + window;

            if (b.getStart() >= windowStart && b.getStart() <= windowEnd) {
                return true;
            } else if (b.getEnd() >= windowStart && b.getEnd() <= windowEnd) {
                return true;
            }

            //if a is fully within b we need another test
            final int bStart = Math.min(b.getStart(), b.getEnd());
            final int bEnd = Math.max(b.getStart(), b.getEnd());

            if (windowStart >= bStart && windowStart <= bEnd) {
                return true;
            }


        }

        return false;


    }

}
