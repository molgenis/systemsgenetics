package nl.systemsgenetics.downstreamer.summarystatistic.filters;

import java.io.File;
import java.io.IOException;
import java.util.Arrays;
import java.util.HashSet;
import java.util.List;
import java.util.Set;

/**
 * The enum Filter type.
 */
public enum FilterType {

    /**
     * Variant id filter filter type.
     */
    VARIANT_ID_FILTER,
    /**
     * Chromosome filter filter type.
     */
    CHROMOSOME_FILTER,
    /**
     * Transition filter filter type.
     */
    TRANSITION_FILTER,
    /**
     * Null record filter filter type.
     */
    NULL_RECORD_FILTER,
    /**
     * Pvalue filter st filter type.
     */
    PVALUE_FILTER_ST,
    /**
     * Pvalue filter gt filter type.
     */
    PVALUE_FILTER_GT,
    /**
     * Region filter filter type.
     */
    REGION_FILTER;

    /**
     * Create filter summary statistics record filter.
     *
     * @param arg the arg
     * @return the summary statistics record filter
     */
    public static SummaryStatisticsRecordFilter createFilter(String arg) {
       /* String filterType = arg.trim().split(";")[0];
        String args = "";

        try {
            switch (FilterType.valueOf(filterType)) {
                case VARIANT_ID_FILTER:
                    args = arg.trim().split(";")[1];
                    Set<String> variantIds;
                    if (new File(args).exists()) {
                        variantIds = GenericSummaryStatisticFileReader.readFileAsSet(args);
                    } else {
                        variantIds = new HashSet<>(Arrays.asList(args.split(",")));
                    }
                    return new VariantIdFilter(variantIds);
                case CHROMOSOME_FILTER:
                    args = arg.trim().split(";")[1];
                    List<String> sequenceNames = Arrays.asList(args.split(","));
                    return new SequenceNameFilter(sequenceNames);
                case TRANSITION_FILTER:
                    return new TransitionSiteFilter();
                case NULL_RECORD_FILTER:
                    return new NullRecordFilter();
                case PVALUE_FILTER_ST:
                    args = arg.trim().split(";")[1];
                    return new PvalueFilterSmaller(Double.parseDouble(args));
                case PVALUE_FILTER_GT:
                    args = arg.trim().split(";")[1];
                    return new PvalueFilterLarger(Double.parseDouble(args));
                case REGION_FILTER:
                    args = arg.trim().split(";")[1];
                    String sequenceName = args.split(":")[0];
                    String locus = args.split(":")[1];
                    long lower = Long.parseLong(locus.split("-")[0]);
                    long upper = Long.parseLong(locus.split("-")[1]);
                    return new RegionFilter(sequenceName, upper, lower);
            }

        } catch (IllegalArgumentException e) {
            throw new IllegalArgumentException("Invalid filter type specified: " + filterType);
        } catch (IOException e) {
            throw new IllegalArgumentException("Invalid filter file specified: " + args);
        }*/

        return null;
    }

}