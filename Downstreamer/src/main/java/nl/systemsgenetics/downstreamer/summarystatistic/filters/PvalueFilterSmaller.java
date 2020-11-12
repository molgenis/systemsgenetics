package nl.systemsgenetics.downstreamer.summarystatistic.filters;

import nl.systemsgenetics.downstreamer.summarystatistic.SummaryStatisticRecord;

/**
 * The type Pvalue filter smaller.
 */
public class PvalueFilterSmaller implements SummaryStatisticsRecordFilter {

    private static FilterType filterType = FilterType.PVALUE_FILTER_ST;
    private double pvalueThreshold;

    /**
     * Instantiates a new Pvalue filter smaller.
     *
     * @param pvalueThreshold the pvalue threshold
     */
    public PvalueFilterSmaller(double pvalueThreshold) {
        this.pvalueThreshold = pvalueThreshold;
    }

    @Override
    public boolean passesFilter(SummaryStatisticRecord record) {

        try {
            return record.getPvalue() <= pvalueThreshold;
        } catch (NullPointerException e) {
            throw new FieldNotAvailibleException("pvalue");
        }
    }

    @Override
    public FilterType getFilterType() {
        return filterType;
    }
}
