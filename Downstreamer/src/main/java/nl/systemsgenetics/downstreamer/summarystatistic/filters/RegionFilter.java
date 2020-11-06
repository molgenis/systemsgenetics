package nl.systemsgenetics.downstreamer.summarystatistic.filters;

import nl.systemsgenetics.downstreamer.summarystatistic.Locus;
import nl.systemsgenetics.downstreamer.summarystatistic.SummaryStatisticRecord;

/**
 * The type Region filter.
 */
public class RegionFilter implements SummaryStatisticsRecordFilter {

    private static FilterType filterType = FilterType.REGION_FILTER;
    private String sequenceName;
    private long upperBound;
    private long lowerBound;


    /**
     * Instantiates a new Region filter.
     *
     * @param locus the locus
     */
    public RegionFilter(Locus locus) {
        this.sequenceName = locus.getSequenceName();
        this.upperBound = locus.getEnd();
        this.lowerBound = locus.getStart();
    }

    /**
     * Instantiates a new Region filter.
     *
     * @param sequenceName the sequence name
     * @param upperBound   the upper bound
     * @param lowerBound   the lower bound
     */
    public RegionFilter(String sequenceName, long upperBound, long lowerBound) {
        this.sequenceName = sequenceName;
        this.upperBound = upperBound;
        this.lowerBound = lowerBound;
    }

    @Override
    public boolean passesFilter(SummaryStatisticRecord record) {
        if (record.getSequenceName().equals(sequenceName)) {
            return record.getPosition() >= lowerBound && record.getPosition() <= upperBound;
        } else {
            return false;
        }
    }

    @Override
    public FilterType getFilterType() {
        return filterType;
    }
}
