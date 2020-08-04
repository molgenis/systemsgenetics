package nl.systemsgenetics.depict2.summarystatistic.filters;

import nl.systemsgenetics.depict2.summarystatistic.SummaryStatisticRecord;

/**
 * The interface Summary statistics record filter.
 */
public interface SummaryStatisticsRecordFilter {

    /**
     * Passes filter boolean.
     *
     * @param record the record
     * @return the boolean
     */
    boolean passesFilter(SummaryStatisticRecord record);

    /**
     * Gets filter type.
     *
     * @return the filter type
     */
    FilterType getFilterType();
}
