package nl.systemsgenetics.downstreamer.summarystatistic.filters;


import nl.systemsgenetics.downstreamer.summarystatistic.SummaryStatisticRecord;

import java.util.Set;

/**
 * The type Variant id filter.
 */
public class VariantIdFilter implements SummaryStatisticsRecordFilter {

    private final FilterType filterType = FilterType.VARIANT_ID_FILTER;
    private Set<String> variantIds;

    /**
     * Instantiates a new Variant id filter.
     *
     * @param variantIds the variant ids
     */
    public VariantIdFilter(Set<String> variantIds) {
        this.variantIds = variantIds;
    }

    @Override
    public boolean passesFilter(SummaryStatisticRecord record) {
        return variantIds.contains(record.getPrimaryVariantId());
    }

    @Override
    public FilterType getFilterType() {
        return filterType;
    }
}
