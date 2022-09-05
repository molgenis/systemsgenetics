package nl.systemsgenetics.downstreamer.summarystatistic;

import htsjdk.samtools.util.Interval;
import nl.systemsgenetics.downstreamer.gene.Gene;
import nl.systemsgenetics.downstreamer.summarystatistic.filters.SummaryStatisticsRecordFilter;

import java.util.*;
import java.util.function.Function;
import java.util.stream.Collectors;

/**
 * The type Region.
 */
public class Locus extends Interval implements OverlappableGenomicRange {

    private Map<String, SummaryStatisticRecord> records;
    private Set<String> indepVariants;
    private SummaryStatisticRecord cachedTopHit;
    private ldMatrix ldMatrix;
    private List<Gene> genes;

    /**
     * Instantiates a new Region.
     */
    public Locus() {
        super(null, 0, 0);
    }

    /**
     * Instantiates a new Locus.
     *
     * @param record the record
     * @param window the window
     */
    public Locus(SummaryStatisticRecord record, int window) {

        //this.start = record.getPosition() - window;
        //this.end = record.getPosition() + window;
        //this.sequenceName = record.getSequenceName();
        super(record.getSequenceName(), record.getPosition() - window, record.getPosition() + window);
        this.records = new HashMap<>();
        this.records.put(record.getPrimaryVariantId(), record);
        this.indepVariants = new HashSet<>();
        this.cachedTopHit = null;
        this.genes = new ArrayList<>();
    }

    /**
     * Instantiates a new Locus.
     *
     * @param other the locus to clone
     */
    public Locus(Locus other) {
        //this.start = other.getStart();
        //this.end = other.getEnd();
        // this.sequenceName = other.getSequenceName();
        super(other);
        this.records = new HashMap<>();

        for (String record : other.getRecords().keySet()) {
            SummaryStatisticRecord old = other.getRecords().get(record);
            records.put(record, new SummaryStatisticRecord(old));
        }

        this.indepVariants = new HashSet<>();

        for (String variantId : other.getIndepVariants()) {
            indepVariants.add(new String(variantId));
        }

        this.cachedTopHit = null;
        this.ldMatrix = new ldMatrix(other.getLdMatrix());
    }

    public ldMatrix getLdMatrix() {
        return ldMatrix;
    }

    public void setLdMatrix(ldMatrix ldMatrix) {
        this.ldMatrix = ldMatrix;
    }

    /**
     * Add record.
     *
     * @param record the record
     */
    public void addRecord(SummaryStatisticRecord record) {

        if (isInRegion(record)) {
            records.put(record.getPrimaryVariantId(), record);
        } else {
            throw new IllegalArgumentException("Records position not in locus");
        }

    }

    /**
     * Add an (LD) independent variant. If the variant is not in the record map
     * it will NOT be added
     *
     * @param variantId the variant id
     */
    public void addIndepVariant(String variantId) {
        if (records.keySet().contains(variantId)) {
            this.indepVariants.add(variantId);
        }
    }

    /**
     * Add multiple (LD) independent variants. If the variants are not present the record map
     * they will NOT be added
     *
     * @param variantIds the variant ids
     */
    public void addIndepVariants(Collection<String> variantIds) {
        for (String variantId : variantIds) {
            addIndepVariant(variantId);
        }
    }

    /**
     * Gets the cached (LD) independent variants variants.
     *
     * @return the independent variants
     */
    public Set<String> getIndepVariants() {
        return indepVariants;
    }

    /**
     * Add a gene to the locus.
     *
     * @param gene
     */
    public void addGene(Gene gene) {
        genes.add(gene);
    }


    /**
     * Add all genes to the locus.
     *
     * @param genesToAdd
     */
    public void addGenes(Collection<Gene> genesToAdd) {
        genes.addAll(genesToAdd);
    }

    /**
     * Gets the (LD) independent summary statistic records.
     *
     * @return the independent summary statistic records
     */
    public Collection<SummaryStatisticRecord> getIndepSummaryStatisticRecords() {
        return indepVariants.stream()
                .filter(records::containsKey)
                .collect(Collectors.toMap(Function.identity(), records::get))
                .values();
    }

    /**
     * Gets record.
     *
     * @param id the id
     * @return the record
     */
    public SummaryStatisticRecord getRecord(String id) {
        return records.get(id);
    }

    /**
     * Is in locus boolean.
     *
     * @param record the record
     * @return the boolean
     */
    public boolean isInRegion(SummaryStatisticRecord record) {

        if (record.getSequenceName().equals(getContig())) {
            return record.getPosition() >= getStart() && record.getPosition() <= getEnd();
        } else {
            return false;
        }

    }

    /**
     * Gets the record with the min p-value or the current cached top hit.
     *
     * @return the min pval
     */
    public SummaryStatisticRecord getMinPvalRecord() {

        //TODO: look at the caching later
        //if (cachedTopHit == null){
        deterimineCachedTopHit();
        //}

        return cachedTopHit;
    }


    /**
     * Get the minimal pvalue in the locus.
     *
     * @return the minimal pvalue in the locus
     */
    public double getMinPvalInLocus() {
        return this.getMinPvalRecord().getPvalue();
    }

    /**
     * Determine the cached top hit.
     */
    public void deterimineCachedTopHit() {
        double minPval = 1;
        SummaryStatisticRecord finalRecord = null;

        for (SummaryStatisticRecord record : records.values()) {
            if (record.getPvalue() < minPval) {
                minPval = record.getPvalue();
                finalRecord = record;
            }
        }
        cachedTopHit = finalRecord;
    }


    /**
     * Gets records.
     *
     * @return the records
     */
    public Map<String, SummaryStatisticRecord> getRecords() {
        return records;
    }


    public void filterRecords(List<SummaryStatisticsRecordFilter> filters) {
        SummaryStatisticRecordUtils.filterSummaryStatistics(this.records, filters);
    }

    @Override
    public String getSequenceName() {
        return getContig();
    }

    @Override
    @Deprecated
    public boolean isOverlapping(OverlappableGenomicRange other) {
        return LocusUtils.partialGenomicRangeOverlap(this, other);
    }

    @Override
    @Deprecated
    public boolean isOverlapping(OverlappableGenomicRange other, int window) {
        return LocusUtils.partialGenomicRangeOverlapWindow(this, other, window);
    }

    public List<Gene> getGenes() {
        return genes;
    }
}
