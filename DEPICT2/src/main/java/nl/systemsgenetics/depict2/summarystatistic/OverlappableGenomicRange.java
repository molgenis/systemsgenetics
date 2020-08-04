package nl.systemsgenetics.depict2.summarystatistic;

public interface OverlappableGenomicRange {

    int getStart();
    int getEnd();
    String getSequenceName();
    boolean isOverlapping(OverlappableGenomicRange other);
}
