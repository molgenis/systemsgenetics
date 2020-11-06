package nl.systemsgenetics.downstreamer.summarystatistic;

public interface OverlappableGenomicRange {

    int getStart();
    int getEnd();
    String getSequenceName();
    boolean isOverlapping(OverlappableGenomicRange other);
    boolean isOverlapping(OverlappableGenomicRange other, int window);

}
