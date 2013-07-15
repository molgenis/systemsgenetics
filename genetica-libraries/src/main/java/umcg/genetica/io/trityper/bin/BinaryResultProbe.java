/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package umcg.genetica.io.trityper.bin;

/**
 *
 * @author harmjan
 */
public class BinaryResultProbe {
    private int id;
    private String name;
    private Integer start;
    private Integer stop;
    private Byte chr;
    private Integer midpoint;
    private String annotation;

    /**
     * @return the id
     */
    public int getId() {
        return id;
    }

    /**
     * @return the name
     */
    public String getName() {
        return name;
    }

    /**
     * @return the start
     */
    public int getStart() {
        return start;
    }

    /**
     * @return the stop
     */
    public int getStop() {
        return stop;
    }

    /**
     * @return the chr
     */
    public byte getChr() {
        return chr;
    }

    /**
     * @return the midpoint
     */
    public int getMidpoint() {
        return midpoint;
    }

    /**
     * @return the annotation
     */
    public String getAnnotation() {
        return annotation;
    }

    /**
     * @param id the id to set
     */
    public void setId(int id) {
        this.id = id;
    }

    /**
     * @param name the name to set
     */
    public void setName(String name) {
        this.name = name;
    }

    /**
     * @param start the start to set
     */
    public void setStart(int start) {
        this.start = start;
    }

    /**
     * @param stop the stop to set
     */
    public void setStop(int stop) {
        this.stop = stop;
    }

    /**
     * @param chr the chr to set
     */
    public void setChr(byte chr) {
        this.chr = chr;
    }

    /**
     * @param midpoint the midpoint to set
     */
    public void setMidpoint(int midpoint) {
        this.midpoint = midpoint;
    }

    /**
     * @param annotation the annotation to set
     */
    public void setAnnotation(String annotation) {
        this.annotation = annotation;
    }

    public void clearMetaData() {
	throw new UnsupportedOperationException("Not yet implemented");
    }

}
