/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package umcg.genetica.io.ucsc;

import java.util.ArrayList;
import java.util.Collections;

/**
 *
 * @author harmjan
 */
public class UCSCTrack implements Comparable<UCSCTrack> {

    private String name;
    private ArrayList<UCSCDataObject> data;
    private double min, max;
    private double sd;
    private double mean;

    public UCSCTrack(String name) {
        this.name = name;
        data = new ArrayList<UCSCDataObject>();
    }

    public void addDataObject(UCSCDataObject o) {
        data.add(o);
    }

    public void setMax(double max) {
        this.max = max;
    }

    public void setMin(double min) {
        this.min = min;
    }

    public void sort() {
        Collections.sort(data);
    }

    public void setSD(double sd) {
        this.sd = sd;
    }

    public void setMean(double mean) {
        this.mean = mean;
    }

    @Override
    public int compareTo(UCSCTrack o) {
        if (o.equals(this)) {
            return 0;
        } else {
            return name.compareTo(o.name);
        }
    }

    public boolean equals(UCSCTrack o) {
        if (o.name.equals(name)) {
            return true;
        } else {
            return false;
        }
    }
    
    public String getName(){
        return name;
    }
    
    public ArrayList<UCSCDataObject> getData(){
        return data;
    }
    
    public double getMean(){
        return mean;
    }
    
    public double getSD(){
        return sd;
    }
}
