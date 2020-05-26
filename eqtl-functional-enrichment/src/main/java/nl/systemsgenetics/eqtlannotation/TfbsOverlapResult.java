/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package nl.systemsgenetics.eqtlannotation;

/**
 *
 * @author MarcJan
 */
public class TfbsOverlapResult implements Comparable<TfbsOverlapResult>{
    private final String name;
    private final double pvalue;
    private final int nrOfTimesTfForSnp;
    private final int nrOfTimesSnpNoneTf;
    private final int nrOfTimesTfFOtherSnps;
    private final int currentOther;
    
    public TfbsOverlapResult(String name,double pvalue, int nrOfTimesTfForSnp,int nrOfTimesSnpNoneTf,int nrOfTimesTfFOtherSnps,int currentOther){
        this.name = name;
        this.pvalue = pvalue;
        this.nrOfTimesTfForSnp = nrOfTimesTfForSnp;
        this.nrOfTimesSnpNoneTf = nrOfTimesSnpNoneTf;
        this.nrOfTimesTfFOtherSnps = nrOfTimesTfFOtherSnps;
        this.currentOther = currentOther;
    }

    public String getName() {
        return name;
    }

    public double getPvalue() {
        return pvalue;
    }

    public int getNrOfTimesTfForSnp() {
        return nrOfTimesTfForSnp;
    }

    public int getNrOfTimesSnpNoneTf() {
        return nrOfTimesSnpNoneTf;
    }

    public int getNrOfTimesTfFOtherSnps() {
        return nrOfTimesTfFOtherSnps;
    }

    public int getCurrentOther() {
        return currentOther;
    }
    
    public double getRatio(){
        double percentage1  = nrOfTimesTfForSnp/(double)(nrOfTimesTfForSnp+nrOfTimesSnpNoneTf);
        double percentage2= nrOfTimesTfFOtherSnps/(double)(nrOfTimesTfFOtherSnps+currentOther);
        return percentage1/percentage2;
    }
    
    @Override
    public String toString(){
        double percentage1  = nrOfTimesTfForSnp/(double)(nrOfTimesTfForSnp+nrOfTimesSnpNoneTf);
        double percentage2= nrOfTimesTfFOtherSnps/(double)(nrOfTimesTfFOtherSnps+currentOther);
        double ratio = percentage1/percentage2;
        return(name+'\t'+nrOfTimesTfForSnp+'\t'+nrOfTimesSnpNoneTf+'\t'+nrOfTimesTfFOtherSnps+'\t'+currentOther+'\t'+pvalue+'\t'+(percentage1*100)+'\t'+(percentage2*100)+'\t'+ratio);
    }

    @Override
    public int compareTo(TfbsOverlapResult other) {
        if (other.getPvalue()> this.getPvalue()) {
            return -1;
        } else if (other.getPvalue() < this.getPvalue()) {
            return 1;
        } else {
            return 0;
        }
    }
    
}
