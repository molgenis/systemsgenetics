/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package nl.systemsgenetics.simplegeneticriskscorecalculator;

/**
 *
 * @author MarcJan
 */
class RiskEntry implements Comparable<RiskEntry>  {
    private final String rsName;
    private final double or;
    private final int pos;
    private final int chr;
    private final char allele;
    private final double pValue;

    RiskEntry(String rsName, int chr, int pos, char allele, String or, double pValue) {
        this.rsName = rsName;
        this.chr = chr;
        this.pos = pos;
        this.allele = allele;
        this.pValue = pValue;
        this.or = Double.parseDouble(or);
    }
    
    RiskEntry(String rsName, String chr, int pos, char allele, String or, double pValue) {
        this.rsName = rsName;
        this.pos = pos;
        this.allele = allele;
        this.pValue = pValue;
        this.or = Double.parseDouble(or);
        if(chr.equals("X")){
            this.chr = 23;
        } else if(chr.equals("Y")){
            this.chr = 24;
        } else {
            this.chr=Integer.parseInt(chr);
        }
    }

    public String getRsName() {
        return rsName;
    }

    public double getOr() {
        return or;
    }

    public int getPos() {
        return pos;
    }

    public String getChr() {
        return String.valueOf(chr);
    }

    public char getAllele() {
        return allele;
    }

    public double getpValue() {
        return pValue;
    }

    @Override
    public int compareTo(RiskEntry o) {
        if(this.pValue<o.getpValue()){
            return -1;
        } else if (this.pValue == o.getpValue()) {
            if(Math.abs(this.or)>Math.abs(o.getOr())){
                return -1;
            } else if (Math.abs(this.or) == Math.abs(o.getOr())) {
                return 0;
            } else {
                return 1;
            }
        } else {
            return 1;
        }
    }

    String InfoToString() {
        StringBuilder s = new StringBuilder();
        s.append(this.rsName).append("\t");
        s.append(this.chr).append("\t");
        s.append(this.pos).append("\t");
        s.append(this.allele).append("\t");
        s.append(this.or).append("\t");
        s.append(this.pValue);
        return(s.toString());
    }
    
    
    
}
