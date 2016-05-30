/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package nl.systemsgenetics.eqtlannotation;

import java.util.regex.Pattern;

/**
 *
 * @author MarcJan
 */
public class SNP {
    private final String primaryId;
    private final String chr;
    private final long position;
    private final String[] alleles = new String[2];
    private final Integer[] moduleAssociation;
    private final String[] direction;
    static final Pattern SEP_PATTERN = Pattern.compile("/");
    
    SNP(String[] split) {
        primaryId = split[1]+"-"+split[4];
        chr = split[2];
        position = Long.parseLong(split[3]);
           
        String[] allelesTmp = SEP_PATTERN.split(split[8]);
        if(!allelesTmp[1].equals(split[9])){
            alleles[0] = allelesTmp[1];
            alleles[1] = allelesTmp[0];
        } else {
            alleles[0] = allelesTmp[0];
            alleles[1] = allelesTmp[1];
        }
        
        String[] modules  = split[12].split(",");
        moduleAssociation = new Integer[modules.length];
        direction = new String[modules.length];
        
        for(int i=0; i<modules.length; i++){
            if(modules[i].charAt(0)=='-'){
                direction[i]="-";
                modules[i] = modules[i].replaceFirst("-","");
            } else {
                direction[i]="+";
            }
            moduleAssociation[i] = Integer.parseInt(modules[i]);
        }
    }

    public String getPrimaryId() {
        return primaryId;
    }

    public String getChr() {
        return chr;
    }

    public long getPosition() {
        return position;
    }

    public String[] getAlleles() {
        return alleles;
    }
    
    public String getPositionId() {
        return (chr+":"+position);
    }

    public Integer[] getModuleAssociation() {
        return moduleAssociation;
    }

    public String[] getDirection() {
        return direction;
    }

    
}
