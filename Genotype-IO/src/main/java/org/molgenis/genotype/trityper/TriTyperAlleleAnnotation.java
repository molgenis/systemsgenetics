/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package org.molgenis.genotype.trityper;

import org.molgenis.genotype.Allele;

/**
 *
 * @author harmjan
 */
public class TriTyperAlleleAnnotation {

    private static final byte A = 65;
    private static final byte C = 67;
    private static final byte T = 84;
    private static final byte G = 71;
    private static final byte U = 85;
    private static final byte I = 73;
    private static final byte N = 78;
    
    public static Allele convertByteToAllele(byte b){
        switch(b){
            case A:
                return Allele.A;
            case C:
                return Allele.C;
            case T:
                return Allele.T;
            case G: 
                return Allele.G;
            case U:
                return Allele.create('U');
            case I:
                return Allele.create('I');
            case N:
                return Allele.create('N');
            default:
                return Allele.ZERO;
        }
    }
}
