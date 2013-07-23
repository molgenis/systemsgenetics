/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package org.molgenis.genotype.annotation;

/**
 *
 * @author harmjan
 */
public enum CaseControlAnnotation {

    CASE((byte) 2), CONTROL((byte) 1), UNKNOWN((byte) 0);
    private final byte plinkCase;

    private CaseControlAnnotation(byte plinkCase) {
        this.plinkCase = plinkCase;
    }

    /**
     * @return the plinkSex
     */
    public byte getPlinkCase() {
        return plinkCase;
    }

    public static CaseControlAnnotation getCaseAnnotationForPlink(byte plinkCase) {
        switch (plinkCase) {
            case 1:
                return CaseControlAnnotation.CONTROL;
            case 2:
                return CaseControlAnnotation.CASE;
            default:
                return CaseControlAnnotation.UNKNOWN;
        }
    }

    public static CaseControlAnnotation getCaseAnnotationForTriTyper(String ttCase) {
        if (ttCase == null) {
            return CaseControlAnnotation.UNKNOWN;
        } else if (ttCase.toLowerCase().equals("control")) {
            return CONTROL;
        } else if (ttCase.toLowerCase().equals("case")) {
            return CASE;
        } else {
            return CaseControlAnnotation.UNKNOWN;
        }
    }
}
