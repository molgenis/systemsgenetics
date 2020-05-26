/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package umcg.genetica.enums;

/**
 *
 * @author Harm-Jan
 */
public enum ChromosomeUCSC {

    ONE(1, "Chr1", 248956422),
    TWO(2, "Chr2", 242193529),
    THREE(3, "Chr3", 198295559),
    FOUR(4, "Chr4", 19021455),
    FIVE(5, "Chr5", 181538259),
    SIX(6, "Chr6", 170805979),
    SEVEN(7, "Chr7", 159345973),
    EIGHT(8, "Chr8", 145138636),
    NINE(9, "Chr9", 138394717),
    TEN(10, "Chr10", 133797422),
    ELEVEN(11, "Chr11", 135086622),
    TWELVE(12, "Chr12", 133275309),
    THIRTEEN(13, "Chr13", 114364328),
    FOURTEEN(14, "Chr14", 107043718),
    FIFTEEN(15, "Chr15", 101991189),
    SIXTEEN(16, "Chr16", 90338345),
    SEVENTEEN(17, "Chr17", 83257441),
    EIGHTEEN(18, "Chr18", 80373285),
    NINETEEN(19, "Chr19", 58617616),
    TWENTY(20, "Chr20", 64444167),
    TWENTYONE(21, "Chr21", 46709983),
    TWENTYTWO(22, "Chr22", 50818468),
    X(23, "ChrX", 156040895),
    Y(24, "ChrY", 57227415),
    MT(25, "ChrMT", 1),
    NA(0, "N/A", 1);

    private final int number;
    private final String name;
    private final int length;

    private ChromosomeUCSC(int num, String name, int length) {
        this.number = num;
        this.name = name;
        this.length = length;
    }

    public String getName() {
        return name;
    }

    public int getLength() {
        return length;
    }
    
     public int getNumber() {
        return number;
    }

    public static ChromosomeUCSC parseChr(String chrStr) {
        chrStr = chrStr.toLowerCase().trim();
        if (chrStr.equals("chr1")) {
            return ChromosomeUCSC.ONE;
        }
        if (chrStr.equals("chr2")) {
            return ChromosomeUCSC.TWO;
        }
        if (chrStr.equals("chr3")) {
            return ChromosomeUCSC.THREE;
        }
        if (chrStr.equals("chr4")) {
            return ChromosomeUCSC.FOUR;
        }
        if (chrStr.equals("chr5")) {
            return ChromosomeUCSC.FIVE;
        }
        if (chrStr.equals("chr6")) {
            return ChromosomeUCSC.SIX;
        }
        if (chrStr.equals("chr7")) {
            return ChromosomeUCSC.SEVEN;
        }
        if (chrStr.equals("chr8")) {
            return ChromosomeUCSC.EIGHT;
        }
        if (chrStr.equals("chr9")) {
            return ChromosomeUCSC.NINE;
        }
        if (chrStr.equals("chr10")) {
            return ChromosomeUCSC.TEN;
        }
        if (chrStr.equals("chr11")) {
            return ChromosomeUCSC.ELEVEN;
        }
        if (chrStr.equals("chr12")) {
            return ChromosomeUCSC.TWELVE;
        }
        if (chrStr.equals("chr13")) {
            return ChromosomeUCSC.THIRTEEN;
        }
        if (chrStr.equals("chr14")) {
            return ChromosomeUCSC.FOURTEEN;
        }
        if (chrStr.equals("chr15")) {
            return ChromosomeUCSC.FIFTEEN;
        }
        if (chrStr.equals("chr16")) {
            return ChromosomeUCSC.SIXTEEN;
        }
        if (chrStr.equals("chr17")) {
            return ChromosomeUCSC.SEVENTEEN;
        }
        if (chrStr.equals("chr18")) {
            return ChromosomeUCSC.EIGHTEEN;
        }
        if (chrStr.equals("chr19")) {
            return ChromosomeUCSC.NINETEEN;
        }
        if (chrStr.equals("chr20")) {
            return ChromosomeUCSC.TWENTY;
        }
        if (chrStr.equals("chr21")) {
            return ChromosomeUCSC.TWENTYONE;
        }
        if (chrStr.equals("chr22")) {
            return ChromosomeUCSC.TWENTYTWO;
        }
        if (chrStr.equals("chry")) {
            return ChromosomeUCSC.Y;
        }
        if (chrStr.equals("chrx")) {
            return ChromosomeUCSC.X;
        }

        return ChromosomeUCSC.NA;
    }
}
