/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package umcg.genetica.io.trityper.util;

/**
 *
 * @author harm-jan
 */
public class ChrAnnotation {
    public static byte parseChr(String chrStr){
        byte chr = -1;
        try{
             chr = Byte.parseByte(chrStr);
        } catch (NumberFormatException e){
            if (chrStr.toUpperCase().equals("X")) {
                chr = 23;   //X chromosome is 23
            } else if (chrStr.toUpperCase().equals("Y")) {
                chr = 24;   //Y chromosome is 24
            } else if (chrStr.toUpperCase().equals("XY")) {
                chr = 25;   //XY chromosome is 25
            } else if (chrStr.toUpperCase().equals("M") || chrStr.toUpperCase().equals("MT") ) {
                chr = 26;   //M chromosome is 26
//            } else if (chrStr.toUpperCase().equals("CHRM")) {
//                chr = 26;   //Mitochondrial chromosome is 26
//            } else if (chrStr.toUpperCase().equals("MT")) {
//                chr = 26;   //Mitochondrial chromosome is 26
//            } else if (chrStr.equals("c6_COX")) {
//                chr = 6; //HLA region on chromosome 6
            } else {
                chr = -1;
                if(!chrStr.equals("-")){
                    System.err.println("Could not parse chr:\t"+chrStr);
                }
            }

        }
        return chr;
    }

    public static String parseByte(byte chr){
        if(chr < 23){
            return ""+chr;
        } else {
            if(chr == 23){
                return "X";
            }
            if(chr == 24){
                return "Y";
            }

            if(chr == 25){
                return "XY";
            }

            if(chr == 26){
                return "MT";
            }
        }

        return "-1";
    }
}
