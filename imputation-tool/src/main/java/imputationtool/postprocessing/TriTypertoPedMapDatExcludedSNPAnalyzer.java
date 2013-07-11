/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package imputationtool.postprocessing;

import java.io.IOException;
import umcg.genetica.io.text.TextFile;

/**
 *
 * @author harmjan
 */
public class TriTypertoPedMapDatExcludedSNPAnalyzer {

    public void printSNPsWithCallRateHigherThan(double threshold, String location) throws IOException {

        for (int i = 1; i < 23; i++) {
            TextFile in = new TextFile(location + "chr" + i + ".excludedsnps.txt", TextFile.R);
            String line = "";

            while ((line = in.readLine()) != null) {
                String[] elems = line.split("\t");

                if (elems.length > 1) {

                    String snpName = elems[0];
                    String snpProperties = elems[2];
                    // System.out.println(snpName +" \t " + snpProperties );
                    String[] separateSNPProperties = snpProperties.split(" ");

                    if (separateSNPProperties.length > 0) {
                        for (int q = 0; q < separateSNPProperties.length; q++) {
                            String[] properties = snpProperties.split("=");
                            if (properties[0].equals("CallRate")) {
                                double callrate = Double.parseDouble(properties[1]);
                                if (callrate >= threshold) {
                                    System.out.println(snpName);
                                }
                            }
                        }

                    } else {

                        String[] properties = snpProperties.split("=");
                        if (properties[0].equals("CallRate")) {
                            double callrate = Double.parseDouble(properties[1]);
                            if (callrate >= threshold) {
                                System.out.println(snpName);
                            }
                        }


                    }

                }

            }
            in.close();
        }


    }
}
