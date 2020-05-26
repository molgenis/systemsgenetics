/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package eqtlmappingpipeline.binarymeta;

import eqtlmappingpipeline.binarymeta.util.Filter;
import eqtlmappingpipeline.binarymeta.meta.IndividualAnalysis;
import eqtlmappingpipeline.binarymeta.meta.MetaAnalyze;
import eqtlmappingpipeline.binarymeta.meta.cis.CisAnalysis;
import eqtlmappingpipeline.binarymeta.util.SNPAlleleCheck;
import eqtlmappingpipeline.util.LDCalculator;
import java.io.IOException;

/**
 *
 * @author harm-jan
 */
public class Main {

    /**
     * @param args the command line arguments
     */
    public static void main(String[] args) {
//        try{
//        MetaAnalyze m = new MetaAnalyze();
//
//        m.reannotateprobesInEQTLFile("D:\\Work\\MetaTest\\2011-10-06-ProbeTranslationTable+H8HT12Conversion.log_reannotatedHG18_96PercIdentity.txt", "Z:\\MarjoleinHomeAccount\\marjolein\\Results\\2011-10-25-METATEST.Rotterdam+EGCUT+Groningen+SHIP-40PCs.YES-GVR.YES-GWAS-PCS\\eQTLsFDR0.05.txt.gz");
//        } catch (Exception e){
//            e.printStackTrace();
//        }


        String texttoreplace = null;
        String replacetextwith = null;
        boolean individualAnalysis = false;
        String mode = null;
        String settings = null;
        String probetranslation = null;
        String out = null;

        String eqtlfile = null;

        String annot = null;
        String in = null;

        for (int i = 0; i < args.length; i++) {
            String val = null;
            if (i + 1 < args.length) {
                val = args[i + 1];
            }
            if (args[i].equals("--metamode")) {
                mode = val;
            } else if (args[i].equals("--settings")) {
                settings = val;
            } else if (args[i].equals("--settings")) {
            } else if (args[i].equals("--probetranslation")) {
                probetranslation = val;
            } else if (args[i].equals("--out")) {
                out = val;
            } else if (args[i].equals("--eqtlfile")) {
                eqtlfile = val;
            } else if (args[i].equals("--in")) {
                in = val;
            } else if (args[i].equals("--individual")) {
                individualAnalysis = true;
            } else if (args[i].equals("--texttoreplace")) {
                texttoreplace = val;
            } else if (args[i].equals("--replacetextwith")) {
                replacetextwith = val;
            } 


        }

        if (mode == null) {
            System.out.println("Specify metamode (meta, cismeta, individual, allelecheck, filter or ld)");
        } else if (mode.equals("meta")) {
            if (settings == null) {
                System.out.println("Specify settings");
            } else {
                try {

                    MetaAnalyze m2 = new MetaAnalyze();
                     m2.init(settings, texttoreplace, replacetextwith);
                     m2.analyze();

//                    System.gc();
//                    System.gc();
//                    System.gc();
//                    System.gc();
//                    System.gc();


                    if (individualAnalysis) {
                        IndividualAnalysis m = new IndividualAnalysis();
                        m.init(settings, texttoreplace, replacetextwith);
                        m.analyze();
                        m = null;
                    }

                } catch (Exception e) {
                    e.printStackTrace();
                }
            }
        } else if (mode.equals("cismeta")) {
            if (settings == null) {
                System.out.println("Specify settings");
            } else {
                try {

                    CisAnalysis c = new CisAnalysis(settings);
                    c.analyze();

                } catch (Exception e) {
                    e.printStackTrace();
                }
            }
        } else if (mode.equals("individual")) {
            try {
                IndividualAnalysis m = new IndividualAnalysis();
                m.init(settings, texttoreplace, replacetextwith);
                m.analyze();
                m = null;
            } catch (Exception e) {
                e.printStackTrace();
            }

        } else if (mode.equals("allelecheck")) {
            if (settings == null) {
                System.out.println("Specify settings");
            } else {
                try {
                    if (true) {
                        SNPAlleleCheck m2 = new SNPAlleleCheck();
                        m2.init(settings, texttoreplace, replacetextwith);
                        m2.analyze();

//                        System.gc();
//                        System.gc();
//                        System.gc();
//                        System.gc();
//                        System.gc();
                    }
                } catch (Exception e) {
                    e.printStackTrace();
                }
            }
        } else if (mode.equals("filter")) {
            if (in == null || probetranslation == null || annot == null) {
                System.out.println("Please specify --in --probetranslation and --annot");
            } else {
                try {
                    Filter f = new Filter();
                    f.run(in, probetranslation, annot);
                } catch (IOException e) {
                    e.printStackTrace();
                }
            }
        }  else if (mode.equals("ld")) {
            if (in == null || eqtlfile == null || out == null) {
                System.out.println("Please specify --in and --eqtlfile and --outdir");
            } else {
                try {
                    LDCalculator ld = new LDCalculator();
                    ld.calculatePairwiseLD(eqtlfile, in, out);
                } catch (IOException e) {
                    e.printStackTrace();
                } catch (Exception e) {
                    e.printStackTrace();
                }
            }

        }  else {
            System.out.print("Invalid option, valid options are:");
            System.out.println("meta, cismeta, individual, allelecheck, filter or ld");
        }

        System.exit(0);
    }
}
