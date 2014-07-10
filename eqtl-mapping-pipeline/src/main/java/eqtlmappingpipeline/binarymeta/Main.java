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
import eqtlmappingpipeline.metaqtl3.FDR;
import eqtlmappingpipeline.util.NoLdSnpProbeListCreator;
import java.io.FileNotFoundException;
import java.io.IOException;
import java.io.UnsupportedEncodingException;
import java.util.Arrays;
import java.util.logging.Level;
import java.util.logging.Logger;

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
//        m.reannotateprobesInEQTLFile("D:\\Work\\MetaTest\\2011-10-06-ProbeTranslationTable+H8HT12Conversion.log_reannotatedHG18_96PercIdentity.txt", "Z:\\MarjoleinHomeAccount\\marjolein\\Results\\2011-10-25-METATEST.Rotterdam+EGCUT+Groningen+SHIP-40PCs.YES-GVR.YES-GWAS-PCS\\eQTLsFDR0.05.txt");
//        } catch (Exception e){
//            e.printStackTrace();
//        }


        String texttoreplace = null;
        String replacetextwith = null;
        boolean individualAnalysis = false;
        String mode = null;
        String settings = null;
        String gwascatalog = null;
        String cormat = null;
        String transfile = null;
        String cisfile = null;
        String probetranslation = null;
        String zscore = null;
        String out = null;

        String eqtlfile = null;
        String probeannot = null;

        String annot = null;
        String in = null;

        Integer nrPerm = 0;
        Integer nrEQTLs = null;
        Double cutoff = null;
        String query = null;

        String ldfile = null;
        String dbsnp = null;
        String pw = null;
        String snpselectionlist = null;
        String snpprobeselectionlist = null;
        boolean createQQPlot = true;

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
            } else if (args[i].equals("--gwascatalog")) {
                gwascatalog = val;
            } else if (args[i].equals("--trans")) {
                transfile = val;
            } else if (args[i].equals("--cis")) {
                cisfile = val;
            } else if (args[i].equals("--cormat")) {
                cormat = val;
            } else if (args[i].equals("--probetranslation")) {
                probetranslation = val;
            } else if (args[i].equals("--zscore")) {
                zscore = val;
            } else if (args[i].equals("--out")) {
                out = val;
            } else if (args[i].equals("--probeannot")) {
                probeannot = val;
            } else if (args[i].equals("--eqtlfile")) {
                eqtlfile = val;
            } else if (args[i].equals("--in")) {
                in = val;
            } else if (args[i].equals("--annot")) {
                annot = val;
            } else if (args[i].equals("--nrperm")) {
                nrPerm = Integer.parseInt(val);
            } else if (args[i].equals("--cutoff")) {
                cutoff = Double.parseDouble(val);
            } else if (args[i].equals("--nreqtls")) {
                nrEQTLs = Integer.parseInt(val);
            } else if (args[i].equals("--query")) {
                query = val;
            } else if (args[i].equals("--ldfile")) {
                ldfile = val;
            } else if (args[i].equals("--dbsnp")) {
                dbsnp = val;
            } else if (args[i].equals("--pw")) {
                pw = val;
            } else if (args[i].equals("--individual")) {
                individualAnalysis = true;
            } else if (args[i].equals("--texttoreplace")) {
                texttoreplace = val;
            } else if (args[i].equals("--replacetextwith")) {
                replacetextwith = val;
            } else if (args[i].equals("--skipqqplot")) {
                createQQPlot = false;
            } else if (args[i].equals("--snpselectionlist")) {
                snpselectionlist = val;
            } else if (args[i].equals("--snpprobeselectionlist")) {
                snpprobeselectionlist = val;
            }


        }

        if (mode == null) {
            System.out.println("Specify metamode (meta, cismeta, summary, splitzscoretable, allelecheck, filter or fdr, crosshybparser, determineSnpProbList)");
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
                System.out.println("Specify mode (meta or summary or splitzscoretable, or reannotate)");
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
        } else if (mode.equals("fdr")) {
            if (in == null || nrEQTLs == null || cutoff == null) {
                System.out.println("Please specify --in --nrperm and --cutoff and --nreqtls [--skipqqplot]");
            } else {
                if(snpselectionlist!=null){
                    try {
                        FDR.calculateFDR2(in, nrPerm, nrEQTLs, cutoff, createQQPlot, null, null, FDR.FDRMethod.ALL, true, snpselectionlist, snpprobeselectionlist);
                    } catch (IOException e) {
                        e.printStackTrace();
                    }
                } else {
                    try {
                        FDR.calculateFDR(in, nrPerm, nrEQTLs, cutoff, createQQPlot, null, null, FDR.FDRMethod.ALL, true);
                    } catch (IOException e) {
                        e.printStackTrace();
                    }
                }
            }
        } else if (mode.equals("ld")) {
            if (in == null || eqtlfile == null || out == null) {
                System.out.println("Please specify --in and --eqtlfile and --outdir");
            } else {
                try {
                    LDCalc ld = new LDCalc();
                    ld.LDCalc(eqtlfile, in, out);
                } catch (IOException e) {
                    e.printStackTrace();
                } catch (Exception e) {
                    e.printStackTrace();
                }
            }

        } else if (mode.equals("determineSnpProbList")) {
            try {
                NoLdSnpProbeListCreator.main(Arrays.copyOfRange(args, 2, args.length));
            } catch (UnsupportedEncodingException ex) {
                Logger.getLogger(Main.class.getName()).log(Level.SEVERE, null, ex);
            } catch (FileNotFoundException ex) {
                Logger.getLogger(Main.class.getName()).log(Level.SEVERE, null, ex);
            } catch (IOException ex) {
                Logger.getLogger(Main.class.getName()).log(Level.SEVERE, null, ex);
            } catch (Exception ex) {
                Logger.getLogger(Main.class.getName()).log(Level.SEVERE, null, ex);
            }
        } else {
            System.out.print("Invalid option, valid options are:");
            System.out.println("fdr, ld, filter, allelecheck, individual, cismeta, meta, determineSnpProbList");
        }

        System.exit(0);
    }
}
