/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package eqtlmappingpipeline.trityperconverter;

import umcg.genetica.io.trityper.converters.MachImputedToTriTyper;
import umcg.genetica.io.trityper.converters.MachImputedTransposedDosageToTriTyper;
import java.io.IOException;
import java.util.logging.Level;
import java.util.logging.Logger;
import umcg.genetica.console.ConsoleGUIElems;
import umcg.genetica.io.Gpio;
import umcg.genetica.io.text.TextFile;

/**
 *
 * @author harmjan
 */
public class TriTyperGenotypeConverter {

    public TriTyperGenotypeConverter(String[] args) {
        
        String in = null;
        String out = null;
        String type = "machtransposed";
        
        String header = null;
        for (int i = 0; i < args.length; i++) {
            String arg = args[i];
            String val = null;

            if (i + 1 < args.length) {
                val = args[i + 1];
            }

            if (arg.equals("--in")) {
                in = val;
            } else if (arg.equals("--out")) {
                out = val;
            } else if (arg.equals("--header")) {
                header = val;
            }
        }

        if (in == null || out == null) {
            System.err.println("Error: please supply both --in and --out locations");
            printUsage();
        } else {
            if (!Gpio.exists(in)) {
                System.err.println("Error: directory " + in + " does not exist!");
                printUsage();
            } else {
                if (type == null) {
                    try {
                        String filetype = determineFileType(in);

                        if (filetype == null) {
                            System.err.println("Error: no files in directory.");
                        } else if (filetype.equals("multi")) {
                            System.err.println("Error: multiple useable filetypes in directory. Please specify type to import using --type");
                        } else if (filetype.equals("beagle")) {
                            System.err.println("Error: BEAGLE files not supported (yet)");
                        } else if (filetype.equals("mach")) {
                            MachImputedToTriTyper mach = new MachImputedToTriTyper(in, out);
                        } else if (filetype.equals("machtransposed")) {
                            MachImputedTransposedDosageToTriTyper mach = new MachImputedTransposedDosageToTriTyper(in, out, header);
                        } else if (filetype.equals("Imputev2")) {
                            System.err.println("Error: ImputeV2 files not supported (yet)");
                        } else if (filetype.equals("ped")) {
                            System.err.println("Error: PED files not supported (yet)");
                        } else {
                            System.out.println("Filetype not recognized.");
                            printUsage();
                        }

                    } catch (IOException ex) {
                        Logger.getLogger(TriTyperGenotypeConverter.class.getName()).log(Level.SEVERE, null, ex);
                    }
                } else {
                    try {
                        
                        if(type.equals("machtransposed")){
                            MachImputedTransposedDosageToTriTyper mach = new MachImputedTransposedDosageToTriTyper(in, out, header);
                        } else {
                            System.out.println("Filetype not recognized.");
                            printUsage();
                        }
                    } catch (IOException ex) {
                        Logger.getLogger(TriTyperGenotypeConverter.class.getName()).log(Level.SEVERE, null, ex);
                    }
                }
            }
        }
    }

    /*
     * prints usage
     */
    private void printUsage() {
        System.out.println("");
        System.out.print("Command line options:\n"+ConsoleGUIElems.LINE);
        System.out.println("--in\t\tdir\tLocation of (imputed) genotype data");
        System.out.println("--out\t\tdir\tLocation to save converted data");
        System.out.println("--header\tstring\tFull path to header data file");
//        System.out.println("[--type] [pedmap|mach|machtransposed|beagle|impute]\tLocation to save converted data\n");
        System.out.println("");
    }

    private String determineFileType(String in) throws IOException {
        String[] fileList = Gpio.getListOfFiles(in);
        if (fileList.length == 0) {
            return null;
        } else {
            boolean multi = false;
            String fileType = null;
            for (String filename : fileList) {
                if (filename.endsWith("mldose")) {
                    if (fileType == null) {
                        fileType = "mach";
                    } else if (!fileType.equals("mach")) {
                        multi = true;
                    }
                } else if (filename.endsWith("ped")) {
                    if (fileType == null) {
                        fileType = "ped";
                    } else if (!fileType.equals("ped")) {
                        multi = true;
                    }
                } else if (filename.endsWith("gprobs.gz")) {
                    if (fileType == null) {
                        fileType = "beagle";
                    } else if (!fileType.equals("beagle")) {
                        multi = true;
                    }
                } else {
                    // try to parse the first line of the file and determine the filetype..
                    TextFile t = new TextFile(filename, TextFile.R);
                    
                }
            }

            if (multi) {
                return "multi";
            } else {
                return fileType;
            }
        }

    }
}
