package nl.systemsgenetics.geneticriskscorecalculator;

import gnu.trove.map.hash.TObjectDoubleHashMap;
import java.io.File;
import java.io.IOException;
import java.text.DateFormat;
import java.text.SimpleDateFormat;
import java.util.ArrayList;
import java.util.Date;
import java.util.Iterator;
import java.util.List;
import java.util.ResourceBundle;
import java.util.regex.Pattern;
import org.apache.commons.cli.ParseException;
import org.apache.log4j.FileAppender;
import org.apache.log4j.Level;
import org.apache.log4j.Logger;
import org.apache.log4j.SimpleLayout;
import org.molgenis.genotype.Alleles;
import org.molgenis.genotype.GenotypeDataException;
import org.molgenis.genotype.RandomAccessGenotypeData;
import org.molgenis.genotype.multipart.IncompatibleMultiPartGenotypeDataException;
import org.molgenis.genotype.variant.GeneticVariant;
import umcg.genetica.io.text.TextFile;

/**
 * Hello world!
 *
 */
public class Main {

    //private static final String VERSION = ResourceBundle.getBundle("verion").getString("application.version");
    private static final String VERSION = "0.0.1";

    private static final String HEADER
            = "  /---------------------------------------\\\n"
            + "  |     Genetic Risk Score Calculator     |\n"
            + "  |                                       |\n"
            + "  |              Rudi Alberts             |\n"
            + "  |                                       |\n"
            + "  |             Patrick Deelen            |\n"
            + "  |                                       |\n"
            + "  |                                       |\n"
            + "  | Dasha Zhernakova, Marijke v/d Sijde,  |\n"
            + "  |   Marc Jan Bonder, Harm-Jan Westra,   |\n"
            + "  |      Lude Franke, Morris Swertz       |\n"
            + "  |                                       |\n"
            + "  |     Genomics Coordication Center      |\n"
            + "  |        Department of Genetics         |\n"
            + "  |  University Medical Center Groningen  |\n"
            + "  \\---------------------------------------/";
    private static final Logger LOGGER = Logger.getLogger(Main.class);
    private static final DateFormat DATE_TIME_FORMAT = new SimpleDateFormat("yyyy-MM-dd HH:mm:ss");
    private static final Date currentDataTime = new Date();
    private static final Pattern TAB_PATTERN = Pattern.compile("\\t");

    public static void main(String[] args) {

        System.out.println(HEADER);
        System.out.println();
        System.out.println("          --- Version: " + VERSION + " ---");
        System.out.println();
        System.out.println("More information: http://molgenis.org/systemsgenetics");
        System.out.println();

        System.out.println("Current date and time: " + DATE_TIME_FORMAT.format(currentDataTime));
        System.out.println();

        System.out.flush(); //flush to make sure header is before errors
        try {
            Thread.sleep(25); //Allows flush to complete
        } catch (InterruptedException ex) {
        }

        if (args.length == 0) {
            Configuration.printHelp();
            System.exit(1);
        }

        final Configuration configuration;
        try {
            configuration = new Configuration(args);
        } catch (ParseException ex) {
            System.err.println("Invalid command line arguments: ");
            System.err.println(ex.getMessage());
            System.err.println();
            Configuration.printHelp();
            System.exit(1);
            return;
        }

        final File outputFile = configuration.getOutputFile();
        final File outputFolder = outputFile.getParentFile();
        if (outputFolder != null && !outputFolder.exists()) {
            if (!outputFolder.mkdirs()) {
                System.err.println("Failed to create ouput folder at: " + outputFolder.getAbsolutePath());
                System.exit(1);
            }
        }

        final File logFile = new File(outputFile.getAbsolutePath() + ".log");
        startLogging(logFile, true);

        final RandomAccessGenotypeData inputGenotypes;
        try {
            inputGenotypes = configuration.getInputDataType().createGenotypeData(configuration.getInputPaths(), 100);
            System.out.println("Loading reference data complete");
            LOGGER.info("Loading reference data complete");
        } catch (IOException ex) {
            System.err.println("Unable to load reference genotypes file.");
            LOGGER.fatal("Unable to load reference genotypes file.", ex);
            System.exit(1);
            return;
        } catch (IncompatibleMultiPartGenotypeDataException ex) {
            System.err.println("Unable to load reference genotypes file.");
            LOGGER.fatal("Unable to load reference genotypes file.", ex);
            System.exit(1);
            return;
        } catch (GenotypeDataException ex) {
            System.err.println("Unable to load reference genotypes file.");
            LOGGER.fatal("Unable to load reference genotypes file.", ex);
            System.exit(1);
            return;
        }
        int debug = 0;

        int index = 0;
        int index2 = 0;
        for (String mySample : inputGenotypes.getSampleNames()) {
            if (index < 5) {
                System.out.println("Samplename: " + mySample);
            }
            index++;
        }

        String fileLine;
        String[] fileLineData;
        TextFile riskSnpsFile;

        //List<Integer> chr; // = new List<Integer>();
        List<Integer> chr = new ArrayList<Integer>();
        List<Integer> pos = new ArrayList<Integer>();

        List<String> rsid = new ArrayList<String>();
        List<String> riskallele = new ArrayList<String>();

        index = 0;
        try {
            riskSnpsFile = new TextFile(configuration.getRisksnpsFile().toString(), false);
            riskSnpsFile.readLine(); // reading the header
            while ((fileLine = riskSnpsFile.readLine()) != null) {
                fileLineData = TAB_PATTERN.split(fileLine);
                if (index < 5) {
                    System.out.println("--" + fileLineData[1] + "--" + fileLineData[2] + "--" + fileLineData[3]);
                }
                if (index < 5) {
                    System.out.println("RiskAllele: " + fileLineData[4]);
                }
                chr.add(Integer.valueOf(fileLineData[1]));
                pos.add(Integer.valueOf(fileLineData[2]));
                rsid.add(fileLineData[3]);
                riskallele.add(fileLineData[4]);
                index++;

            }
            riskSnpsFile.close();
        } catch (IOException ex) {
            System.err.println("Unable to load risk snps file.");
            LOGGER.fatal("Unable to load risk snps file.", ex);
            System.exit(1);
            return;
        }
        //catch (NumberFormatException ex) {
//             System.err.println("Check your number format.");           
//            
//        }

        for (index = 0; index < 5; index++) {
            System.out.println("Risk snps:");
            System.out.println(" chr: " + chr.get(index) + " pos: " + pos.get(index) + " rs id: " + rsid.get(index) + " riskallele: " + riskallele.get(index));
        }

        System.out.println("amount of risk alleles read: " + chr.size());

        // we get alleles for all samples per snp
        // this is the main loop over all snps
        // per snp all allels per individual are collected
        //GeneticVariant snpVariantByPos;
        // test an SNP
        System.out.println("test snp ---");
        GeneticVariant snpVariantByPos = inputGenotypes.getSnpVariantByPos("1", 161012760);

        System.out.println("MINORALLELE  " + snpVariantByPos.getMinorAllele() + "  ALLALLELES  " + snpVariantByPos.getVariantAlleles() + "  IS GC AT snp: " + snpVariantByPos.getVariantAlleles().isAtOrGcSnp() + "  MAF  "
                + snpVariantByPos.getMinorAlleleFrequency());
        System.out.println("test snp ---");
        System.out.println();

        System.out.println("Amount of risk snps: " + chr.size());

        List<Integer> score = new ArrayList<Integer>();
        for (index = 0; index < inputGenotypes.getSampleNames().length; index++) {
            score.add(0);
        }

        int nrGCAT = 0;
        int nrnonGCAT = 0;
        
        for (index = 0; index < chr.size(); index++) {

            Integer currentChr = chr.get(index);
            int currentPos = Integer.valueOf(pos.get(index)); // necessary?

            System.out.print("SNP" + index + " " + rsid.get(index) + "  " + chr.get(index) + "  " + pos.get(index) + "  " + riskallele.get(index) + "  ");

            try {
                //snpVariantByPos = inputGenotypes.getSnpVariantByPos("1", 161012760);
                snpVariantByPos = inputGenotypes.getSnpVariantByPos(currentChr.toString(), currentPos);

                System.out.print("MINORALLELE  " + snpVariantByPos.getMinorAllele() + "  ALLALLELES  " + snpVariantByPos.getVariantAlleles() + "  Is GC AT snp: " + snpVariantByPos.getVariantAlleles().isAtOrGcSnp() + "  MAF  "
                        + snpVariantByPos.getMinorAlleleFrequency() + " ");

                if (snpVariantByPos.getVariantAlleles().isAtOrGcSnp()) {
                    nrGCAT++;
                } else {
                    nrnonGCAT++;
                }
                
                List<Alleles> alleles = snpVariantByPos.getSampleVariants();
                float[] dosages = snpVariantByPos.getSampleDosages();
                List<String> myids = snpVariantByPos.getAllIds();

                System.out.print(myids.get(0) + " ");
                index2 = 0;
                for (Alleles allele : alleles) {
                    if (index2 < 20) {
                        System.out.print(allele.toString() + " ");  // only alleles
                    }
                    if (snpVariantByPos.getVariantAlleles().isAtOrGcSnp()==false && ( riskallele.get(index).equals(allele.getAllelesAsString().get(0)) || riskallele.get(index).equals(allele.getAllelesAsString().get(1))) ) {
                        score.set(index2, score.get(index2) + 1);
                    }
                    index2++;
                }
                System.out.println();

            } catch (NullPointerException ex) {
                System.out.println("null pointer exception");

            }

        }  
        
        System.out.println("Number of non GC AT snps: " + nrnonGCAT + "  number of GC AT snps: " + nrGCAT);
        for (index = 0; index < inputGenotypes.getSampleNames().length; index++) {
            System.out.println(inputGenotypes.getSampleNames()[index] + " " + score.get(index));

        }


        GwasCatalogLoader gwasCatalogLoader = new GwasCatalogLoader();
        List<GeneticRiskScoreCalculator> geneticRiskScoreCalculators = gwasCatalogLoader.getGeneticRiskScoreCalculators(configuration.getRisksnpsFile().toString());

		List<String> phenotypes = new ArrayList<String>();
		for(GeneticRiskScoreCalculator calculator : geneticRiskScoreCalculators){
			phenotypes.add(calculator.getPhenotype());
		}
		
		RiskScoreMatrix riskScoreMatrix = new RiskScoreMatrix(inputGenotypes.getSampleNames(), phenotypes);
		
		for(GeneticRiskScoreCalculator calculator : geneticRiskScoreCalculators){
			TObjectDoubleHashMap<String> riskScores = calculator.calculateRiskScores(inputGenotypes);
			for(String sample : inputGenotypes.getSampleNames()){
				riskScoreMatrix.setRiskScore(sample, calculator.getPhenotype(), riskScores.get(sample));
			}
		}
		try {
			riskScoreMatrix.save(outputFile);
		} catch (IOException ex) {
			System.err.println("Could not save output file to: " + outputFile.getAbsolutePath());
			LOGGER.fatal("Could not save output file to: " + outputFile.getAbsolutePath(), ex);
			System.exit(1);
			return;
		}
        System.out.println("Risk score calculation complete");
        LOGGER.info("Risk score calculation complete");

    }
// get sample dosages
    // 0,1,2   ....    1.5   1.4   1.33434   
    // imputation not sure .... in between aa ab bb
    
    
    private static void startLogging(File logFile, boolean debugMode) {

        try {
            FileAppender logAppender = new FileAppender(new SimpleLayout(), logFile.getCanonicalPath(), false);
            Logger.getRootLogger().removeAllAppenders();
            Logger.getRootLogger().addAppender(logAppender);
            if (debugMode) {
                LOGGER.setLevel(Level.DEBUG);
            } else {
                LOGGER.setLevel(Level.INFO);
            }
        } catch (IOException e) {
            System.err.println("Failed to create logger: " + e.getMessage());
            System.exit(1);
        }

        LOGGER.info(
                "\n" + HEADER);
        LOGGER.info("Version: " + VERSION);
        LOGGER.info("Current date and time: " + DATE_TIME_FORMAT.format(currentDataTime));
        LOGGER.info("Log level: " + LOGGER.getLevel());

        System.out.println("Started logging");
        System.out.println();
    }

}


//        index = 0;
//        for (float mydosage : dosages) {
//            if (index < 5) {
//                System.out.println("dosage:" + mydosage); // only dosages
//            }
//            index++;
//        }
//
//        for (String myid : myids) {
//            System.out.println("ID : " + myid);// one rs id  
//        }
        // print allele info etc
//        if (debug == 1) {
//            String onesample = inputGenotypes.getSampleNames()[0];
//            System.out.println("First sample: " + onesample);
//
//            System.out.println("MINORALLELE  " + snpVariantByPos.getMinorAllele() + "  ALLALLELES  " + snpVariantByPos.getVariantAlleles() + "  MAF  "
//                    + snpVariantByPos.getMinorAlleleFrequency());
//
//            for (String mySequence : inputGenotypes.getSeqNames()) {
//                System.out.println("Sequencename: " + mySequence);   // chromosomes
//            }
//
//            index = 0;
//            for (Alleles allele : alleles) {
//                if (index < 20) {
//                    System.out.println(allele.toString());  // only alleles
//                    System.out.println(allele.getAllelesAsString().get(0));
//                    System.out.println(allele.getAllelesAsString().get(1));
//                    System.out.println(allele.isAtOrGcSnp());
//                    System.out.println("-------------------------------------");
//                }
//                index++;
//            }
//
//        }
        //System.out.println("Total amount of samples: " + inputGenotypes.getSampleNames().length);
        //System.out.println("Total amount of alleles: " + alleles.size());