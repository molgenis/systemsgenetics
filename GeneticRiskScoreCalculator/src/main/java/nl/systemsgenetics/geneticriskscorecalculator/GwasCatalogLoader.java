package nl.systemsgenetics.geneticriskscorecalculator;

import java.util.ArrayList;
import java.util.List;
import umcg.genetica.io.text.TextFile;
import java.util.regex.Pattern;
import java.io.File;
import java.io.IOException;

/**
 *
 * @author Patrick Deelen
 */
public class GwasCatalogLoader {

    private static final Pattern TAB_PATTERN = Pattern.compile("\\t");

    public List<GeneticRiskScoreCalculator> getGeneticRiskScoreCalculators(String myFile) {

        List<GeneticRiskScoreCalculator> geneticRiskScoreCalculators = new ArrayList<GeneticRiskScoreCalculator>();

        List<String> phenotypes = new ArrayList<String>();
        //phenotypes.add("Height");
        phenotypes.add("CD");
        phenotypes.add("UC");
        phenotypes.add("IBD");
        //phenotypes.add("Breast cancer");
        //phenotypes.add("Ulcerative colitis");
        //phenotypes.add("Crohns disease");
        //phenotypes.add("Obesity");
       // phenotypes.add("Body mass index");
        

        
         
        for (String currentPhenotype : phenotypes) {

            SimpleGeneticRiskScoreCalculator firstCalculator = new SimpleGeneticRiskScoreCalculator();

            firstCalculator.setPhenotype(currentPhenotype);

            String fileLine;
            String[] fileLineData;
            TextFile riskSnpsFile;

            int index = 0;
            try {
                riskSnpsFile = new TextFile(myFile, false);
                riskSnpsFile.readLine(); // reading the header
                while ((fileLine = riskSnpsFile.readLine()) != null) {
                    fileLineData = TAB_PATTERN.split(fileLine);
                    if (index < 5) {
                        // System.out.println("--" + fileLineData[1] + "--" + fileLineData[2] + "--" + fileLineData[3]);
                    }
                    if (index < 5) {
                        // System.out.println("RiskAllele: " + fileLineData[4]);
                    }

                    try {

                        if (fileLineData[0].equals(currentPhenotype)) {
                            
                            
                            firstCalculator.addChr(Integer.valueOf(fileLineData[1]));  // like this or via setter?
                            firstCalculator.addPos(Integer.valueOf(fileLineData[2]));
                            firstCalculator.addRsid(fileLineData[3]);
                            firstCalculator.addRiskallele(fileLineData[4]);
                            firstCalculator.addOtherallele(fileLineData[5]);
                            firstCalculator.addPvalue(Double.valueOf(fileLineData[7]));
                        }

                    } catch (NumberFormatException e) {
                        System.err.println("no integer in chr or pos field");

                    }
                    
                    // try to read in OR or beta separately
                    try {

                        if (fileLineData[0].equals(currentPhenotype)) {
                            firstCalculator.addOrorbeta(Double.valueOf(fileLineData[8]));
                        }

                    } catch (NumberFormatException e) {
                        firstCalculator.addOrorbeta(Double.valueOf("0"));
                        System.err.println("effect size missing: set to zero");

                    }                   
                    
                    
                    
                    index++;
                }
                riskSnpsFile.close();
            } catch (IOException ex) {
                System.err.println("Unable to load risk snps file.");
                //   LOGGER.fatal("Unable to load risk snps file.", ex);
                System.exit(1);
                //  return;
            }

            geneticRiskScoreCalculators.add(firstCalculator);

        }
        return geneticRiskScoreCalculators;
    }

}
