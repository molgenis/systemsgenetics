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

        SimpleGeneticRiskScoreCalculator firstCalculator = new SimpleGeneticRiskScoreCalculator();

        firstCalculator.setPhenotype("height");

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
                    System.out.println("--" + fileLineData[1] + "--" + fileLineData[2] + "--" + fileLineData[3]);
                }
                if (index < 5) {
                    System.out.println("RiskAllele: " + fileLineData[4]);
                }
                firstCalculator.chr.add(Integer.valueOf(fileLineData[1]));  // like this or via setter?
                firstCalculator.pos.add(Integer.valueOf(fileLineData[2]));
                firstCalculator.rsid.add(fileLineData[3]);
                firstCalculator.riskallele.add(fileLineData[4]);
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

        return geneticRiskScoreCalculators;
    }

}
