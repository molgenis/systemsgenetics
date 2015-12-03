package nl.systemsgenetics.geneticriskscorecalculator;

import java.io.File;
import java.io.IOException;
import java.util.Arrays;
import java.util.List;
import umcg.genetica.math.matrix2.DoubleMatrixDataset;

/**
 *
 * @author Patrick Deelen
 */
public class RiskScoreMatrix {

    private final DoubleMatrixDataset<String, String> riskScores;

    public RiskScoreMatrix(String[] samples, List<String> phenotypes) {
        this(Arrays.asList(samples), phenotypes);
    }

    public RiskScoreMatrix(List<String> samples, List<String> phenotypes) {

        riskScores = new DoubleMatrixDataset<String, String>(samples, phenotypes);

    }

    public void setRiskScore(String sample, String phenotype, double score) {
        // System.out.println("MATRIX" + score + "MATRIX");
        riskScores.setElement(sample, phenotype, score);
    }

    public void save(File file) throws IOException {
        riskScores.save(file);
    }

    public double getRiskScore(String sample, String phenotype) {
        return riskScores.getElement(sample, phenotype);
    }

    public int rows() {
        return riskScores.rows();
    }

    public int cols() {
        return riskScores.columns();
    }
}
