package nl.systemsgenetics.depict2;

import nl.systemsgenetics.depict2.pathway.PathwayEnrichments;
import umcg.genetica.math.matrix2.DoubleMatrixDataset;

import java.util.List;

public class Depict2Step2Results {


    private DoubleMatrixDataset<String, String> genePvalues;
    private List<PathwayEnrichments> pathwayEnrichments;


    public Depict2Step2Results(List<PathwayEnrichments> pathwayEnrichments, DoubleMatrixDataset<String, String> genePvalues) {
        this.genePvalues = genePvalues;
        this.pathwayEnrichments = pathwayEnrichments;
    }

    public DoubleMatrixDataset<String, String> getGenePvalues() {
        return genePvalues;
    }

    public void setGenePvalues(DoubleMatrixDataset<String, String> genePvalues) {
        this.genePvalues = genePvalues;
    }

    public List<PathwayEnrichments> getPathwayEnrichments() {
        return pathwayEnrichments;
    }

    public void setPathwayEnrichments(List<PathwayEnrichments> pathwayEnrichments) {
        this.pathwayEnrichments = pathwayEnrichments;
    }
}
