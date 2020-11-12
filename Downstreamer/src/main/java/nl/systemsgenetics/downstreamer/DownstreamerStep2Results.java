package nl.systemsgenetics.downstreamer;

import nl.systemsgenetics.downstreamer.pathway.PathwayEnrichments;
import umcg.genetica.math.matrix2.DoubleMatrixDataset;

import java.util.List;

public class DownstreamerStep2Results {


    private DoubleMatrixDataset<String, String> genePvalues;
    private DoubleMatrixDataset<String, String> normalizedGenePvalues;
    private List<PathwayEnrichments> pathwayEnrichments;

    public DownstreamerStep2Results(List<PathwayEnrichments> pathwayEnrichments, DoubleMatrixDataset<String, String> genePvalues) {
        this.genePvalues = genePvalues;
        this.pathwayEnrichments = pathwayEnrichments;
    }

    public DownstreamerStep2Results(List<PathwayEnrichments> pathwayEnrichments, DoubleMatrixDataset<String, String> genePvalues, DoubleMatrixDataset<String, String> normalizedGenePvalues) {
        this.genePvalues = genePvalues;
        this.normalizedGenePvalues = normalizedGenePvalues;
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

    public DoubleMatrixDataset<String, String> getNormalizedGenePvalues() {
        return normalizedGenePvalues;
    }
}
