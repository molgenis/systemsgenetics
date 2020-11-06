package nl.systemsgenetics.downstreamer;

import nl.systemsgenetics.downstreamer.summarystatistic.Locus;

import java.util.HashMap;
import java.util.List;
import java.util.Map;

public class DownstreamerStep3Results {

    private Map<String, List<Locus>> lociPerTrait;

    public DownstreamerStep3Results() {
        this.lociPerTrait = new HashMap<>();
    }

    public DownstreamerStep3Results(Map<String, List<Locus>> loci) {
        this.lociPerTrait = loci;
    }

    public void addLoci(String trait, List<Locus> loci) {
        lociPerTrait.put(trait, loci);
    }

    public Map<String, List<Locus>> getLoci() {
        return lociPerTrait;
    }

    public void setLoci(Map<String, List<Locus>> loci) {
        this.lociPerTrait = loci;
    }
}
