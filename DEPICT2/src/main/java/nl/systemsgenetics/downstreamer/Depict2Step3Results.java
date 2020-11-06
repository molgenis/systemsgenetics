package nl.systemsgenetics.downstreamer;

import nl.systemsgenetics.downstreamer.summarystatistic.Locus;

import java.util.HashMap;
import java.util.List;
import java.util.Map;

public class Depict2Step3Results {

    private Map<String, List<Locus>> lociPerTrait;

    public Depict2Step3Results() {
        this.lociPerTrait = new HashMap<>();
    }

    public Depict2Step3Results(Map<String, List<Locus>> loci) {
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
