package nl.systemsgenetics.downstreamer;


import java.util.HashMap;
import java.util.List;
import java.util.Map;
import nl.systemsgenetics.downstreamer.containers.GwasLocus;

public class DownstreamerStep3Results {

    private Map<String, List<GwasLocus>> lociPerTrait;

    public DownstreamerStep3Results() {
        this.lociPerTrait = new HashMap<>();
    }

    public DownstreamerStep3Results(Map<String, List<GwasLocus>> loci) {
        this.lociPerTrait = loci;
    }

    public void addLoci(String trait, List<GwasLocus> loci) {
        lociPerTrait.put(trait, loci);
    }

    public Map<String, List<GwasLocus>> getLoci() {
        return lociPerTrait;
    }

    public void setLoci(Map<String, List<GwasLocus>> loci) {
        this.lociPerTrait = loci;
    }
}
