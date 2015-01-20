package eqtlmappingpipeline.ase;

import cern.colt.list.tdouble.DoubleArrayList;
import cern.colt.list.tint.IntArrayList;
import java.util.ArrayList;
import org.molgenis.genotype.Allele;
import org.molgenis.genotype.variant.id.GeneticVariantId;

/**
 *
 * @author Patrick Deelen
 */
public class AseVariantBean implements AseVariant{

	private final String chr;
	private final int pos;
	private final GeneticVariantId id;
	private final Allele a1;
	private final Allele a2;
	private final IntArrayList a1Counts;
	private final IntArrayList a2Counts;
	private final DoubleArrayList pValues;
	private final ArrayList<String> sampleIds;
	private final double metaZscore;
	private final double metaPvalue;
	private final double countPearsonR;
    private double effect;
    private double LikelihoodRatioP;
    private double LikelihoodRatioD;
	private AseMleBeta mle;
	private final String genes;

	public AseVariantBean(String chr, int pos, GeneticVariantId id, Allele a1, Allele a2, IntArrayList a1Counts, IntArrayList a2Counts, DoubleArrayList pValues, ArrayList<String> sampleIds, double metaZscore, double metaPvalue, double countPearsonR) {
		this.chr = chr;
		this.pos = pos;
		this.id = id;
		this.a1 = a1;
		this.a2 = a2;
		this.a1Counts = a1Counts;
		this.a2Counts = a2Counts;
		this.pValues = pValues;
		this.sampleIds = sampleIds;
		this.metaZscore = metaZscore;
		this.metaPvalue = metaPvalue;
		this.countPearsonR = countPearsonR;
		genes = "";
	}
    
    public AseVariantBean(String outputLine[]) {
        this.effect = Double.parseDouble(outputLine[2]);
        this.LikelihoodRatioP = Double.parseDouble(outputLine[0]);
        this.LikelihoodRatioD = Double.parseDouble(outputLine[1]);
		this.chr = outputLine[5];
		this.pos = Integer.parseInt(outputLine[6]);
		this.id = GeneticVariantId.createVariantId(outputLine[7]);
		this.a1 = Allele.create(outputLine[9]);
		this.a2 = Allele.create(outputLine[10]);
		
		this.genes = outputLine[12];

        this.a1Counts = new IntArrayList();
        for(String s : outputLine[13].split(",")){
            a1Counts.add(Integer.parseInt(s));
        }
		
		this.a2Counts = new IntArrayList();
        for(String s : outputLine[14].split(",")){
            a2Counts.add(Integer.parseInt(s));
        }
        //outputLine[15]
		this.pValues = new DoubleArrayList();
        for(String s : outputLine[15].split(",")){
            pValues.add(Double.parseDouble(s));
        }
        //outputLine[16]
		this.sampleIds = new ArrayList<String>();
        for(String s : outputLine[15].split(",")){
            sampleIds.add(s);
        }

		this.metaZscore = Double.parseDouble(outputLine[4]);
		this.metaPvalue = Double.parseDouble(outputLine[3]);
		this.countPearsonR  = Double.parseDouble(outputLine[11]);
	}

	@Override
	public String getChr() {
		return chr;
	}

	@Override
	public int getPos() {
		return pos;
	}

	@Override
	public GeneticVariantId getId() {
		return id;
	}

	@Override
	public Allele getA1() {
		return a1;
	}

	@Override
	public Allele getA2() {
		return a2;
	}

	@Override
	public IntArrayList getA1Counts() {
		return a1Counts;
	}

	@Override
	public IntArrayList getA2Counts() {
		return a2Counts;
	}

	@Override
	public DoubleArrayList getPValues() {
		return pValues;
	}

	@Override
	public ArrayList<String> getSampleIds() {
		return sampleIds;
	}

	@Override
	public double getMetaZscore() {
		return metaZscore;
	}

	@Override
	public double getMetaPvalue() {
		return metaPvalue;
	}

	@Override
	public double getCountPearsonR() {
		return countPearsonR;
	}

	@Override
	public AseMleBeta getMle() {
		return mle;
	}

	@Override
	public int getSampleCount() {
		return a1Counts.size();
	}

	@Override
	public int compareTo(AseVariant o) {
		double thisRatioD = this.getMle().getRatioD();
		double otherRatioD = o.getMle().getRatioD();
		
		//Reverse compare. Largest first
		return Double.compare(otherRatioD, thisRatioD);
	}
	
    public DoubleArrayList getpValues() {
        return pValues;
    }

	@Override
    public double getEffect() {
        return effect;
    }

	@Override
    public double getLikelihoodRatioP() {
        return LikelihoodRatioP;
    }

	@Override
    public double getLikelihoodRatioD() {
        return LikelihoodRatioD;
    }

	@Override
	public void calculateStatistics() {
		//empty by design
	}

	public String getGenes() {
		return genes;
	}
	
	
}
