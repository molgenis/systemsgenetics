package umcg.genetica.io.binInteraction;

import gnu.trove.map.hash.TObjectIntHashMap;
import java.util.HashSet;
import umcg.genetica.io.binInteraction.gene.BinaryInteractionGeneCreator;
import umcg.genetica.io.binInteraction.variant.BinaryInteractionVariantCreator;

/**
 *
 * @author Patrick Deelen
 */
public class BinaryInteractionFileCreator {

	private final BinaryInteractionVariantCreator[] variants;
	private final BinaryInteractionGeneCreator[] genes;
	private final BinaryInteractionCohort[] cohorts;
	private final String[] covariats;
	private final boolean allCovariants;
	private final boolean metaAnalysis;
	private final boolean normalQtlStored;
	private final boolean flippedZscoreStored;
	private String description = "";
	private int[][] covariatesTested;
	private long interactions = 0;
	private long startQtlBlock;
	private long startInteractionBlock;
	
	private final TObjectIntHashMap<String> variantMap;
	private final TObjectIntHashMap<String> genesMap;

	public BinaryInteractionFileCreator(BinaryInteractionVariantCreator[] variants, BinaryInteractionGeneCreator[] genes, BinaryInteractionCohort[] cohorts, String[] covariats, boolean allCovariants, boolean metaAnalysis, boolean normalQtlStored, boolean flippedZscoreStored) throws BinaryInteractionFileException {
		this.variants = variants;
		this.genes = genes;
		this.cohorts = cohorts;
		this.covariats = covariats;
		this.allCovariants = allCovariants;
		this.metaAnalysis = metaAnalysis;
		this.normalQtlStored = normalQtlStored;
		this.flippedZscoreStored = flippedZscoreStored;
		
		variantMap = new TObjectIntHashMap<String>(variants.length, 0.75f, -1);
		genesMap = new TObjectIntHashMap<String>(genes.length, 0.75f, -1);
		
		for(int i = 0 ; i < variants.length ; ++i){
			if(variantMap.put(variants[i].getName(), i) == -1){
				throw new BinaryInteractionFileException("Cannot store the same variant twice (" + variants[i].getName() + ")");
			}
		}
		
		for(int i = 0 ; i < genes.length ; ++i){
			if(genesMap.put(genes[i].getName(), i) == -1){
				throw new BinaryInteractionFileException("Cannot store the same gene twice (" + genes[i].getName() + ")");
			}
		}
		
		HashSet<String> tmp = new HashSet<String>(covariats.length);
		for(int i = 0 ; i < covariats.length ; ++i){
			if(!tmp.add(covariats[i])){
				throw new BinaryInteractionFileException("Cannot store the same covariate twice (" + covariats[i] + ")");
			}
		}
		
		tmp.clear();
		for(int i = 0 ; i < cohorts.length ; ++i){
			if(!tmp.add(cohorts[i].getName())){
				throw new BinaryInteractionFileException("Cannot store the same cohort twice (" + cohorts[i] + ")");
			}
		}
		
	}

	public void addTestedVariantGene(String variant, String gene) throws BinaryInteractionFileException{
		
		//TODO
		
	}
	
	public void addTestedInteraction(String variant, String gene, String covariate) throws BinaryInteractionFileException{
		
		
		
		//TODO
		
	}

	public void setDescription(String description) {
		this.description = description;
	}
	
	
	
}
