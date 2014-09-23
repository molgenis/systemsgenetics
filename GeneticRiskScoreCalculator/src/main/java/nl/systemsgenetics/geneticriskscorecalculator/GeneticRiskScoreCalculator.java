package nl.systemsgenetics.geneticriskscorecalculator;

import gnu.trove.map.hash.TObjectDoubleHashMap;
import org.molgenis.genotype.RandomAccessGenotypeData;

/**
 *
 * @author Patrick Deelen
 */
public interface GeneticRiskScoreCalculator {

	TObjectDoubleHashMap<String> calculateRiskScores(RandomAccessGenotypeData genotypeData);

	TObjectDoubleHashMap<String> calculateRiskScores(RandomAccessGenotypeData genotypeData, PhenotypeData phenotypeData);
	
	String getPhenotype();
        
	
}
