package nl.systemsgenetics.geneticriskscorecalculator;

import gnu.trove.map.hash.TObjectDoubleHashMap;
import org.molgenis.genotype.RandomAccessGenotypeData;

/**
 *
 * @author Patrick Deelen
 */
public class SimpleGeneticRiskScoreCalculator implements GeneticRiskScoreCalculator{

	@Override
	public TObjectDoubleHashMap<String> calculateRiskScores(RandomAccessGenotypeData genotypeData) {
		throw new UnsupportedOperationException("Not supported yet.");
	}

	@Override
	public TObjectDoubleHashMap<String> calculateRiskScores(RandomAccessGenotypeData genotypeData, PhenotypeData phenotypeData) {
		throw new UnsupportedOperationException("Not supported yet.");
	}

	@Override
	public String getPhenotype() {
		throw new UnsupportedOperationException("Not supported yet.");
	}

}
