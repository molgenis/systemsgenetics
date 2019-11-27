/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package nl.systemsgenetics.depict2;

import java.io.File;
import java.util.HashSet;
import org.apache.log4j.Logger;
import umcg.genetica.math.matrix2.DoubleMatrixDataset;

/**
 *
 * @author patri
 */
public class GwasSummStats {

	private static final Logger LOGGER = Logger.getLogger(GwasSummStats.class);

	private final String trait;
	private final File summStatsFile;
	private final String pvalueColumn;

	public GwasSummStats(String trait, String summStatsPath, String pvalueColumn) {
		this.trait = trait;
		this.summStatsFile = new File(summStatsPath);
		this.pvalueColumn = pvalueColumn;
	}

	public String getTrait() {
		return trait;
	}

	public File getSummStatsFile() {
		return summStatsFile;
	}

	public String getPvalueColumn() {
		return pvalueColumn;
	}

	public DoubleMatrixDataset<String, String> loadSubsetSummStats(HashSet<String> variantsToLoad) throws Exception {
		HashSet<String> pvalueColumnSet = new HashSet<>();
		pvalueColumnSet.add(pvalueColumn);

		DoubleMatrixDataset<String, String> summStats = DoubleMatrixDataset.loadSubsetOfTextDoubleData(summStatsFile.getAbsolutePath(), '\t', variantsToLoad, pvalueColumnSet);

		LOGGER.debug("Loaded summ stats for " + trait + " from: " + summStatsFile.getAbsolutePath() + ". Variants: " + summStats.rows());

		if (summStats.columns() != 1) {
			throw new RuntimeException("Summary statistics files read here should only contain on p-value column");
		}

		summStats.getHashCols().remove(pvalueColumn);
		summStats.getHashCols().put(trait, 0);

		return summStats;

	}

	public DoubleMatrixDataset<String, String> loadSummStats() throws Exception {
		return loadSubsetSummStats(null);
	}

}
