/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package PathwayEnrichmentRestuls;

import java.io.File;
import java.util.HashSet;
import java.util.Set;
import nl.systemsgenetics.downstreamer.DownstreamerDeprecated;
import nl.systemsgenetics.downstreamer.pathway.PathwayDatabase;
import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import umcg.genetica.math.matrix2.DoubleMatrixDataset;

/**
 *
 * @author patri
 */
public class PathwayEnrichmentResults {

	private static final Logger LOGGER = LogManager.getLogger(DownstreamerDeprecated.class);

	private final PathwayDatabase pathwayDatabase;
	private final HashSet<String> hlaGenesToExclude;
	private final boolean ignoreGeneCorrelations;
	private DoubleMatrixDataset<String, String> betas;
	private DoubleMatrixDataset<String, String> pValues;
	private DoubleMatrixDataset<String, String> qValues;

	//private DoubleMatrixDataset<String, String> standardErrors;
	//private DoubleMatrixDataset<String, String> zscores;
	//private final DoubleMatrixDataset<String, String> pValuesNull;
	//private final DoubleMatrixDataset<String, String> betasNull;
	private final int numberOfPathways;
	private final File intermediateFolder;
	private final Set<String> excludeGenes;
	//private final List<Gene> genes;
	private final String outputBasePath;

	public PathwayEnrichmentResults(PathwayDatabase pathwayDatabase, DoubleMatrixDataset<String, String> betas, DoubleMatrixDataset<String, String> pValues, DoubleMatrixDataset<String, String> qValues) {
		this.pathwayDatabase = pathwayDatabase;
		this.hlaGenesToExclude = null;
		this.ignoreGeneCorrelations = false;
		this.betas = betas;
		this.pValues = pValues;
		this.qValues = qValues;
		this.numberOfPathways = pValues.rows();
		this.intermediateFolder = null;
		this.excludeGenes = null;
		this.outputBasePath = null;
	}
}
