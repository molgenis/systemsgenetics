/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */

package umcg.genetica.io.trityper;

import gnu.trove.set.hash.THashSet;

/**
 * @author harmjan
 */
public class TriTyperGeneticalGenomicsDatasetSettings {
	public String name;
	public String genotypeLocation;
	public String snpFileLocation;
	public String snpmapFileLocation;
	public String expressionLocation;
	public String genotypeToExpressionCoupling;
	public boolean logtransform;
	public boolean quantilenormalize;
	public THashSet<String> tsProbesConfine = null;
	public boolean confineProbesToProbesMappingToAnyChromosome = false;
	public Integer confineProbesToProbesThatMapToChromosome = null;
	public String expressionplatform;
	public String probeannotation;
	public boolean cisAnalysis, transAnalysis;
	public String covariateFile;

	@Override
	public String toString() {
		return "TriTyperGeneticalGenomicsDatasetSettings{" +
				"name='" + name + '\'' +
				"\n, genotypeLocation='" + genotypeLocation + '\'' +
				"\n, snpFileLocation='" + snpFileLocation + '\'' +
				"\n, snpmapFileLocation='" + snpmapFileLocation + '\'' +
				"\n, expressionLocation='" + expressionLocation + '\'' +
				"\n, genotypeToExpressionCoupling='" + genotypeToExpressionCoupling + '\'' +
				"\n, logtransform=" + logtransform +
				"\n, quantilenormalize=" + quantilenormalize +
				"\n, tsProbesConfine=" + tsProbesConfine +
				"\n, confineProbesToProbesMappingToAnyChromosome=" + confineProbesToProbesMappingToAnyChromosome +
				"\n, confineProbesToProbesThatMapToChromosome=" + confineProbesToProbesThatMapToChromosome +
				"\n, expressionplatform='" + expressionplatform + '\'' +
				"\n, probeannotation='" + probeannotation + '\'' +
				"\n, cisAnalysis=" + cisAnalysis +
				"\n, transAnalysis=" + transAnalysis +
				"\n, covariateFile='" + covariateFile + '\'' +
				'}';
	}
}
