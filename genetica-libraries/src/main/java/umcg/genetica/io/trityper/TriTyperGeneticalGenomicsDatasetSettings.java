/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */

package umcg.genetica.io.trityper;

import java.util.HashSet;

/**
 *
 * @author harmjan
 */
public class TriTyperGeneticalGenomicsDatasetSettings {
    public String name;
    public String genotypeLocation;
    public String expressionLocation;
    public String genotypeToExpressionCoupling;
    public boolean logtransform;
    public boolean quantilenormalize;
    public HashSet<String> tsProbesConfine = null;
    public boolean confineProbesToProbesMappingToAnyChromosome = false;
    public Integer confineProbesToProbesThatMapToChromosome = null;
    public String expressionplatform;
    public String probeannotation;
    public boolean cisAnalysis, transAnalysis;
}
