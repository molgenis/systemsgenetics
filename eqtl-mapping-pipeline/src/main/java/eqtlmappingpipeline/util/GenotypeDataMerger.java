/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package eqtlmappingpipeline.util;

import java.io.IOException;
import umcg.genetica.io.trityper.util.TriTyperGenotypeDataMerger;

/**
 *
 * @author harmjan
 */
public class GenotypeDataMerger {

    public void merge(String in, String in2, String out, String snps) throws IOException{

	TriTyperGenotypeDataMerger merger = new TriTyperGenotypeDataMerger();
	merger.combinePrioritizerDatasetsMergeCommonSNPs(in, in2, out, snps);

    }
}
