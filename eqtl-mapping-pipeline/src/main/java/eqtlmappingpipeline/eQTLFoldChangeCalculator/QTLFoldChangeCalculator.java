/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package eqtlmappingpipeline.eQTLFoldChangeCalculator;

import eqtlmappingpipeline.metaqtl3.MetaQTL3;
import umcg.genetica.io.text.TextFile;
import umcg.genetica.io.trityper.SNP;
import umcg.genetica.io.trityper.SNPLoader;
import umcg.genetica.io.trityper.util.BaseAnnot;

/**
 *
 * @author harmjan
 */
public class QTLFoldChangeCalculator extends MetaQTL3 {

    public void run() {
    }

    void calculateFoldChanges(String eqtlfile) throws Exception {


	for (int d = 0; d < m_gg.length; d++) {
	    TextFile tf = new TextFile(eqtlfile, TextFile.R);
	    tf.readLineElemsReturnReference(TextFile.tab);
	    String[] elems = tf.readLineElemsReturnReference(TextFile.tab);

	    SNPLoader loader = m_gg[d].getGenotypeData().createSNPLoader();
	    int[] indWGA = m_gg[d].getExpressionToGenotypeIdArray();
	    int numgg = m_gg[d].getTotalGGSamples();
	    while (elems != null) {
		String snp = elems[1];
		String probe = elems[4];
		String assessedAllele = elems[9];
		Integer snpId = m_gg[d].getGenotypeData().getSnpToSNPId().get(snp);
		Integer probeId = m_gg[d].getExpressionData().getProbeToId().get(probe);

		if (snpId != -9 && probeId != -9) {
		    SNP snpObject = m_gg[d].getGenotypeData().getSNPObject(snpId);
		    loader.loadGenotypes(snpObject);
		    double[] expression = m_gg[d].getExpressionData().getMatrix()[probeId];
		    byte[] genotypes = snpObject.getGenotypes();
		    double[] genotypemeans = new double[3];
		    int numAA = 0;
		    int numBB = 0;

		    double sumAA = 0;
		    double sumBB = 0;

		    for (int i = 0; i < numgg; i++) {
			if (indWGA[i] != -1) {
			    if (genotypes[indWGA[i]] == 0) {
				sumAA += (expression[i] * 2);
				numAA += 2;
			    }
			    if (genotypes[indWGA[i]] == 1) {
				sumAA += (expression[i]);
				sumBB += (expression[i]);
				numAA++;
				numBB++;
			    }
			    if (genotypes[indWGA[i]] == 2) {
				sumBB += (expression[i] * 2);
				numBB += 2;
			    }
			}

		    }





		    sumAA /= (double) numAA;
		    sumBB /= (double) numBB;
		    double foldchange = 0;
		    if (assessedAllele.equals(BaseAnnot.toString(snpObject.getAlleles()[0]))) {
			foldchange = sumAA / sumBB;
			System.out.println(snp + "\t" + probe + "\t" + BaseAnnot.toString(snpObject.getAlleles()[0]) + "/" + BaseAnnot.toString(snpObject.getAlleles()[1]) + "\t" + foldchange);
		    } else {
			foldchange = sumBB / sumAA;
			System.out.println(snp + "\t" + probe + "\t" + BaseAnnot.toString(snpObject.getAlleles()[1]) + "/" + BaseAnnot.toString(snpObject.getAlleles()[0]) + "\t" + foldchange);
		    }

		}



		elems = tf.readLineElemsReturnReference(TextFile.tab);
	    }
	    tf.close();
	}

    }
}
