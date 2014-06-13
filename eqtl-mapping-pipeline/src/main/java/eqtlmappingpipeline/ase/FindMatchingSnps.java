package eqtlmappingpipeline.ase;

import java.util.ArrayList;
import java.util.TreeMap;
import org.molgenis.genotype.RandomAccessGenotypeData;
import org.molgenis.genotype.RandomAccessGenotypeDataReaderFormats;
import org.molgenis.genotype.variant.GeneticVariant;

/**
 *
 * @author Patrick Deelen
 */
public class FindMatchingSnps {

	/**
	 * @param args the command line arguments
	 */
	public static void main(String[] args) throws Exception {

		String referenceVcfFilePath = "";
		String queryVcfFilePath = "";

		RandomAccessGenotypeData referenceGenotypeData = RandomAccessGenotypeDataReaderFormats.VCF.createGenotypeData(referenceVcfFilePath);

		TreeMap<Double, ArrayList<ChrPos>> mafMap = new TreeMap<Double, ArrayList<ChrPos>>();
	
		for (GeneticVariant variant : referenceGenotypeData) {
			
//			mafMap.put(variant.getMinorAlleleFrequency(), new ChrPos(variant.getSequenceName(), variant.getStartPos()));
			
		}
		
		for (GeneticVariant variant : referenceGenotypeData) {
			
			
			
		}
		



	}

	private static class ChrPos {

		private final String chr;
		private final int pos;

		public ChrPos(String chr, int pos) {
			this.chr = chr;
			this.pos = pos;
		}

		public int getPos() {
			return pos;
		}

		public String getChr() {
			return chr;
		}
	}
}
