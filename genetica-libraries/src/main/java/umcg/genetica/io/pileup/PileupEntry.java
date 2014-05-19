package umcg.genetica.io.pileup;

import gnu.trove.map.hash.TObjectDoubleHashMap;
import gnu.trove.map.hash.TObjectIntHashMap;
import org.molgenis.genotype.Allele;

/**
 *
 * @author Patrick Deelen
 */
public class PileupEntry {

	private final String chr;
	private final int pos;
	private final Allele refAllele;
	private final int readDepth;
	private final TObjectIntHashMap<Allele> alleleCounts;
	private final TObjectDoubleHashMap<Allele> alleleAverageQualities;
	private final int minimumBaseQuality;

	/**
	 * Only handles SNPs
	 *
	 * @param chr
	 * @param pos
	 * @param refAllele
	 * @param readDepth
	 * @param basesString
	 */
	public PileupEntry(String chr, int pos, Allele refAllele, int readDepth, String basesString, int minimumBaseQuality) throws PileupParseException {
		this(chr, pos, refAllele, readDepth, basesString, null, minimumBaseQuality);
	}

	public PileupEntry(String chr, int pos, Allele refAllele, int readDepth, String basesString, String basesQualityString, int minimumBaseQuality) throws PileupParseException {
		this.chr = chr;
		this.pos = pos;
		this.refAllele = refAllele;
		this.readDepth = readDepth;
		this.alleleCounts = new TObjectIntHashMap<Allele>();
		this.alleleAverageQualities = new TObjectDoubleHashMap<Allele>();
		this.minimumBaseQuality = minimumBaseQuality;

		alleleCounts.put(Allele.A, 0);
		alleleCounts.put(Allele.C, 0);
		alleleCounts.put(Allele.G, 0);
		alleleCounts.put(Allele.T, 0);

		alleleAverageQualities.put(Allele.A, 0);
		alleleAverageQualities.put(Allele.C, 0);
		alleleAverageQualities.put(Allele.G, 0);
		alleleAverageQualities.put(Allele.T, 0);

		if (!alleleCounts.containsKey(refAllele)) {
			throw new PileupParseException("Error parsing pipeup entry");
		}

		int[] basesQuality = basesQualityString == null ? null : parseBasesQualityString(basesQualityString);
		parseBasesString(basesString, basesQuality);
	}

	/**
	 * Parse bases string and put in allele counts
	 *
	 * @param basesString
	 */
	private void parseBasesString(final String basesString, final int[] basesQualities) throws PileupParseException {

		char[] basesChars = basesString.toCharArray();

		int basesQualityI = 0;
		for (int i = 0; i < basesChars.length; ++i) {

			switch (basesChars[i]) {
				case '.':
				case ',':
					if(basesQualities == null || basesQualities[basesQualityI] >= minimumBaseQuality){
						alleleCounts.increment(refAllele);
						if (basesQualities != null) {
							alleleAverageQualities.adjustValue(refAllele, basesQualities[basesQualityI]);
						}
					}
					++basesQualityI;
					break;
				case '-':
				case '+':
					//insertion and deletion. First determine size then skip
					StringBuilder indelSizeString = new StringBuilder();
					++i;
					while (Character.isDigit(basesChars[i])) {
						indelSizeString.append(basesChars[i]);
						++i;
					}
					//indel size -1 since for loop will do ++i
					i += Integer.valueOf(indelSizeString.toString()) - 1;
					break;
				case '^':
					++i; //skip read mapping quality
					break;
				case '$':
				case '>':
				case '<':
				case 'N':
				case 'n':
				case '*':
					break;
				case 'a':
				case 'c':
				case 't':
				case 'g':
					basesChars[i] = Character.toUpperCase(basesChars[i]);
				//No break rest is same as capital chars
				case 'A':
				case 'C':
				case 'T':
				case 'G':
					if(basesQualities == null || basesQualities[basesQualityI] >= minimumBaseQuality){
						Allele allele = Allele.create(basesChars[i]);
						alleleCounts.increment(allele);
						if (basesQualities != null) {
							alleleAverageQualities.adjustValue(allele, basesQualities[basesQualityI]);
						}
					}
					++basesQualityI;
					break;
				default:
					throw new PileupParseException("Unexpected char in pileup bases string: " + basesChars[i] + " from : " + basesString);

			}

		}

		for (Allele allele : alleleAverageQualities.keySet()) {
			if (basesQualities == null) {
				alleleAverageQualities.put(allele, Double.NaN);
			} else {
				alleleAverageQualities.put(allele, alleleAverageQualities.get(allele) / alleleCounts.get(allele));
			}
		}

	}

	public String getChr() {
		return chr;
	}

	public int getPos() {
		return pos;
	}

	public Allele getRefAllele() {
		return refAllele;
	}

	public int getReadDepth() {
		return readDepth;
	}

	public TObjectIntHashMap<Allele> getAlleleCounts() {
		return alleleCounts;
	}

	public int getAlleleCount(Allele allele) {
		return alleleCounts.get(allele);
	}

	public TObjectDoubleHashMap<Allele> getAlleleAverageQualities() {
		return alleleAverageQualities;
	}
	
	public double getAlleleAverageQuality(Allele allele) {
		return alleleAverageQualities.get(allele);
	}

	public int getMinimumBaseQuality() {
		return minimumBaseQuality;
	}
	
	private int[] parseBasesQualityString(String basesQualityString) {

		char[] basesQualityChars = basesQualityString.toCharArray();

		int[] basesQuality = new int[basesQualityChars.length];

		for (int i = 0; i < basesQualityChars.length; ++i) {
			basesQuality[i] = ((int) basesQualityChars[i]) - 33;
		}

		return basesQuality;

	}
}
