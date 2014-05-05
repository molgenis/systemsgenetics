package eqtlmappingpipeline.ase.pileup;

import gnu.trove.map.hash.TObjectIntHashMap;
import java.util.ArrayList;
import org.molgenis.genotype.Allele;

/**
 *
 * @author Patrick Deelen
 */
public class pileupEntry {

	private final String chr;
	private final int pos;
	private final Allele refAllele;
	private final int readDepth;
	private final TObjectIntHashMap<Allele> alleleCounts;

	/**
	 * Only handle SNPs
	 *
	 * @param chr
	 * @param pos
	 * @param refAllele
	 * @param readDepth
	 * @param basesString
	 */
	public pileupEntry(String chr, int pos, Allele refAllele, int readDepth, String basesString) throws PileupParseException {
		this.chr = chr;
		this.pos = pos;
		this.refAllele = refAllele;
		this.readDepth = readDepth;
		this.alleleCounts = new TObjectIntHashMap<Allele>();

		alleleCounts.put(Allele.A, 0);
		alleleCounts.put(Allele.C, 0);
		alleleCounts.put(Allele.G, 0);
		alleleCounts.put(Allele.T, 0);

		if (!alleleCounts.containsKey(refAllele)) {
			throw new PileupParseException("Error parsing pipeup entry");
		}

		parseBasesString(basesString);

	}

	/**
	 * Parse bases string and put in allele counts
	 *
	 * @param basesString
	 */
	private void parseBasesString(String basesString) throws PileupParseException {

		char[] basesChars = basesString.toCharArray();

		for (int i = 0; i < basesChars.length; ++i) {

			switch (basesChars[i]) {
				case '.':
				case ',':
					alleleCounts.increment(refAllele);
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
					Allele allele = Allele.create(basesChars[i]);
					alleleCounts.increment(allele);
					break;
				default:
					throw new PileupParseException("Unexpected char in pileup bases string");

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
}
