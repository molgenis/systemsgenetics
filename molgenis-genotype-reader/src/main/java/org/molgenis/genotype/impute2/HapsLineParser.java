package org.molgenis.genotype.impute2;

import java.util.ArrayList;
import java.util.List;
import java.util.StringTokenizer;

/**
 * Parses a line in a haps file.
 * 
 * This assumes the file is tab separated !!! Official haps file is space separated.
 * 
 * We use tabix indexed files, wich must be tab separated
 * 
 * 
 * @author erwin
 * 
 */
public class HapsLineParser
{
	public static HapsEntry parse(String line)
	{
		StringTokenizer tokenizer = new StringTokenizer(line, "	");
		int chrom = Integer.parseInt(tokenizer.nextToken());
		String snpId = tokenizer.nextToken();
		int position = Integer.parseInt(tokenizer.nextToken());
		String firstAllele = tokenizer.nextToken();
		String secondAllele = tokenizer.nextToken();

		List<String[]> sampleAlleles = new ArrayList<String[]>();
		while (tokenizer.hasMoreTokens())
		{
			String[] sampleAllele = new String[]
			{ tokenizer.nextToken(), tokenizer.nextToken() };
			sampleAlleles.add(sampleAllele);
		}

		return new HapsEntry(chrom, snpId, position, firstAllele, secondAllele, sampleAlleles);
	}
}
