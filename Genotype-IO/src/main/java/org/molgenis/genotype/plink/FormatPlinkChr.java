package org.molgenis.genotype.plink;

import java.util.regex.Matcher;
import java.util.regex.Pattern;

/**
 *
 * @author Patrick Deelen
 */
public class FormatPlinkChr {

	private static final Pattern CHR_PATTERN = Pattern.compile("^chr(.*)$", Pattern.CASE_INSENSITIVE);
	
	public static String formatChr(String chrName){
		
		Matcher chrMatcher = CHR_PATTERN.matcher(chrName);
		if (chrMatcher.find()) {
			chrName = chrMatcher.group(1);
		}
		
		switch(chrName){
			case "X": return "23";
			case "Y": return "24";
			case "XY": return "25";
			case "MT": return "26";
			default: return chrName;
		}
		
	}
	
}
