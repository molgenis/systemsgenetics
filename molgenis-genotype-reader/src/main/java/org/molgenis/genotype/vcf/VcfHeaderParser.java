package org.molgenis.genotype.vcf;

import org.molgenis.util.tuple.KeyValueTuple;
import org.molgenis.util.tuple.Tuple;

/**
 * Parser for vcf headers.
 * 
 * A VCF header can be key=value, for example ##fileformat=VCFv4.1
 * 
 * Or can be key=tuple, for example ##INFO=<ID=NS,Number=1,Type
 * =Integer,Description="Number of Samples With Data"> The key value pairs are
 * inside < >
 * 
 * Question can '<' be escaped? 
 * For now we assume not, so if '<' is encountered we assume a complex header.
 * 
 * @author erwin
 * 
 */
public class VcfHeaderParser
{
	private String rawHeader;

	public VcfHeaderParser(String header)
	{
		if (header == null) throw new IllegalArgumentException("Header can not be null");
		if (!header.startsWith("##")) throw new IllegalArgumentException("A VCF header must start with '##'");
		if (!header.contains("=")) throw new IllegalArgumentException(
				"A VCF header must at least contain one '=', is key value based");
		if (header.contains("<") && !header.contains(">")) throw new IllegalArgumentException("Missing '>'");
		if (header.contains(">") && !header.contains("<")) throw new IllegalArgumentException("Missing '<'");

		rawHeader = header;
	}

	/**
	 * Parse the raw text header The tuple returned can contain one key and one
	 * stringvalue in case of a simple header or can contain a key and a tuple
	 * as value in case of complex headers like ##INFO
	 * 
	 * 
	 * @return the tuple containing the header
	 */
	public Tuple parse()
	{
		KeyValueTuple result = new KeyValueTuple();

		// header block starts with < and ends with >
		boolean inHeaderBlock = false;
		// values can be quoted to allow for = and ,
		boolean inQuotes = false;
		// header is divided in key and value using '='
		boolean inKey = true;
		// to store a key while parsing
		String key = "";
		// to store a value while parsing
		String value = "";

		for (Character c : rawHeader.toCharArray())
		{
			if (!inHeaderBlock)
			{
				if ('<' == c)
				{
					inHeaderBlock = true;
				}
			}
			else
			{
				// check for seperator between key and value
				if (inKey && '=' == c)
				{
					inKey = false;
				}
				// parse quotes
				else if (!inKey && '"' == c)
				{
					if (inQuotes)
					{
						inQuotes = false;
						result.set(key, value);
						key = "";
						value = "";
						inKey = true;
					}
					else
						inQuotes = true;
				}
				// close the key/value pair
				else if (!inQuotes && (',' == c || '>' == c))
				{
					result.set(key, value);
					key = "";
					value = "";
					inKey = true;
				}
				// otherwise just add key/value char
				else
				{
					if (inKey) key += c;
					else
						value += c;
				}
			}
		}

		return result;
	}

}
