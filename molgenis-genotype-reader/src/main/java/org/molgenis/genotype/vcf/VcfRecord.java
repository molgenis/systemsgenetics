package org.molgenis.genotype.vcf;

import static org.molgenis.genotype.vcf.VcfUtils.checkNullValue;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;
import java.util.LinkedHashMap;
import java.util.List;
import java.util.Map;
import java.util.regex.Pattern;

import org.apache.commons.lang3.builder.ToStringBuilder;
import org.molgenis.util.tuple.KeyValueTuple;
import org.molgenis.util.tuple.Tuple;
import org.molgenis.util.tuple.WritableTuple;

/**
 * Class that represents one row in a VCF file
 * 
 * @author erwin
 * 
 */
public class VcfRecord
{
	public static final String GENOTYPE_FORMAT = "GT";
	private final Tuple record;
	private final static Pattern TAB_PATTERN = Pattern.compile("\\t");
	private final static Pattern COMMA_PATTERN = Pattern.compile(",");

	// cache for the info map
	private Map<String, List<String>> infoMap = null;

	public VcfRecord(Tuple record)
	{
		this.record = record;
	}

	public VcfRecord(String line, List<String> columnNames)
	{
		String[] values = TAB_PATTERN.split(line);
		if (values.length != columnNames.size())
		{
			throw new IllegalArgumentException("The number of columns does not match the number of columnnames");
		}

		record = new KeyValueTuple();

		for (int i = 0; i < values.length; i++)
		{
			//Create new string that is not backed by whole line. Otherwise we keep to much of the genotype data in memory
			((WritableTuple) record).set(columnNames.get(i), new String(values[i]));
		}

	}

	public String getChrom()
	{
		return record.getString("#CHROM");
	}

	public Integer getPos()
	{
		return record.getInt("POS");
	}

	public List<String> getId()
	{
		String id = checkNullValue(record.getString("ID"));
		if (id == null)
		{
			return Collections.emptyList();
		}

		return Collections.unmodifiableList(Arrays.asList(id.split(";")));
	}

	public String getRef()
	{
		return record.getString("REF");
	}

	public List<String> getAlt()
	{
		String alt = checkNullValue(record.getString("ALT"));
		if (alt == null)
		{
			return Collections.emptyList();
		}

		return Collections.unmodifiableList(Arrays.asList(COMMA_PATTERN.split(alt)));
	}

	/**
	 * Get all posible alleles, the first in the list is the reference
	 * 
	 * @return
	 */
	public List<String> getAlleles()
	{
		List<String> alleles = new ArrayList<String>();
		alleles.add(getRef());
		alleles.addAll(getAlt());

		return alleles;
	}

	public Double getQual()
	{
		String qual = checkNullValue(record.getString("QUAL"));
		if (qual == null)
		{
			return null;
		}

		return record.getDouble("QUAL");
	}

	public List<String> getFilter()
	{
		String filter = checkNullValue(record.getString("FILTER"));
		if (filter == null)
		{
			return null;
		}

		return Collections.unmodifiableList(Arrays.asList(filter.split(";")));
	}

	public String getInfo()
	{
		return checkNullValue(record.getString("INFO"));
	}

	public List<String> getFormat()
	{
		String format = checkNullValue(record.getString("FORMAT"));
		if (format == null)
		{
			return Collections.emptyList();
		}

		return Collections.unmodifiableList(Arrays.asList(format.split(":")));
	}

	public String getSampleValue(String sampleName, String key)
	{
		// first get the position from the key
		int index = getFormat().indexOf(key);

		if (index == -1)
		{
			return null;
		}

		// then get the sample from the tuple
		String sampleRecord = record.getString(sampleName);

		// and parse out the value
		if (sampleRecord != null)
		{
			String[] values = sampleRecord.split(":");
			if (index < values.length) return values[index];
		}

		return null;
	}

	/**
	 * See javadoc in VcfSampleGenotype.java
	 * 
	 * @param sampleName
	 * @return
	 */
	public VcfSampleGenotype getSampleGenotype(String sampleName)
	{
		String value = getSampleValue(sampleName, GENOTYPE_FORMAT);
		if (value == null)
		{
			return null;
		}

		return new VcfSampleGenotypeParser(value).parse();
	}

	public List<String> getInfo(String key)
	{
		if (infoMap == null)
		{
			infoMap = new LinkedHashMap<String, List<String>>();

			String infoString = getInfo();
			
			if(infoString != null){

				String[] keyvalues = infoString.split(";");
				for (String keyvalue : keyvalues)
				{
					String[] kv = keyvalue.split("=");
					if (kv.length == 1)
					{
						infoMap.put(kv[0], Arrays.asList(new String[]
						{ "TRUE" }));
					}
					else
					{
						if (VcfUtils.checkNullValue(kv[1]) != null)
						{
							infoMap.put(kv[0], Arrays.asList(COMMA_PATTERN.split(kv[1])));
						}
					}
				}
			}
		}

		if (!infoMap.keySet().contains(key))
		{
			return Collections.emptyList();
		}

		return Collections.unmodifiableList(infoMap.get(key));
	}

	@Override
	public String toString()
	{
		return ToStringBuilder.reflectionToString(this);
	}

}
