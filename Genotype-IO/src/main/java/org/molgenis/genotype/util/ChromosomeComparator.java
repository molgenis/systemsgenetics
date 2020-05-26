package org.molgenis.genotype.util;

import java.io.Serializable;
import java.util.Comparator;

public class ChromosomeComparator implements Comparator<String>, Serializable
{
	
	public static final ChromosomeComparator comparator = new ChromosomeComparator();
	
	public static boolean chrASmallerChrB(String chrA, String chrB){
		return comparator.compare(chrA, chrB) < 0;
	}
	
	public static boolean chrASmallerEqualChrB(String chrA, String chrB){
		return comparator.compare(chrA, chrB) <= 0;
	}
	
	public static boolean chrALargerChrB(String chrA, String chrB){
		return comparator.compare(chrA, chrB) > 0;
	}

	@Override
	public int compare(String arg0, String arg1)
	{

		if (arg0.equals(arg1))
		{
			return 0;
		}

		if (arg0.equals("MT"))
		{
			return 1;
		}

		if (arg1.equals("MT"))
		{
			return -1;
		}

		if (arg0.length() == 1)
		{
			if (arg1.length() == 1)
			{
				return arg0.charAt(0) - arg1.charAt(0);
			}
			else if (arg0.charAt(0) == 'X')
			{
				return 1;
			}
			else if (arg0.charAt(0) == 'Y')
			{
				return 1;
			}
			else
			{
				return -1;
			}
		}

		if (arg1.length() == 1)
		{
			if (arg1.charAt(0) == 'X')
			{
				return -1;
			}
			else if (arg1.charAt(0) == 'Y')
			{
				return -1;
			}
			else
			{
				return 1;
			}
		}

		return arg0.compareTo(arg1);

	}
}
