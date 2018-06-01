package org.molgenis.genotype;

import java.io.File;
import java.net.URISyntaxException;

public class ResourceTest
{
	public File getTestResourceFile(String name) throws URISyntaxException
	{
		return new File(this.getClass().getResource(name).toURI());
	}
    
    protected File getTriTyperFolder() throws URISyntaxException
	{
		return getTestResourceFile("/TriTyper");
	}

	protected File getTestVcfGz() throws URISyntaxException
	{
		return getTestResourceFile("/test.vcf.gz");
	}

		protected File getTestVcfGz2() throws URISyntaxException
	{
		return getTestResourceFile("/chr21.imputed.head.vcf.gz");
	}

	protected File getTestVcf1() throws URISyntaxException
	{
		return getTestResourceFile("/test1.vcf");
	}

	protected File getTestVcfGzTbi() throws URISyntaxException
	{
		return getTestResourceFile("/test.vcf.gz.tbi");
	}
	
		protected File getTestVcfGz2Tbi() throws URISyntaxException
	{
		return getTestResourceFile("/chr21.imputed.head.vcf.gz.tbi");
	}

	protected File getTestMap() throws URISyntaxException
	{
		return getTestResourceFile("/test.map");
	}

	protected File getTestPed() throws URISyntaxException
	{
		return getTestResourceFile("/test.ped");
	}

	protected File getTestBed9() throws URISyntaxException
	{
		return getTestResourceFile("/test.bed");
	}

	protected File getTestBim9() throws URISyntaxException
	{
		return getTestResourceFile("/test.bim");
	}

	protected File getTestFam9() throws URISyntaxException
	{
		return getTestResourceFile("/test.fam");
	}
	
	protected File getTestBed6() throws URISyntaxException
	{
		return getTestResourceFile("/test6.bed");
	}

	protected File getTestBim6() throws URISyntaxException
	{
		return getTestResourceFile("/test6.bim");
	}

	protected File getTestFam6() throws URISyntaxException
	{
		return getTestResourceFile("/test6.fam");
	}
	
	protected File getTestBed7() throws URISyntaxException
	{
		return getTestResourceFile("/test7.bed");
	}

	protected File getTestBim7() throws URISyntaxException
	{
		return getTestResourceFile("/test7.bim");
	}

	protected File getTestFam7() throws URISyntaxException
	{
		return getTestResourceFile("/test7.fam");
	}
	
	protected File getTestBed8() throws URISyntaxException
	{
		return getTestResourceFile("/test8.bed");
	}

	protected File getTestBim8() throws URISyntaxException
	{
		return getTestResourceFile("/test8.bim");
	}

	protected File getTestFam8() throws URISyntaxException
	{
		return getTestResourceFile("/test8.fam");
	}

	protected File getTestImpute2Gz() throws URISyntaxException
	{
		return getTestResourceFile("/test.haps.tab.gz");
	}

	protected File getTestImpute2GzTbi() throws URISyntaxException
	{
		return getTestResourceFile("/test.haps.tab.gz.tbi");
	}

	protected File getTestImpute2Sample() throws URISyntaxException
	{
		return getTestResourceFile("/test.sample");
	}

	protected File getTestImpute2Haps() throws URISyntaxException
	{
		return getTestResourceFile("/test.haps");
	}

	protected File getLdTestVcf() throws URISyntaxException
	{
		return getTestResourceFile("/ldTest.vcf.gz");
	}

	protected File getLdTestVcfTbi() throws URISyntaxException
	{
		return getTestResourceFile("/ldTest.vcf.gz.tbi");
	}
	
	protected File getTest2Gen() throws URISyntaxException
	{
		return getTestResourceFile("/test2.gen");
	}
	
	protected File getTest2Sample() throws URISyntaxException
	{
		return getTestResourceFile("/test2.sample");
	}
	
	protected File getTest3Bgen() throws URISyntaxException
	{
		return getTestResourceFile("/test3.bgen");
	}
	
	protected File getTest3Sample() throws URISyntaxException
	{
		return getTestResourceFile("/test3.sample");
	}
    
    protected File getTest10Gen() throws URISyntaxException
	{
		return getTestResourceFile("/test10.gen");
	}
}
