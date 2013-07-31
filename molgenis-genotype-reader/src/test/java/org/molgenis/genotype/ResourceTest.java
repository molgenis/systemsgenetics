package org.molgenis.genotype;

import java.io.File;
import java.net.URISyntaxException;

public class ResourceTest
{
	protected File getTestResourceFile(String name) throws URISyntaxException
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

	protected File getTestVcf1() throws URISyntaxException
	{
		return getTestResourceFile("/test1.vcf");
	}

	protected File getTestVcfGzTbi() throws URISyntaxException
	{
		return getTestResourceFile("/test.vcf.gz.tbi");
	}

	protected File getTestMap() throws URISyntaxException
	{
		return getTestResourceFile("/test.map");
	}

	protected File getTestPed() throws URISyntaxException
	{
		return getTestResourceFile("/test.ped");
	}

	protected File getTestBed() throws URISyntaxException
	{
		return getTestResourceFile("/test.bed");
	}

	protected File getTestBim() throws URISyntaxException
	{
		return getTestResourceFile("/test.bim");
	}

	protected File getTestFam() throws URISyntaxException
	{
		return getTestResourceFile("/test.fam");
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

}
