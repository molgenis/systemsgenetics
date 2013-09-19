package org.molgenis.genotype.plink.drivers;

import java.io.File;
import java.net.URI;
import java.net.URISyntaxException;

public abstract class AbstractResourceTest
{
	protected File getTestResource(String name) throws URISyntaxException
	{
		URI resource = this.getClass().getResource(name).toURI();
		return new File(resource);
	}
}
