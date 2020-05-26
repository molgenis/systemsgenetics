package org.molgenis.genotype.util;


import java.util.Arrays;
import java.util.List;
import static org.testng.Assert.assertEquals;

import org.testng.annotations.Test;

public class UtilsTest {

	@Test
	public void iteratorToList() {
		List<String> list = Arrays.asList("1", "2", null);
		assertEquals(Utils.iteratorToList(list.iterator()), list);
	}
	
}

