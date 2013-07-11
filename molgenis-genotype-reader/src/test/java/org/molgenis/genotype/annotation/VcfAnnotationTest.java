package org.molgenis.genotype.annotation;

import static org.testng.Assert.assertEquals;
import static org.testng.Assert.assertFalse;
import static org.testng.Assert.assertNotNull;
import static org.testng.Assert.assertNull;
import static org.testng.Assert.assertTrue;

import org.molgenis.genotype.vcf.VcfInfo;
import org.testng.annotations.Test;

public class VcfAnnotationTest
{
	@Test
	public void fromVcfInfo()
	{
		VcfInfo info = new VcfInfo("ID", VcfInfo.Type.CHARACTER, "A", "Desc");
		VcfAnnotation annotation = VcfAnnotation.fromVcfInfo(info);
		assertNotNull(annotation);
		assertEquals(annotation.getId(), info.getId());
		assertEquals(annotation.getName(), info.getId());
		assertEquals(annotation.getDescription(), info.getDescription());
		assertNull(annotation.getNumber());
		assertEquals(annotation.getType(), Annotation.Type.CHAR);
		assertTrue(annotation.isPerAltAllele());
		assertFalse(annotation.isPerGenotype());
		assertTrue(annotation.isList());
		assertFalse(annotation.isUnbounded());

		info = new VcfInfo("ID", VcfInfo.Type.FLAG, "G", "Desc");
		annotation = VcfAnnotation.fromVcfInfo(info);
		assertNotNull(annotation);
		assertEquals(annotation.getId(), info.getId());
		assertEquals(annotation.getName(), info.getId());
		assertEquals(annotation.getDescription(), info.getDescription());
		assertNull(annotation.getNumber());
		assertEquals(annotation.getType(), Annotation.Type.BOOLEAN);
		assertFalse(annotation.isPerAltAllele());
		assertTrue(annotation.isPerGenotype());
		assertTrue(annotation.isList());
		assertFalse(annotation.isUnbounded());

		info = new VcfInfo("ID", VcfInfo.Type.STRING, ".", "Desc");
		annotation = VcfAnnotation.fromVcfInfo(info);
		assertNotNull(annotation);
		assertEquals(annotation.getId(), info.getId());
		assertEquals(annotation.getName(), info.getId());
		assertEquals(annotation.getDescription(), info.getDescription());
		assertNull(annotation.getNumber());
		assertEquals(annotation.getType(), Annotation.Type.STRING);
		assertFalse(annotation.isPerAltAllele());
		assertFalse(annotation.isPerGenotype());
		assertTrue(annotation.isList());
		assertTrue(annotation.isUnbounded());

		info = new VcfInfo("ID", VcfInfo.Type.UNKNOWN, "1", "Desc");
		annotation = VcfAnnotation.fromVcfInfo(info);
		assertNotNull(annotation);
		assertEquals(annotation.getId(), info.getId());
		assertEquals(annotation.getName(), info.getId());
		assertEquals(annotation.getDescription(), info.getDescription());
		assertNotNull(annotation.getNumber());
		assertEquals(annotation.getNumber().intValue(), 1);
		assertEquals(annotation.getType(), Annotation.Type.UNKOWN);
		assertFalse(annotation.isPerAltAllele());
		assertFalse(annotation.isPerGenotype());
		assertFalse(annotation.isList());
		assertFalse(annotation.isUnbounded());
	}
}
