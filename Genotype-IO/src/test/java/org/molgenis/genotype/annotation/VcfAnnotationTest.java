package org.molgenis.genotype.annotation;

import static org.testng.Assert.assertEquals;
import static org.testng.Assert.assertFalse;
import static org.testng.Assert.assertNotNull;
import static org.testng.Assert.assertNull;
import static org.testng.Assert.assertTrue;

import java.util.HashMap;
import java.util.Map;

import org.molgenis.vcf.meta.VcfMetaInfo;
import org.testng.annotations.Test;

public class VcfAnnotationTest
{
	@Test
	public void fromVcfInfo_CharacterList()
	{
		Map<String, String> properties = new HashMap<String, String>();
		properties.put("ID", "ID");
		properties.put("Type", "Character");
		properties.put("Number", "A");
		properties.put("Description", "Desc");
		VcfMetaInfo info = new VcfMetaInfo(properties);
		
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
	}
	
	@Test
	public void fromVcfInfo_FlagList()
	{
		Map<String, String> properties = new HashMap<String, String>();
		properties.put("ID", "ID");
		properties.put("Type", "Flag");
		properties.put("Number", "G");
		properties.put("Description", "Desc");
		VcfMetaInfo info = new VcfMetaInfo(properties);
		
		VcfAnnotation annotation = VcfAnnotation.fromVcfInfo(info);
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
	}
	
	@Test
	public void fromVcfInfo_StringList()
	{
		Map<String, String> properties = new HashMap<String, String>();
		properties.put("ID", "ID");
		properties.put("Type", "String");
		properties.put("Number", ".");
		properties.put("Description", "Desc");
		VcfMetaInfo info = new VcfMetaInfo(properties);
		
		VcfAnnotation annotation = VcfAnnotation.fromVcfInfo(info);
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
	}
	
	@Test
	public void fromVcfInfo_Integer()
	{
		Map<String, String> properties = new HashMap<String, String>();
		properties.put("ID", "ID");
		properties.put("Type", "Integer");
		properties.put("Number", "1");
		properties.put("Description", "Desc");
		VcfMetaInfo info = new VcfMetaInfo(properties);
		
		VcfAnnotation annotation = VcfAnnotation.fromVcfInfo(info);
		assertNotNull(annotation);
		assertEquals(annotation.getId(), info.getId());
		assertEquals(annotation.getName(), info.getId());
		assertEquals(annotation.getDescription(), info.getDescription());
		assertEquals(annotation.getNumber(), Integer.valueOf(info.getNumber()));
		assertEquals(annotation.getType(), Annotation.Type.INTEGER);
		assertFalse(annotation.isPerAltAllele());
		assertFalse(annotation.isPerGenotype());
		assertFalse(annotation.isList());
		assertFalse(annotation.isUnbounded());
	}
}
