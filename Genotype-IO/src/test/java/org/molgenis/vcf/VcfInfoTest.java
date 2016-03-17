package org.molgenis.vcf;

import static org.mockito.Mockito.mock;
import static org.mockito.Mockito.when;
import static org.testng.Assert.assertEquals;
import static org.testng.Assert.assertNull;

import java.util.Arrays;

import org.molgenis.vcf.meta.VcfMeta;
import org.molgenis.vcf.meta.VcfMetaInfo;
import org.molgenis.vcf.meta.VcfMetaInfo.Type;
import org.testng.annotations.Test;

public class VcfInfoTest
{

	@Test(expectedExceptions = IllegalArgumentException.class)
	public void VcfInfoVcfMeta()
	{
		new VcfInfo(null);
	}

	@Test(expectedExceptions = IllegalArgumentException.class)
	public void VcfInfoVcfMetaStringString()
	{
		new VcfInfo(null, null, null);
	}

	@Test
	public void getKey()
	{
		assertEquals(new VcfInfo(mock(VcfMeta.class) , "key", "val").getKey(), "key");
	}

	@Test
	public void getVal_noMetaInfo()
	{
		VcfInfo vcfInfo = new VcfInfo(mock(VcfMeta.class), "key", "treat-as-string");
		assertEquals(vcfInfo.getVal(), "treat-as-string");
	}
	
	@Test
	public void getVal_Character()
	{
		VcfMeta vcfMeta = mock(VcfMeta.class);
		VcfMetaInfo metaInfo = mock(VcfMetaInfo.class);
		when(metaInfo.getType()).thenReturn(Type.CHARACTER);
		when(metaInfo.getNumber()).thenReturn("1");
		when(vcfMeta.getInfoMeta("key")).thenReturn(metaInfo);
		VcfInfo vcfInfo = new VcfInfo(vcfMeta, "key", "c");
		assertEquals(vcfInfo.getVal(), Character.valueOf('c'));
	}

	@Test
	public void getVal_CharacterList()
	{
		VcfMeta vcfMeta = mock(VcfMeta.class);
		VcfMetaInfo metaInfo = mock(VcfMetaInfo.class);
		when(metaInfo.getType()).thenReturn(Type.CHARACTER);
		when(metaInfo.getNumber()).thenReturn("3");
		when(vcfMeta.getInfoMeta("key")).thenReturn(metaInfo);
		VcfInfo vcfInfo = new VcfInfo(vcfMeta, "key", "a,b,c");
		assertEquals(vcfInfo.getVal(), Arrays.asList('a', 'b', 'c'));
	}
	
	@Test
	public void getVal_Integer()
	{
		VcfMeta vcfMeta = mock(VcfMeta.class);
		VcfMetaInfo metaInfo = mock(VcfMetaInfo.class);
		when(metaInfo.getType()).thenReturn(Type.INTEGER);
		when(metaInfo.getNumber()).thenReturn("1");
		when(vcfMeta.getInfoMeta("key")).thenReturn(metaInfo);
		VcfInfo vcfInfo = new VcfInfo(vcfMeta, "key", "1");
		assertEquals(vcfInfo.getVal(), Integer.valueOf(1));
	}
	
	@Test
	public void getVal_IntegerList()
	{
		VcfMeta vcfMeta = mock(VcfMeta.class);
		VcfMetaInfo metaInfo = mock(VcfMetaInfo.class);
		when(metaInfo.getType()).thenReturn(Type.INTEGER);
		when(metaInfo.getNumber()).thenReturn("3");
		when(vcfMeta.getInfoMeta("key")).thenReturn(metaInfo);
		VcfInfo vcfInfo = new VcfInfo(vcfMeta, "key", "1,2,3,.");
		assertEquals(vcfInfo.getVal(), Arrays.asList(1, 2, 3, null));
	}
	
	@Test
	public void getVal_Flag()
	{
		VcfMeta vcfMeta = mock(VcfMeta.class);
		VcfMetaInfo metaInfo = mock(VcfMetaInfo.class);
		when(metaInfo.getType()).thenReturn(Type.FLAG);
		when(metaInfo.getNumber()).thenReturn("0");
		VcfInfo vcfInfo = new VcfInfo(vcfMeta, "key", null);
		assertNull(vcfInfo.getVal());
	}
	
	@Test
	public void getVal_Float()
	{
		VcfMeta vcfMeta = mock(VcfMeta.class);
		VcfMetaInfo metaInfo = mock(VcfMetaInfo.class);
		when(metaInfo.getType()).thenReturn(Type.FLOAT);
		when(metaInfo.getNumber()).thenReturn("1");
		when(vcfMeta.getInfoMeta("key")).thenReturn(metaInfo);
		VcfInfo vcfInfo = new VcfInfo(vcfMeta, "key", "1.23");
		assertEquals(vcfInfo.getVal(), Float.valueOf(1.23f));
	}
	
	@Test
	public void getVal_FloatList()
	{
		VcfMeta vcfMeta = mock(VcfMeta.class);
		VcfMetaInfo metaInfo = mock(VcfMetaInfo.class);
		when(metaInfo.getType()).thenReturn(Type.FLOAT);
		when(metaInfo.getNumber()).thenReturn("4");
		when(vcfMeta.getInfoMeta("key")).thenReturn(metaInfo);
		VcfInfo vcfInfo = new VcfInfo(vcfMeta, "key", "1.23,4.56,7.89,.,na,NaN,Na");
		assertEquals(vcfInfo.getVal(), Arrays.asList(Float.valueOf(1.23f), Float.valueOf(4.56f), Float.valueOf(7.89f), null, Float.NaN, Float.NaN, Float.NaN));
	}
	
	@Test
	public void getVal_String()
	{
		VcfMeta vcfMeta = mock(VcfMeta.class);
		VcfMetaInfo metaInfo = mock(VcfMetaInfo.class);
		when(metaInfo.getType()).thenReturn(Type.STRING);
		when(metaInfo.getNumber()).thenReturn("1");
		when(vcfMeta.getInfoMeta("key")).thenReturn(metaInfo);
		VcfInfo vcfInfo = new VcfInfo(vcfMeta, "key", "aaa");
		assertEquals(vcfInfo.getVal(), "aaa");
	}
	
	@Test
	public void getVal_StringList()
	{
		VcfMeta vcfMeta = mock(VcfMeta.class);
		VcfMetaInfo metaInfo = mock(VcfMetaInfo.class);
		when(metaInfo.getType()).thenReturn(Type.STRING);
		when(metaInfo.getNumber()).thenReturn("3");
		when(vcfMeta.getInfoMeta("key")).thenReturn(metaInfo);
		VcfInfo vcfInfo = new VcfInfo(vcfMeta, "key", "aaa,bbb,ccc");
		assertEquals(vcfInfo.getVal(), Arrays.asList("aaa", "bbb", "ccc"));
	}
	
	@Test
	public void reset()
	{
		VcfInfo vcfInfo = new VcfInfo(mock(VcfMeta.class), "key1", "val1");
		assertEquals(vcfInfo.getKey(), "key1");
		assertEquals(vcfInfo.getVal(), "val1");
		vcfInfo.reset("key2", "val2");
		assertEquals(vcfInfo.getKey(), "key2");
		assertEquals(vcfInfo.getVal(), "val2");
	}
}
