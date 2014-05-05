package org.molgenis.vcf.meta;

import java.util.Map;

public class VcfMetaFormat extends VcfMetaEntry
{
	private static final String KEY_ID = "ID";
	private static final String KEY_NUMBER = "Number";
	private static final String KEY_TYPE = "Type";
	private static final String KEY_DESCRIPTION = "Description";
	
	public enum Type {
		INTEGER("Integer"), FLOAT("Float"), CHARACTER("Character"), STRING("String");
		
		private String type;
		
		private Type(String type) {
			this.type = type;
		}
		
		public static Type from(String str) {
			for(Type type : values()) {
				if(type.toString().equals(str))
					return type;
			}
			return null;
		}
		
		@Override
		public String toString() {
			return type;
		}
	}

	public VcfMetaFormat(Map<String, String> properties) {
		super(properties);
	}

	@Override
	public String getName()
	{
		return "FORMAT";
	}
	
	public String getId()
	{
		return properties.get(KEY_ID);
	}

	public String getNumber()
	{
		return properties.get(KEY_NUMBER);
	}

	public Type getType()
	{
		return Type.from(properties.get(KEY_TYPE));
	}

	public String getDescription()
	{
		return properties.get(KEY_DESCRIPTION);
	}
}
