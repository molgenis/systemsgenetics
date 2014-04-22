package org.molgenis.vcf.meta;

import java.util.Map;

public class VcfMetaInfo extends VcfMetaEntry
{
	private static final String KEY_ID = "ID";
	private static final String KEY_NUMBER = "Number";
	private static final String KEY_TYPE = "Type";
	private static final String KEY_DESCRIPTION = "Description";
	private static final String KEY_SOURCE = "Source";
	private static final String KEY_VERSION = "Version";
	
	public enum Type {
		INTEGER("Integer"), FLOAT("Float"), FLAG("Flag"), CHARACTER("Character"), STRING("String");
		
		private final String type;
		
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

	public VcfMetaInfo(Map<String, String> properties) {
		super(properties);
	}

	@Override
	public String getName()
	{
		return "INFO";
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

	public String getSource()
	{
		return properties.get(KEY_SOURCE);
	}

	public String getVersion()
	{
		return properties.get(KEY_VERSION);
	}
}
