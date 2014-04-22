package org.molgenis.vcf.meta;

import java.util.Map;

public abstract class VcfMetaEntry
{
	protected final Map<String, String> properties;

	public VcfMetaEntry(Map<String, String> properties) {
		if(properties == null) throw new IllegalArgumentException("properties is null");
		this.properties = properties;
	}
	
	public abstract String getName();
	
	public String get(String key) {
		return properties.get(key);
	}
	
	@Override
	public String toString() {
		String prefix = getName() + '\t';
		StringBuilder strBuilder = new StringBuilder(prefix);
		for(Map.Entry<String, String> entry : properties.entrySet()) {
			if(strBuilder.length() > prefix.length()) strBuilder.append(' ');
			strBuilder.append(entry.getKey()).append('=').append('[').append(entry.getValue()).append(']');
		}
		strBuilder.append('\n');
		return strBuilder.toString();
	}

	@Override
	public int hashCode()
	{
		final int prime = 31;
		int result = 1;
		result = prime * result + ((properties == null) ? 0 : properties.hashCode());
		return result;
	}

	@Override
	public boolean equals(Object obj)
	{
		if (this == obj) return true;
		if (obj == null) return false;
		if (getClass() != obj.getClass()) return false;
		VcfMetaEntry other = (VcfMetaEntry) obj;
		if (properties == null)
		{
			if (other.properties != null) return false;
		}
		else if (!properties.equals(other.properties)) return false;
		return true;
	}
}
