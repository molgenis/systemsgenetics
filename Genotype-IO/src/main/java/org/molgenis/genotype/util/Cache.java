package org.molgenis.genotype.util;

import java.util.LinkedHashMap;
import java.util.Map;

/**
 * Least Recently Used (LRU) cache. The items accessed last is removed when
 * capacity is reached Code from:
 * http://java-planet.blogspot.nl/2005/08/how-to-set
 * -up-simple-lru-cache-using.html
 * 
 * @author Patrick Deelen
 * 
 */
public class Cache<K, V> extends LinkedHashMap<K, V>
{

	private static final long serialVersionUID = 1L;
	private final int capacity;

	public Cache(int capacity)
	{
		super(capacity + 1, 1.1f, true);
		this.capacity = capacity;
	}

	@Override
	synchronized public V get(Object key)
	{
		return super.get(key);
	}
	
	@Override
	synchronized public V put(K key, V value) {
		return super.put(key, value);
	}
	
	@Override
	protected boolean removeEldestEntry(Map.Entry<K, V> eldest)
	{
		return size() > capacity;
	}

}