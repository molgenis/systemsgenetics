/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package org.molgenis.genotype.variant;

import java.util.Collections;
import java.util.HashMap;
import java.util.Map;
import org.molgenis.genotype.GenotypeDataException;

/**
 *
 * @author Patrick Deelen
 */
public class GeneticVariantMetaMap implements GeneticVariantMeta {

	private final Map<String, Type> metaMap;
	private static final Map<String, Type> RESEVERED_IDS;
	private static final GeneticVariantMetaMap GT_META_MAP = new GeneticVariantMetaMap(Collections.singletonMap("GT", Type.ALLELES));
	private static final GeneticVariantMetaMap GP_META_MAP = new GeneticVariantMetaMap(Collections.singletonMap("GP", Type.FLOAT_LIST));

	static {
		HashMap<String, Type> reservedIdsTmp = new HashMap<String, Type>();
		reservedIdsTmp.put("GT", Type.ALLELES);
		reservedIdsTmp.put("GP", Type.FLOAT_LIST);
		reservedIdsTmp.put("DS", Type.FLOAT);
		RESEVERED_IDS = Collections.unmodifiableMap(reservedIdsTmp);
	}

	private GeneticVariantMetaMap(Map<String, Type> metaMap) {
		this.metaMap = metaMap;
	}

	public static GeneticVariantMeta getGeneticVariantMetaGt() {
		return GT_META_MAP;
	}
	
	public static GeneticVariantMeta getGeneticVariantMetaGp() {
		return GP_META_MAP;
	}

	public static GeneticVariantMeta createGeneticVariantMeta(Map<String, Type> metaMap) {
		for (Map.Entry<String, Type> entry : metaMap.entrySet()) {
			if (checkNotOverwriteReserved(entry.getKey(), entry.getValue())) {
				throw new GenotypeDataException("Using illegal genotype field: " + entry.getKey() + " is reserved for: " + RESEVERED_IDS.get(entry.getKey()));
			}
		}
		return new GeneticVariantMetaMap(Collections.unmodifiableMap(metaMap));
	}

	private static boolean checkNotOverwriteReserved(String id, Type type) {
		if (RESEVERED_IDS.containsKey(id)) {
			if (RESEVERED_IDS.get(id).equals(type)) {
				return false;
			} else {
				return true;
			}
		} else {
			return false;
		}
	}

	@Override
	public Iterable<String> getRecordIds() {
		return metaMap.keySet();
	}

	@Override
	public Type getRecordType(String recordId) {
		return metaMap.get(recordId);
	}
}
