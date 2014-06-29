package org.molgenis.genotype.snpeff;

import java.util.ArrayList;
import java.util.List;
import org.apache.commons.lang3.StringUtils;

/**
 *
 * @author Patrick Deelen
 */
public final class SnpEffAnnotationTranslator {
	
	private static final int EFFECT_TYPE_FIELD = 0;
	private static final int EFFECT_INPACT_FIELD = 1;
	private static final int FUNCTIONAL_CLASS_FIELD = 2;

	public static List<SnpEffEffect> translateSnpEffVcfField(String field) throws Exception{
		
		String[] effectStrings = StringUtils.split(field, ',');
		
		ArrayList<SnpEffEffect> effects = new ArrayList<SnpEffEffect>(effectStrings.length);
		
		for(String effectString : effectStrings){
			effects.add(translateSnpEffEffectString(effectString));
		}
		
		return effects;
		
	}
	
	private static SnpEffEffect translateSnpEffEffectString(String effectString) throws Exception{
		
		String[] effectFields = StringUtils.split(effectString, "(|)");
		
		if(effectFields.length < 11){
			throw new Exception("Error parsing SnpEff annotation");
		}
		
		return new SnpEffEffect(SnpEffEffect.EffectType.valueOf(effectFields[EFFECT_TYPE_FIELD]), SnpEffEffect.EffectImpact.valueOf(effectFields[EFFECT_INPACT_FIELD]), SnpEffEffect.FunctionalClass.valueOf(effectFields[FUNCTIONAL_CLASS_FIELD]));
		
		
	}
	
}
