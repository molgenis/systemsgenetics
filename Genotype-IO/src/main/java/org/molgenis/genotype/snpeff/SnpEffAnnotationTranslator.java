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

	public static SnpEffEffect[] translateSnpEffVcfField(String field) throws Exception {

		field = field.subSequence(4, field.length()).toString();
		
		String[] effectStrings = StringUtils.split(field, ',');

		SnpEffEffect[] effects = new SnpEffEffect[effectStrings.length];

		int i = 0;
		for (String effectString : effectStrings) {
			effects[i] = translateSnpEffEffectString(effectString);
			++i;
		}

		return effects;

	}

	private static SnpEffEffect translateSnpEffEffectString(String effectString) throws Exception {

		String[] effectFields = StringUtils.splitPreserveAllTokens(effectString, "(|)");

		if (effectFields.length < 11) {
			throw new Exception("Error parsing SnpEff annotation field count: " + effectFields.length + "\n\t" + effectString);
		}

		SnpEffEffect.EffectType effectType;
		try {
			effectType = SnpEffEffect.EffectType.valueOf(effectFields[EFFECT_TYPE_FIELD]);
		} catch (IllegalArgumentException ex) {
			System.err.println("Error parsing " + effectFields[EFFECT_TYPE_FIELD] + " as EffectType, Field: " + effectString);
			throw ex;
		}
		SnpEffEffect.EffectImpact effectImpact;
		try {
			effectImpact = SnpEffEffect.EffectImpact.valueOf(effectFields[EFFECT_INPACT_FIELD]);
		} catch (IllegalArgumentException ex) {
			System.err.println("Error parsing " + effectFields[EFFECT_INPACT_FIELD] + " as EffectImpact, Field: " + effectString);
			throw ex;
		}
		SnpEffEffect.FunctionalClass effectFunctionalClass;
		if (effectFields[FUNCTIONAL_CLASS_FIELD].length() == 0) {
			effectFunctionalClass = SnpEffEffect.FunctionalClass.NONE;
		} else {
			try {
				effectFunctionalClass = SnpEffEffect.FunctionalClass.valueOf(effectFields[FUNCTIONAL_CLASS_FIELD]);
			} catch (IllegalArgumentException ex) {
				System.err.println("Error parsing " + effectFields[FUNCTIONAL_CLASS_FIELD] + " as FunctionalClass, Field: " + effectString);
				throw ex;
			}
		}

		return new SnpEffEffect(effectType, effectImpact, effectFunctionalClass);


	}
}
