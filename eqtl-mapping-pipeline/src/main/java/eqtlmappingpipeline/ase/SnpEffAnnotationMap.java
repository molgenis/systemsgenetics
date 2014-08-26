package eqtlmappingpipeline.ase;

import java.io.BufferedReader;
import java.io.FileInputStream;
import java.io.IOException;
import java.io.InputStreamReader;
import java.util.List;
import org.apache.commons.lang3.StringUtils;
import org.molgenis.genotype.snpeff.SnpEffAnnotationTranslator;
import org.molgenis.genotype.snpeff.SnpEffEffect;
import umcg.genetica.collections.ChrPosMap;

/**
 *
 * @author Patrick Deelen
 */
public class SnpEffAnnotationMap {


	public static ChrPosMap<SnpEffEffect[]> loadSnpEffAnnotationMap(String vcfPath) throws IOException, Exception {
		
		ChrPosMap<SnpEffEffect[]> effects = new ChrPosMap<SnpEffEffect[]>();

		BufferedReader reader = new BufferedReader(new InputStreamReader(new FileInputStream(vcfPath), "UTF-8"));

		String line;
		while ((line = reader.readLine()) != null) {
			if(line.charAt(0) == '#'){
				continue;
			}
			String[] elements = StringUtils.split(line, '\t');
			SnpEffEffect[] effect = SnpEffAnnotationTranslator.translateSnpEffVcfField(new String(elements[7]));
			effects.put(new String(elements[0]), Integer.parseInt(elements[1]), effect);
		}
		
		return effects;
		
	}
	
	
}
