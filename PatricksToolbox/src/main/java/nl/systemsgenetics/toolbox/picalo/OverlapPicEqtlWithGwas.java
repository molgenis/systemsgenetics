/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package nl.systemsgenetics.toolbox.picalo;

import java.io.IOException;
import java.util.Set;
import nl.systemsgenetics.toolbox.Options;
import nl.systemsgenetics.toolbox.Utils;
import umcg.genetica.io.gwascatalog.GWASCatalog;

/**
 *
 * @author patri
 */
public class OverlapPicEqtlWithGwas {
	
	public static void overlap(Options options)  throws IOException, Exception {
		
		GWASCatalog gwasCatalog = new GWASCatalog(options.getGwasCatalogFile(), 5e-8);
		
		Set<String> gwasVariants = gwasCatalog.getAllSnpsIds();
		
		Utils.loadGenotypes(options, gwasVariants);
		
		
	}
	
	
}
