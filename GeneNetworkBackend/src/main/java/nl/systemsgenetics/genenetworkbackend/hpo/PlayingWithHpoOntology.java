/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package nl.systemsgenetics.genenetworkbackend.hpo;

import umcg.genetica.io.hpo.HpoOntology;
import java.io.File;
import java.util.List;
import org.biojava.nbio.ontology.Term;

/**
 *
 * @author patri
 */
public class PlayingWithHpoOntology {

	/**
	 * @param args the command line arguments
	 */
	public static void main(String[] args) throws Exception {
		
		final File hpoOboFile = new File("C:\\UMCG\\Genetica\\Projects\\GeneNetwork\\HPO\\135\\hp.obo");
		

		HpoOntology hpoOntology = new HpoOntology(hpoOboFile);


		List<Term> mainPhenotypes = hpoOntology.getChildern(hpoOntology.nameToTerm("HP:0000118"));
		
		for(Term mainPhenotype : mainPhenotypes){
			System.out.println(mainPhenotype.getName());
		}
		
		
	}
	
}
