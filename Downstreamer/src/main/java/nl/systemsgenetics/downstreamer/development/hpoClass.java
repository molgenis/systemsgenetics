/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package nl.systemsgenetics.downstreamer.development;

import java.io.File;
import java.io.FileNotFoundException;
import java.io.IOException;
import java.text.ParseException;
import java.util.ArrayList;
import java.util.List;
import org.biojava.nbio.ontology.Ontology;
import org.biojava.nbio.ontology.Term;
import umcg.genetica.io.hpo.HpoOntology;

/**
 *
 * @author patri
 */
public class hpoClass {

	/**
	 * @param args the command line arguments
	 */
	public static void main(String[] args) throws IOException, FileNotFoundException, ParseException {

		final File hpoOboFile = new File("D:\\UMCG\\Genetica\\Projects\\GeneNetwork\\HPO\\135\\hp.obo");

		HpoOntology hpo = new HpoOntology(hpoOboFile);
		
		Ontology ontology = hpo.getHpoOntologyData();

		Term phenotypicAbnormalityTerm = hpo.nameToTerm("HP:0000118");

		List<Term> phenotypicAbnormalityChilderen = hpo.getChildern(phenotypicAbnormalityTerm);

		for (Term phenotypicAbnormalityChild : phenotypicAbnormalityChilderen) {

			System.out.println(phenotypicAbnormalityChild.getName() + "\t" + phenotypicAbnormalityChild.getDescription());

		}
		
		ArrayList<hpoTermWithTopDiseaseClass> terms = new ArrayList();
		
		for (Term term : hpo.getAllTerms()){
			
			hpoTermWithTopDiseaseClass termInfo = new hpoTermWithTopDiseaseClass(term);
			for(Term topTerm : phenotypicAbnormalityChilderen){
				if(hpo.isTermABelowTermB(term, topTerm)){
					termInfo.addDiseaseClass(topTerm);
				}
			}
			terms.add(termInfo);
		}
		
		for(hpoTermWithTopDiseaseClass term : terms){
			
			if(term.getTerm().getDescription() != null && (term.getTerm().getDescription().equals("is_a-relationship") || term.getTerm().getDescription().equals("subset-relationship"))){
				continue;
			}
			
			System.out.println(term.getTerm().getName() + "\t" + term.getTerm().getDescription() + "\t" + term.getTopDiseaseCallesDescription());
			
		}
		
	}

	private static class hpoTermWithTopDiseaseClass {
		
		final private Term term;
		final private ArrayList<Term> topDiseaseClasses = new ArrayList<>();

		public hpoTermWithTopDiseaseClass(Term term) {
			this.term = term;
		}
		
		public void addDiseaseClass(Term diseaseClass){
			topDiseaseClasses.add(diseaseClass);
		}

		public Term getTerm() {
			return term;
		}

		public ArrayList<Term> getTopDiseaseClasses() {
			return topDiseaseClasses;
		}
		
		public String getTopDiseaseCallesDescription(){
			
			StringBuilder res = new StringBuilder();
			
			for(Term topTerm : topDiseaseClasses){
				if(res.length() > 0){
					res.append(';');
				}
				res.append(topTerm.getDescription());
			}
			return res.toString();
			
		}
		
	}
	
}
