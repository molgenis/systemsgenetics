/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package umcg.genetica.io.hpo;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileInputStream;
import java.io.FileNotFoundException;
import java.io.IOException;
import java.io.InputStreamReader;
import java.text.ParseException;
import java.util.ArrayList;
import java.util.Collections;
import java.util.List;
import java.util.Set;
import org.biojava.nbio.ontology.Ontology;
import org.biojava.nbio.ontology.Term;
import org.biojava.nbio.ontology.Triple;
import org.biojava.nbio.ontology.io.OboParser;

/**
 *
 * @author patri
 */
public class HpoOntology {

	private final Ontology hpoOntologyData;
	private final Term is_a;

	public HpoOntology(File hpoOboFile) throws IOException, FileNotFoundException, ParseException {

		hpoOntologyData = loadHpoOntology(hpoOboFile);
		is_a = hpoOntologyData.getTerm("is_a");

	}
	
	public Set<Term> getAllTerms(){
		return hpoOntologyData.getTerms();
	}

	public List<Term> getChildern(Term queryHpoTerm) {
		ArrayList<Term> hpoChildern = new ArrayList<>();

		for (Triple t : hpoOntologyData.getTriples(null, queryHpoTerm, is_a)) {
			hpoChildern.add(t.getSubject());
		}
		return Collections.unmodifiableList(hpoChildern);
	}
	
	public List<Term> getAllDescendants(Term queryHpoTerm){
		
		ArrayList<Term> descendants = new ArrayList<>();
		
		for (Triple t : hpoOntologyData.getTriples(null, queryHpoTerm, is_a)) {
			Term child = t.getSubject();
			descendants.add(child);
			descendants.addAll(descendants);
		}
		
		return Collections.unmodifiableList(descendants);
	}

	public List<Term> getParents(Term queryHpoTerm) {
		ArrayList<Term> hpoParents = new ArrayList<>();

		for (Triple t : hpoOntologyData.getTriples(queryHpoTerm, null, is_a)) {
			hpoParents.add(t.getObject());
		}
		return Collections.unmodifiableList(hpoParents);
	}

	public Term nameToTerm(String name) {
		if (hpoOntologyData.containsTerm(name)) {
			return hpoOntologyData.getTerm(name);
		} else {
			return null;
		}
	}

	public boolean isAchildOfB(Term termA, Term termB) {
		return !hpoOntologyData.getTriples(termA, termB, is_a).isEmpty();
	}

	public boolean isTermABelowTermB(Term termA, Term termB) {

		if (isAchildOfB(termA, termB)) {
			return true;
		}

		for (Term parentsOfA : getParents(termA)) {

			if (isTermABelowTermB(parentsOfA, termB)) {
				return true;
			}

		}

		return false;

	}

	public static Ontology loadHpoOntology(final File hpoOboFile) throws FileNotFoundException, IOException, ParseException {
		OboParser parser = new OboParser();
		BufferedReader oboFileReader = new BufferedReader(new InputStreamReader(new FileInputStream(hpoOboFile)));
		Ontology hpoOntology = parser.parseOBO(oboFileReader, "HPO", "HPO");
		return hpoOntology;
	}

	public Ontology getHpoOntologyData() {
		return hpoOntologyData;
	}

}
