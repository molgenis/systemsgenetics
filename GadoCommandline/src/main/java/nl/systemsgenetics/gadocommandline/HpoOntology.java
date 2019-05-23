/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package nl.systemsgenetics.gadocommandline;

import java.io.File;
import java.io.FileNotFoundException;
import java.io.IOException;
import java.text.ParseException;
import java.util.ArrayList;
import java.util.Collections;
import java.util.List;
import org.biojava.nbio.ontology.Ontology;
import org.biojava.nbio.ontology.Term;
import org.biojava.nbio.ontology.Triple;

/**
 *
 * @author patri
 */
public class HpoOntology {

	private final Ontology hpoOntologyData;
	private final Term is_a;

	public HpoOntology(File hpoOboFile) throws IOException, FileNotFoundException, ParseException {

		hpoOntologyData = HpoFinder.loadHpoOntology(hpoOboFile);
		is_a = hpoOntologyData.getTerm("is_a");

	}

	public List<Term> getChildern(Term queryHpoTerm) {
		ArrayList<Term> hpoChildern = new ArrayList<>();

		for (Triple t : hpoOntologyData.getTriples(null, queryHpoTerm, is_a)) {
			hpoChildern.add(t.getSubject());
		}
		return Collections.unmodifiableList(hpoChildern);
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

}
