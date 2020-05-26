/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package nl.systemsgenetics.eqtlpermutationtranscriptionfactoranalysis;

import java.util.regex.Pattern;

/**
 *
 * @author Matthieu
 */
public class eQtlsInEncodePseudogenes {
	private static final Pattern TAB_PATTERN = Pattern.compile("\t");
	private static final Pattern SEMICOLON_PATTERN = Pattern.compile(";");
	private static final Pattern SPACE_PATTERN = Pattern.compile("\\s{1}");
	
	public void getGeneName(String annotationLine){
		String[] annotationData = SEMICOLON_PATTERN.split(annotationLine);
		String geneNameLine = new String(annotationData[4]);
		//TODO:
		
	}
	
}
