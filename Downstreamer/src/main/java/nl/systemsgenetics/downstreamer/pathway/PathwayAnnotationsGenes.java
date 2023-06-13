/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package nl.systemsgenetics.downstreamer.pathway;

import java.io.File;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Collections;
import java.util.HashMap;
import java.util.List;
import nl.systemsgenetics.downstreamer.gene.Gene;
import nl.systemsgenetics.downstreamer.io.IoUtils;

/**
 *
 * @author patri
 */
public class PathwayAnnotationsGenes implements PathwayAnnotations{

	private final HashMap<String, Gene> annotations;
	private static final List<String> header;
	
	static {
		
		ArrayList<String> headerTmp = new ArrayList<>(5);
		headerTmp.add("Gene symbol");
		headerTmp.add("Chr");
		headerTmp.add("Band");
		headerTmp.add("Start");
		headerTmp.add("End");
		header = Collections.unmodifiableList(headerTmp);
		
	}
	
	public PathwayAnnotationsGenes(final File geneInfoFile) throws IOException {
		this.annotations= IoUtils.readGenesMap(geneInfoFile);
	}
	
	@Override
	public List<String> getAnnotationHeaders() {
		return header;
	}

	@Override
	public List<String> getAnnotationsForPathway(String pathway) {
		Gene gene = annotations.get(pathway);
		ArrayList<String> annotationsTmp = new ArrayList<>(5);
		annotationsTmp.add(gene.getGeneSymbol());
		annotationsTmp.add(gene.getContig());
		annotationsTmp.add(gene.getBand());
		annotationsTmp.add(String.valueOf(gene.getStart()));
		annotationsTmp.add(String.valueOf(gene.getEnd()));
		return Collections.unmodifiableList(annotationsTmp);
	}

	@Override
	public int getMaxNumberOfAnnotations() {
		return 5;
	}

	@Override
	public String getSetName() {
		return "Gene ID";
	}
	
}
