/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package nl.systemsgenetics.genenetworkbackend.hpo;

import com.opencsv.CSVParser;
import com.opencsv.CSVParserBuilder;
import com.opencsv.CSVReader;
import com.opencsv.CSVReaderBuilder;
import java.io.BufferedReader;
import java.io.File;
import java.io.FileNotFoundException;
import java.io.FileReader;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Collection;
import java.util.Collections;
import java.util.HashMap;
import java.util.HashSet;
import java.util.Set;

/**
 *
 * @author patri
 */
public class DiseaseGeneHpoData {
	
	private final HashMap<String, HashSet<String>> geneToHpos;
	private final HashMap<String, HashSet<String>> diseaseToGenes;

	public DiseaseGeneHpoData(final File diseaseGeneHpoFile, HashMap<String, ArrayList<String>> ncbiToEnsgMap, HashMap<String, ArrayList<String>> hgncToEnsgMap, HashSet<String> exludedHpo ) throws FileNotFoundException, IOException {
		
		geneToHpos = new HashMap<>();
		diseaseToGenes = new HashMap<>();
		
		final CSVParser hpoParser = new CSVParserBuilder().withSeparator('\t').withIgnoreQuotations(true).build();
		final CSVReader hpoReader = new CSVReaderBuilder(new BufferedReader(new FileReader(diseaseGeneHpoFile))).withSkipLines(1).withCSVParser(hpoParser).build();

		String[] nextLine;
		while ((nextLine = hpoReader.readNext()) != null) {
			String disease = nextLine[0];
			String hgcnId = nextLine[1];
			String ncbiId = nextLine[2];
			String hpo = nextLine[3];
			
			if(exludedHpo != null && exludedHpo.contains(hpo)){
				continue;
			}
			
			ArrayList<String> ensgIds = ncbiToEnsgMap.get(ncbiId);
			if (ensgIds == null) {
				ensgIds = hgncToEnsgMap.get(hgcnId);
			}
			if (ensgIds == null) {
				System.err.println("Missing mapping for gene: " + ncbiId + " " + hgcnId);
			} else if(ensgIds.size() > 1) {
				System.err.println("Skipping becasue multiple ENSG IDs for gene: " + ncbiId + " " + hgcnId);				
			}else {
				
				String ensgId = ensgIds.get(0);

				HashSet<String> geneHpos = geneToHpos.get(ensgId);
				if (geneHpos == null) {
					geneHpos = new HashSet<>();
					geneToHpos.put(ensgId, geneHpos);
				}
				
				geneHpos.add(hpo);
				
				HashSet<String> diseaseGenes = diseaseToGenes.get(disease);
				if(diseaseGenes == null){
					diseaseGenes = new HashSet<>();
					diseaseToGenes.put(disease, diseaseGenes);
				}
				diseaseGenes.add(ensgId);

			}

		}
		
	}

	/**
	 * Returns null if no phenotypes associated
	 * 
	 * @param ensgId
	 * @return 
	 */
	public Set<String> getEnsgHpos(String ensgId){
		
		HashSet<String> geneHpos = geneToHpos.get(ensgId);
		
		if(geneHpos == null){
			return null;
		} else {
			return Collections.unmodifiableSet(geneHpos);
		}
		
	}
	
	public Set<String> getDiseaseGenes(){
		return Collections.unmodifiableSet(geneToHpos.keySet());
	}
	
	public Set<String> getDiseases(){
		return Collections.unmodifiableSet(diseaseToGenes.keySet());
	}
	
	/**
	 * Returns null if no disease genes are found
	 * 
	 * @param disease
	 * @return 
	 */
	public Set<String> getGenesForDisease(String disease){
		HashSet<String> diseaseGenes = diseaseToGenes.get(disease);
		
		if(diseaseGenes == null){
			return null;
		} else {
			return Collections.unmodifiableSet(diseaseGenes);
		}
	}
	
}
