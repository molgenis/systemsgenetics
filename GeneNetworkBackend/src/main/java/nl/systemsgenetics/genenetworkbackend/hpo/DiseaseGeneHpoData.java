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
import java.util.Collections;
import java.util.HashMap;
import java.util.HashSet;
import java.util.Map;
import java.util.Objects;
import java.util.Random;
import java.util.Set;
import umcg.genetica.math.matrix2.DoubleMatrixDataset;

/**
 *
 * @author patri
 */
public class DiseaseGeneHpoData {

	private final HashMap<String, HashSet<String>> geneToHpos;
	private final HashMap<String, HashSet<String>> diseaseToGenes;
	private final HashMap<DiseaseGene, HashSet<String>> diseaseGeneToHpos; // disease_gene

	public DiseaseGeneHpoData(final File diseaseGeneHpoFile, HashMap<String, ArrayList<String>> ncbiToEnsgMap, HashMap<String, ArrayList<String>> hgncToEnsgMap, HashSet<String> exludedHpo) throws FileNotFoundException, IOException {

		geneToHpos = new HashMap<>();
		diseaseToGenes = new HashMap<>();
		diseaseGeneToHpos = new HashMap<>();

		final CSVParser hpoParser = new CSVParserBuilder().withSeparator('\t').withIgnoreQuotations(true).build();
		final CSVReader hpoReader = new CSVReaderBuilder(new BufferedReader(new FileReader(diseaseGeneHpoFile))).withSkipLines(1).withCSVParser(hpoParser).build();

		String[] nextLine;
		while ((nextLine = hpoReader.readNext()) != null) {
			String disease = nextLine[0];
			String hgcnId = nextLine[1];
			String ncbiId = nextLine[2];
			String hpo = nextLine[3];

			if (exludedHpo != null && exludedHpo.contains(hpo)) {
				continue;
			}

			ArrayList<String> ensgIds = ncbiToEnsgMap.get(ncbiId);
			if (ensgIds == null) {
				ensgIds = hgncToEnsgMap.get(hgcnId);
			}
			if (ensgIds == null) {
				System.err.println("Missing mapping for gene: " + ncbiId + " " + hgcnId);
			} else if (ensgIds.size() > 1) {
				System.err.println("Skipping becasue multiple ENSG IDs for gene: " + ncbiId + " " + hgcnId);
			} else {

				String ensgId = ensgIds.get(0);

				HashSet<String> geneHpos = geneToHpos.get(ensgId);
				if (geneHpos == null) {
					geneHpos = new HashSet<>();
					geneToHpos.put(ensgId, geneHpos);
				}

				geneHpos.add(hpo);

				HashSet<String> diseaseGenes = diseaseToGenes.get(disease);
				if (diseaseGenes == null) {
					diseaseGenes = new HashSet<>();
					diseaseToGenes.put(disease, diseaseGenes);
				}
				diseaseGenes.add(ensgId);

				DiseaseGene diseaseGene = new DiseaseGene(disease, ensgId);

				HashSet<String> diseaseGeneHpos = diseaseGeneToHpos.get(diseaseGene);
				if (diseaseGeneHpos == null) {
					diseaseGeneHpos = new HashSet<>();
					diseaseGeneToHpos.put(diseaseGene, diseaseGeneHpos);
				}
				diseaseGeneHpos.add(hpo);

			}

		}

	}

	public DiseaseGeneHpoData(HashMap<DiseaseGene, HashSet<String>> diseaseGeneToHpos) {

		this.diseaseGeneToHpos = diseaseGeneToHpos;

		geneToHpos = new HashMap<>();
		diseaseToGenes = new HashMap<>();

		for (Map.Entry<DiseaseGene, HashSet<String>> diseaseGeneToHposEntry : diseaseGeneToHpos.entrySet()) {

			DiseaseGene diseaseGene = diseaseGeneToHposEntry.getKey();
			HashSet<String> hpos = diseaseGeneToHposEntry.getValue();

			HashSet<String> geneHpos = geneToHpos.get(diseaseGene.getGene());
			if (geneHpos == null) {
				geneHpos = new HashSet<>();
				geneToHpos.put(diseaseGene.getGene(), geneHpos);
			}

			geneHpos.addAll(hpos);

			HashSet<String> diseaseGenes = diseaseToGenes.get(diseaseGene.getDisease());
			if (diseaseGenes == null) {
				diseaseGenes = new HashSet<>();
				diseaseToGenes.put(diseaseGene.getDisease(), diseaseGenes);
			}
			diseaseGenes.add(diseaseGene.getGene());

		}

	}

	/**
	 * Returns null if no phenotypes associated
	 *
	 * @param ensgId
	 * @return
	 */
	public Set<String> getEnsgHpos(String ensgId) {

		HashSet<String> geneHpos = geneToHpos.get(ensgId);

		if (geneHpos == null) {
			return null;
		} else {
			return Collections.unmodifiableSet(geneHpos);
		}

	}

	public Set<String> getDiseaseGenes() {
		return Collections.unmodifiableSet(geneToHpos.keySet());
	}

	public Set<String> getDiseases() {
		return Collections.unmodifiableSet(diseaseToGenes.keySet());
	}

	public Set<DiseaseGene> getDiseaseGeneHpos() {
		return Collections.unmodifiableSet(diseaseGeneToHpos.keySet());
	}

	/**
	 * Returns null if no disease genes are found
	 *
	 * @param disease
	 * @return
	 */
	public Set<String> getGenesForDisease(String disease) {
		HashSet<String> diseaseGenes = diseaseToGenes.get(disease);

		if (diseaseGenes == null) {
			return null;
		} else {
			return Collections.unmodifiableSet(diseaseGenes);
		}
	}

	/**
	 * Returns null if no phenotypes associated
	 *
	 * @param diseaseGene disease_gene
	 * @return
	 */
	public Set<String> getDiseaseEnsgHpos(DiseaseGene diseaseGene) {

		HashSet<String> hpos = diseaseGeneToHpos.get(diseaseGene);

		if (hpos == null) {
			return null;
		} else {
			return Collections.unmodifiableSet(hpos);
		}

	}

	public DiseaseGeneHpoData getPermutation() {
		return getPermutation(new Random(), null, null, 0);
	}

	public DiseaseGeneHpoData getPermutation(long seed) {
		return getPermutation(new Random(seed), null, null, 0);
	}
	
	public DiseaseGeneHpoData getPermutation(long seed, ArrayList<String> backgroundGenes) {
		return getPermutation(new Random(seed), backgroundGenes, null, 0);
	}
	
	public DiseaseGeneHpoData getPermutation(ArrayList<String> backgroundGenes) {
		return getPermutation(new Random(), backgroundGenes, null, 0);
	}
	
	public DiseaseGeneHpoData getPermutation(long seed, ArrayList<String> backgroundGenes, DoubleMatrixDataset<String, String> predictionMatrixSignificantCorrelationMatrix, double minCorrelationTomatch) {
		return getPermutation(new Random(seed), backgroundGenes, null, 0);
	}

	private DiseaseGeneHpoData getPermutation(Random random, ArrayList<String> backgroundGenes, DoubleMatrixDataset<String, String> predictionMatrixSignificantCorrelationMatrix, double minCorrelationTomatch) {

		if(backgroundGenes == null){
			backgroundGenes = new ArrayList(geneToHpos.keySet());
		}
		

		HashMap<DiseaseGene, HashSet<String>> randomDiseaseGeneToHpos = new HashMap<>();

		for (Map.Entry<DiseaseGene, HashSet<String>> diseaseGeneToHposEntry : this.diseaseGeneToHpos.entrySet()) {

			DiseaseGene diseaseGene = diseaseGeneToHposEntry.getKey();
			HashSet<String> hpos = diseaseGeneToHposEntry.getValue();

			String disease = diseaseGene.getDisease();

			HashSet<String> knownGenesForDisease = this.diseaseToGenes.get(disease);

			String randomReplacementGene;
			DiseaseGene randomDiseaseGene = null;
			boolean hpoOverlap;
			boolean hpoCorrelated;
			
			int i = 0;
			boolean noRandomFound = false;
			
			do {
				
				if(i++ >= 50000){
					System.err.println("No random match found");
					noRandomFound = true;
					break;
				}
				
				randomReplacementGene = backgroundGenes.get(random.nextInt(backgroundGenes.size()));
				randomDiseaseGene = new DiseaseGene(disease, randomReplacementGene);
				HashSet<String> knownHposForRandomGene = this.geneToHpos.get(randomReplacementGene);
				hpoOverlap = false;
				if (knownHposForRandomGene != null) {
					for (String hpo : hpos) {
						if (knownHposForRandomGene.contains(hpo)) {
							hpoOverlap = true;
							break;
						}
					}
				}
				hpoCorrelated = false;
				if(!hpoOverlap && predictionMatrixSignificantCorrelationMatrix != null && knownHposForRandomGene != null){
					//if already hpo overlap no need to do this
					
					hposLoop:
					for (String hpo : hpos) {
						
						for(String randomHpo : knownHposForRandomGene){
							
							if(predictionMatrixSignificantCorrelationMatrix.getElement(hpo, randomHpo) >= minCorrelationTomatch){
								hpoCorrelated = true;
								break hposLoop;
							}
							
						}
						
					}
					
					
				}
			} while (hpoCorrelated | hpoOverlap | knownGenesForDisease.contains(randomReplacementGene) | randomDiseaseGeneToHpos.containsKey(randomDiseaseGene));

			if(!noRandomFound){
				randomDiseaseGeneToHpos.put(randomDiseaseGene, hpos);
			}
			

		}

		return new DiseaseGeneHpoData(randomDiseaseGeneToHpos);

	}

	public class DiseaseGene {

		private final String disease;
		private final String gene;

		public DiseaseGene(String disease, String gene) {
			this.disease = disease;
			this.gene = gene;
		}

		public String getDisease() {
			return disease;
		}

		public String getGene() {
			return gene;
		}

		@Override
		public int hashCode() {
			int hash = 3;
			hash = 97 * hash + Objects.hashCode(this.disease);
			hash = 97 * hash + Objects.hashCode(this.gene);
			return hash;
		}

		@Override
		public boolean equals(Object obj) {
			if (this == obj) {
				return true;
			}
			if (obj == null) {
				return false;
			}
			if (getClass() != obj.getClass()) {
				return false;
			}
			final DiseaseGene other = (DiseaseGene) obj;
			if (!Objects.equals(this.disease, other.disease)) {
				return false;
			}
			if (!Objects.equals(this.gene, other.gene)) {
				return false;
			}
			return true;
		}

		@Override
		public String toString() {
			return disease + "_" + gene;
		}

	}

}
