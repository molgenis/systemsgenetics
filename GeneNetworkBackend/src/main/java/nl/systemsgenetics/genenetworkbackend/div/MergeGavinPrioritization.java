/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package nl.systemsgenetics.genenetworkbackend.div;

import com.opencsv.CSVParser;
import com.opencsv.CSVParserBuilder;
import com.opencsv.CSVReader;
import com.opencsv.CSVReaderBuilder;
import com.opencsv.CSVWriter;
import java.io.BufferedReader;
import java.io.File;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;
import java.util.HashSet;
import java.util.regex.Matcher;
import java.util.regex.Pattern;

/**
 *
 * @author patri
 */
public class MergeGavinPrioritization {

	/**
	 * @param args the command line arguments
	 */
	public static void main(String[] args) throws IOException {

		final File sampleFile = new File("C:\\UMCG\\Genetica\\Projects\\GeneNetwork\\BenchmarkSamples\\PrioritisationsDcm\\samplesWithGeno.txt");
		final File genoFolder = new File("C:\\UMCG\\Genetica\\Projects\\GeneNetwork\\BenchmarkSamples\\PrioritisationsDcm\\gavinRes\\");
		final File prioFolder = new File("C:\\UMCG\\Genetica\\Projects\\GeneNetwork\\BenchmarkSamples\\PrioritisationsDcm\\rankingCandidateGenes");
		final File outputFile = new File("C:\\UMCG\\Genetica\\Projects\\GeneNetwork\\BenchmarkSamples\\PrioritisationsDcm\\mergedResults.txt");

		String DNA_no;
		String Chromosome;
		String Start;
		String Stop;
		String Transcript;
		String OMIM;
		String A1;
		String A2;
		String Genotype;
		String Allelic_Depth_Allele_1;
		String Allelic_Depth_Allele_2;
		String Mother_Allele_1;
		String Mother_Allele_2;
		String Father_Allele_1;
		String Father_Allele_2;
		String Gene;
		String cDNA;
		String HGVS_protein_level_nomenclature;
		String Opmerking;
		String Artefact;
		String Overgeerfd_van;
		String DN;
		String HGMD_accession_number;
		String HGMD_mutation_type;
		String HGMD_variant_type;
		String HGMD_dissease;
		String Pop_Freq_gnomAD;
		String HMZ_Hemi;
		String PhyloP_score;
		String Eiwit_conservering;
		String t_m;
		String Mogelijke_splice_mutatie;
		String CADD_score;
		String RLV;
		String Pseudogen__bekend;
		String Classificatie_volgens_beslisboom;
		String OMIM2;
		String PubMed;
		String Expressie;
		String Geno2MP_HPO;
		String Mouse_phenotype;
		String Conclusie;

		CSVWriter writer = new CSVWriter(new FileWriter(outputFile), '\t', '\0', '\0', "\n");

		String[] outputLine = new String[42];
		int c = 0;
		outputLine[c++] = "DNA_no";
		outputLine[c++] = "Chromosome";
		outputLine[c++] = "Start";
		outputLine[c++] = "Stop";
		outputLine[c++] = "Transcript";
		outputLine[c++] = "OMIM";
		outputLine[c++] = "A1";
		outputLine[c++] = "A2";
		outputLine[c++] = "Genotype";
		outputLine[c++] = "Allelic_Depth_Allele_1";
		outputLine[c++] = "Allelic_Depth_Allele_2";
		outputLine[c++] = "Mother_Allele_1";
		outputLine[c++] = "Mother_Allele_2";
		outputLine[c++] = "Father_Allele_1";
		outputLine[c++] = "Father_Allele_2";
		outputLine[c++] = "Gene";
		outputLine[c++] = "cDNA";
		outputLine[c++] = "HGVS_protein_level_nomenclature";
		outputLine[c++] = "Opmerking";
		outputLine[c++] = "Artefact";
		outputLine[c++] = "Overgeerfd_van";
		outputLine[c++] = "DN";
		outputLine[c++] = "HGMD_accession_number";
		outputLine[c++] = "HGMD_mutation_type";
		outputLine[c++] = "HGMD_variant_type";
		outputLine[c++] = "HGMD_dissease";
		outputLine[c++] = "Pop_Freq_exac";
		outputLine[c++] = "HMZ_Hemi";
		outputLine[c++] = "PhyloP_score";
		outputLine[c++] = "Eiwit_conservering";
		outputLine[c++] = "t_m";
		outputLine[c++] = "Mogelijke_splice_mutatie";
		outputLine[c++] = "CADD_score";
		outputLine[c++] = "RLV";
		outputLine[c++] = "Pseudogen__bekend";
		outputLine[c++] = "Classificatie_volgens_beslisboom";
		outputLine[c++] = "OMIM2";
		outputLine[c++] = "PubMed";
		outputLine[c++] = "Expressie";
		outputLine[c++] = "Geno2MP_HPO";
		outputLine[c++] = "Mouse_phenotype";
		outputLine[c++] = "Conclusie";
		writer.writeNext(outputLine);

		final CSVParser parser = new CSVParserBuilder().withSeparator('\t').withIgnoreQuotations(true).build();
		final CSVReader sampleFileReader = new CSVReaderBuilder(new BufferedReader(new FileReader(sampleFile))).withSkipLines(0).withCSVParser(parser).build();

		Pattern pattern = Pattern.compile("CADD score of (.*) is");

		String[] nextLine;
		while ((nextLine = sampleFileReader.readNext()) != null) {

			String sample = nextLine[0];

			String genoSampleName = new File(nextLine[1]).getName();
			if (!genoSampleName.endsWith(".txt")) {
				genoSampleName += ".txt";
			}

			File genoFile = new File(genoFolder, genoSampleName);
			File prioFile = new File(prioFolder, sample + ".txt");

			System.out.println("------------------------------------------------------------------");
			System.out.println("Sample: " + sample);
			System.out.println("Geno: " + genoFile.getAbsolutePath());
			System.out.println("Prio: " + prioFile.getAbsolutePath());

			HashSet<String> genesWithStrongZ = new HashSet<>();

			final CSVReader gnFileReader = new CSVReaderBuilder(new BufferedReader(new FileReader(prioFile))).withSkipLines(1).withCSVParser(parser).build();
			while ((nextLine = gnFileReader.readNext()) != null) {

				String gene = nextLine[1];
				double zscore = Double.parseDouble(nextLine[3]);
				if (zscore >= 3) {
					genesWithStrongZ.add(gene);
				}

			}

			final CSVReader gavinFileReader = new CSVReaderBuilder(new BufferedReader(new FileReader(genoFile))).withSkipLines(1).withCSVParser(parser).build();
			while ((nextLine = gavinFileReader.readNext()) != null) {

				Gene = nextLine[8];

				if (genesWithStrongZ.contains(Gene)) {

					String[] genotypeInfo = nextLine[5].split(":");
//					
//					String[] depths = genotypeInfo[2].split(",");
//					
//					System.out.println(Gene);
//					System.out.println(genotypeInfo[2]);
					System.out.println(Gene);
					Matcher matcher = pattern.matcher(nextLine[22]);
					
					DNA_no = sample;
					Chromosome = nextLine[0];
					Start = nextLine[1];
					Stop = "";
					Transcript = nextLine[10];
					OMIM = "";
					A1 = nextLine[3];
					A2 = nextLine[4];
					Genotype = genotypeInfo[0];
					Allelic_Depth_Allele_1 = "";
					Allelic_Depth_Allele_2 = "";
					Mother_Allele_1 = "";
					Mother_Allele_2 = "";
					Father_Allele_1 = "";
					Father_Allele_2 = "";
					cDNA = "";
					HGVS_protein_level_nomenclature = "";
					Opmerking = nextLine[11];
					Artefact = "";
					Overgeerfd_van = "";
					DN = "";
					HGMD_accession_number = "";
					HGMD_mutation_type = "";
					HGMD_variant_type = "";
					HGMD_dissease = "";
					Pop_Freq_gnomAD = nextLine[7];
					HMZ_Hemi = "";
					PhyloP_score = "";
					Eiwit_conservering = "";
					t_m = "";
					Mogelijke_splice_mutatie = "";
					CADD_score =  matcher.find() ? matcher.group(1) : "";
					RLV = nextLine[22];
					Pseudogen__bekend = "";
					Classificatie_volgens_beslisboom = "";
					OMIM2 = "";
					PubMed = "";
					Expressie = "";
					Geno2MP_HPO = "";
					Mouse_phenotype = "";
					Conclusie = "";

					c = 0;
					outputLine[c++] = DNA_no;
					outputLine[c++] = Chromosome;
					outputLine[c++] = Start;
					outputLine[c++] = Stop;
					outputLine[c++] = Transcript;
					outputLine[c++] = OMIM;
					outputLine[c++] = A1;
					outputLine[c++] = A2;
					outputLine[c++] = Genotype;
					outputLine[c++] = Allelic_Depth_Allele_1;
					outputLine[c++] = Allelic_Depth_Allele_2;
					outputLine[c++] = Mother_Allele_1;
					outputLine[c++] = Mother_Allele_2;
					outputLine[c++] = Father_Allele_1;
					outputLine[c++] = Father_Allele_2;
					outputLine[c++] = Gene;
					outputLine[c++] = cDNA;
					outputLine[c++] = HGVS_protein_level_nomenclature;
					outputLine[c++] = Opmerking;
					outputLine[c++] = Artefact;
					outputLine[c++] = Overgeerfd_van;
					outputLine[c++] = DN;
					outputLine[c++] = HGMD_accession_number;
					outputLine[c++] = HGMD_mutation_type;
					outputLine[c++] = HGMD_variant_type;
					outputLine[c++] = HGMD_dissease;
					outputLine[c++] = Pop_Freq_gnomAD;
					outputLine[c++] = HMZ_Hemi;
					outputLine[c++] = PhyloP_score;
					outputLine[c++] = Eiwit_conservering;
					outputLine[c++] = t_m;
					outputLine[c++] = Mogelijke_splice_mutatie;
					outputLine[c++] = CADD_score;
					outputLine[c++] = RLV;
					outputLine[c++] = Pseudogen__bekend;
					outputLine[c++] = Classificatie_volgens_beslisboom;
					outputLine[c++] = OMIM2;
					outputLine[c++] = PubMed;
					outputLine[c++] = Expressie;
					outputLine[c++] = Geno2MP_HPO;
					outputLine[c++] = Mouse_phenotype;
					outputLine[c++] = Conclusie;
					writer.writeNext(outputLine);

				}
			}
		}

		writer.close();
	}

}
