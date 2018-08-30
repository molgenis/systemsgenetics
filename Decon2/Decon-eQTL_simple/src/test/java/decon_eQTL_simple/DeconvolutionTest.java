package test.java.decon_eQTL_simple;

import static org.junit.Assert.*;
import java.io.File;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;

import org.apache.commons.io.FileUtils;
import org.apache.commons.io.LineIterator;
import org.hamcrest.CoreMatchers;
import org.junit.After;
import org.junit.Before;
import org.junit.Test;

import main.java.decon_eQTL_simple.DeconvolutionResult;
import main.java.decon_eQTL_simple.CommandLineOptions;
import main.java.decon_eQTL_simple.Deconvolution;
import main.java.decon_eQTL_simple.Main;

public class DeconvolutionTest {
	String outputDir = "src/test/resources/deconvolutionTestResults/";
	String counts;
	String genotypes;
	String geneSnpList;
	String expression;
	CommandLineOptions commandLineOptions;
	
	@Before
	public void init() {
		commandLineOptions = new CommandLineOptions();
		File countsFile = new File("src/test/resources/cellcount_files/cellcounts.txt");
		File genotypesFile = new File("src/test/resources/genotype_files/genotype_dosages.txt");
		File geneSnpListFile = new File("src/test/resources/gene_snp_list_files/gene_snp_list.txt");
		File expressionFile = new File("src/test/resources/expression_files/expression_levels.txt");
		counts = countsFile.getAbsolutePath();
		genotypes = genotypesFile.getAbsolutePath();
		geneSnpList = geneSnpListFile.getAbsolutePath();
		expression = expressionFile.getAbsolutePath();
	}
	
	@After
	public void tearDown() throws Exception {
		//deleteDir(new File(outputDir));
	}	
	
	@Test
	public void mainTest() throws Exception {
		String[] args = {"-o",outputDir+"deconvolutionTestResultsTestRun","-c",counts,
						 "-e",expression, "-g",genotypes, "-sn", geneSnpList};

		Main.main(args);
	}
	
	/*
	 * Give error when expression and genotype file have different names
	 */
	@Test
	public void readInputDataTest() throws Exception {
		File cellCountsSmall = new File("src/test/resources/cellcount_files/cellcounts_small.txt");
		String[] args = {"-o",outputDir+"deconvolutionTestResults","-c",counts,cellCountsSmall.getAbsolutePath(),
						 "-e",expression, "-g", genotypes, 
						 "-sn", geneSnpList};
		commandLineOptions.parseCommandLine(args);
		Deconvolution deconvolution = new Deconvolution(commandLineOptions);

		deconvolution.readInputData();
	}
	
	/*
	
	/*
	 * Give error when expression and genotype file have different names
	 */
	@Test
	public void readInputDataWrongNamesTest() throws Exception {

		File expTableWrongNames = new File("src/test/resources/expression_files/expression_levels_wrong_names.txt");
		String[] args = {"-o",outputDir+"deconvolutionTestResults","-c",counts,"-e",
						 expTableWrongNames.getAbsolutePath(), "-g", genotypes, 
						 "-sn", geneSnpList};
		commandLineOptions.parseCommandLine(args);
		Deconvolution deconvolution = new Deconvolution(commandLineOptions);
		try {
			deconvolution.readInputData();
			fail( "My method didn't throw when I expected it to" );
		} catch (RuntimeException expectedException) {
			assertThat(expectedException.getMessage(), CoreMatchers.containsString("Samplenames not the same in expression and genotype file, or not in the same order"));
		}
	}

	@Test
	public void runDeconPerGeneSnpPairNotExistingGenotypeTest() throws Exception {
		File geneSnpList = new File("src/test/resources/gene_snp_list_files/gene_snp_list_non_existing_genotype.txt");
		String[] args = {"-o",outputDir+"deconvolutionTestResults","-c",counts,
						 "-e",expression, "-g", genotypes, 
						 "-sn", geneSnpList.getAbsolutePath()};
		commandLineOptions.parseCommandLine(args);
		Deconvolution deconvolution = new Deconvolution(commandLineOptions);
		deconvolution.readInputData();
		try {
			deconvolution.runDeconPerGeneSnpPair();
			fail( "My method didn't throw when I expected it to" );
		} catch (RuntimeException expectedException) {
			assertEquals(expectedException.getMessage(),
						"Error: Genotype genotype_NOT_EXISTING included in gene/snp combinations to test, but not available in the expression file!");
		}
	}

	@Test
	public void runDeconPerGeneSnpPairNotExistingGeneTest() throws Exception {
		File geneSnpList = new File("src/test/resources/gene_snp_list_files/gene_snp_list_non_existing_gene.txt");
		String[] args = {"-o",outputDir+"deconvolutionTestResults","-c",counts,
						 "-e",expression, "-g", genotypes, 
						 "-sn", geneSnpList.getAbsolutePath()};
		commandLineOptions.parseCommandLine(args);
		Deconvolution deconvolution = new Deconvolution(commandLineOptions);
		deconvolution.readInputData();
		try {
			deconvolution.runDeconPerGeneSnpPair();
			fail( "My method didn't throw when I expected it to" );
		} catch (RuntimeException expectedException) {
			assertThat(expectedException.getMessage(), CoreMatchers.containsString("included in gene/snp combinations to test, but not available in the expression file"));
		}
	}
	
	@Test
	public void runDeconPerGeneSnpPairGenotypeNotInGenotypeFileTest() throws Exception {
		File geneSnpList = new File("src/test/resources/gene_snp_list_files/gene_snp_list_non_existing_genotype.txt");
		String[] args = {"-o",outputDir+"deconvolutionTestResults","-c",counts,
						 "-e",expression, "-g", genotypes,
						 "-sn", geneSnpList.getAbsolutePath()};
		commandLineOptions.parseCommandLine(args);
		Deconvolution deconvolution = new Deconvolution(commandLineOptions);
		deconvolution.readInputData();
		try {
			deconvolution.runDeconPerGeneSnpPair();
			fail( "My method didn't throw when I expected it to" );
		} catch (RuntimeException expectedException) {
			assertEquals(expectedException.getMessage(), 
						"Error: Genotype genotype_NOT_EXISTING included in gene/snp combinations to test, but not available in the expression file!");
		}
	}
		
	@Test
	public void runDeconPerGeneSnpPairGeneNotInExpressionFileTest() throws Exception {
		File geneSnpList = new File("src/test/resources/gene_snp_list_files/gene_snp_list_non_existing_gene.txt");
		String[] args = {"-o",outputDir+"deconvolutionTestResults","-c",counts,
						 "-e",expression, "-g", genotypes,
						 "-sn", geneSnpList.getAbsolutePath()};
		commandLineOptions.parseCommandLine(args);
		Deconvolution deconvolution = new Deconvolution(commandLineOptions);
		deconvolution.readInputData();
		try {
			deconvolution.runDeconPerGeneSnpPair();
			fail( "My method didn't throw when I expected it to" );
		} catch (RuntimeException expectedException) {
			assertThat(expectedException.getMessage(), CoreMatchers.containsString("included in gene/snp combinations to test, but not available in the expression file"));
		}
	}

	@Test
	public void writeDeconvolutionResultTest() throws Exception {
		System.out.println("TEST");
		String[] args = {"-o",outputDir+"deconvolutionResults","-c",counts,
						 "-e",expression, "-g", genotypes,
						 "-sn", geneSnpList};
		commandLineOptions.parseCommandLine(args);
		Deconvolution deconvolution = new Deconvolution(commandLineOptions);
		deconvolution.readInputData();
		List<DeconvolutionResult> deconvolutionResults = deconvolution.runDeconPerGeneSnpPair();
		deconvolution.writeDeconvolutionResults(deconvolutionResults);

		LineIterator deconResults = FileUtils.lineIterator(new File(outputDir+"deconvolutionResults/deconvolutionResults.csv"), "UTF-8");
		LineIterator deconExpected = FileUtils.lineIterator(new File("src/test/resources/expected_results/deconExpected.txt"), "UTF-8");
		//test if header is same
		assertEquals("File header the same",deconExpected.next(),deconResults.next());
		while (deconResults.hasNext() && deconExpected.hasNext()){
			ArrayList<String> deconResultsStringVector = new ArrayList<String>(Arrays.asList(deconResults.next().split("\t")));
			ArrayList<String> deconExpectedStringVector = new ArrayList<String>(Arrays.asList(deconExpected.next().split("\t")));
			assertEquals("Deconresult same as expected", deconExpectedStringVector, deconResultsStringVector);
			assertEquals("QTL name the same", deconExpectedStringVector.remove(0), deconResultsStringVector.remove(0));
		}
	}
	
	
	void deleteDir(File file) {
	    File[] contents = file.listFiles();
	    if (contents != null) {
	        for (File f : contents) {
	            deleteDir(f);
	        }
	    }
	    file.delete();
	}
}