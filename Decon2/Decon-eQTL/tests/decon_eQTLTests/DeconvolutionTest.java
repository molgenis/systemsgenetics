package decon_eQTLTests;

import static org.junit.Assert.*;

import java.io.File;
import java.nio.file.Files;
import java.nio.file.Path;
import java.nio.file.Paths;
import java.util.ArrayList;
import java.util.Arrays;

import org.apache.commons.io.FileUtils;
import org.apache.commons.io.LineIterator;
import org.hamcrest.CoreMatchers;
import org.junit.Before;
import org.junit.Test;
import decon_eQTL.Main;

public class DeconvolutionTest {
	String outputDir = "tests/resources/deconvolutionTestResults/";
	
	@Before
	public void init() {
		
	}

	/*
	 *  This is more like an integration test because it runs the whole program!
	 */
	@Test
	public void mainTest() throws Exception {
		File counts = new File("tests/resources/cellcounts.txt");
		File expTable = new File("tests/resources/expression_levels.txt");
		File dsgTable = new File("tests/resources/genotype_dosages.txt");
		File geneSnpList = new File("tests/resources/gene_snp_list.txt");
		String[] args = {"-o",outputDir+"deconvolutionTestResultsMain/","-c",counts.getAbsolutePath(),"-e",
				expTable.getAbsolutePath(), "-g", dsgTable.getAbsolutePath(), "-sn", geneSnpList.getAbsolutePath()};
		Main.main(args);

		LineIterator deconResults = FileUtils.lineIterator(new File(outputDir+"deconvolutionTestResultsMain/deconvolutionResults.csv"), "UTF-8");
		LineIterator deconExpected = FileUtils.lineIterator(new File("tests/resources/deconExpected.txt"), "UTF-8");
		//test if header is same
		assertEquals("File header the same",deconExpected.next(),deconResults.next());
		while (deconResults.hasNext() && deconExpected.hasNext()){
			ArrayList<String> deconResultsStringVector = new ArrayList<String>(Arrays.asList(deconResults.next().split("\t")));
			ArrayList<String> deconExpectedStringVector = new ArrayList<String>(Arrays.asList(deconExpected.next().split("\t")));
			assertEquals("Deconresult same as expected", deconExpectedStringVector, deconResultsStringVector);
			assertEquals("QTL name the same", deconExpectedStringVector.remove(0), deconResultsStringVector.remove(0));
		}
	}

	/*
	 *  I should refractor the code because the only way I can run theses tests is by running from main method
	 *  But for now rather have tests methods than waiting until I have it refractored. So technically not 
	 *  unittests
	 */
	@Test
	public void runDeconPerGeneSnpPairWrongNamesTest() throws Exception {
		File counts = new File("tests/resources/cellcounts.txt");
		File expTable = new File("tests/resources/expression_levels_wrong_names.txt");
		File dsgTable = new File("tests/resources/genotype_dosages.txt");
		File geneSnpList = new File("tests/resources/gene_snp_list.txt");
		String[] args = {"-o",outputDir+"deconvolutionTestResults","-c",counts.getAbsolutePath(),"-e",
				expTable.getAbsolutePath(), "-g", dsgTable.getAbsolutePath(), 
				"-sn", geneSnpList.getAbsolutePath()};

		try {
			Main.main(args);
			fail( "My method didn't throw when I expected it to" );
		} catch (RuntimeException expectedException) {
			assertThat(expectedException.getMessage(), CoreMatchers.containsString("Samplenames not the same in expression and genotype file."));
		}
	}

	/*
	 *  I should refractor the code because the only way I can run theses tests is by running from main method
	 *  But for now rather have tests methods than waiting until I have it refractored. So technically not 
	 *  unittests
	 */
	@Test
	public void runDeconPerGeneSnpPairTestRunTest() throws Exception {
		File counts = new File("tests/resources/cellcounts.txt");
		File expTable = new File("tests/resources/expression_levels.txt");
		File dsgTable = new File("tests/resources/genotype_dosages.txt");
		File geneSnpList = new File("tests/resources/gene_snp_list_long.txt");
		String[] args = {"-o",outputDir+"deconvolutionTestResultsTestRun","-c",counts.getAbsolutePath(),"-e",
				expTable.getAbsolutePath(), "-g", dsgTable.getAbsolutePath(), 
				"-sn", geneSnpList.getAbsolutePath(),"-t"};
		Main.main(args);

		//test if header is same
		Path path = Paths.get(outputDir+"deconvolutionTestResultsTestRun/deconvolutionResults.csv");
		long lineCount = Files.lines(path).count();
		assertEquals("100 example lines written", lineCount, 101);
	}

	/*
	 *  I should refractor the code because the only way I can run theses tests is by running from main method
	 *  But for now rather have tests methods than waiting until I have it refractored. So technically not 
	 *  unittests
	 */
	@Test
	public void runDeconPerGeneSnpPairNotExistingGenotypeTest() throws Exception {
		File counts = new File("tests/resources/cellcounts.txt");
		File expTable = new File("tests/resources/expression_levels.txt");
		File dsgTable = new File("tests/resources/genotype_dosages.txt");
		File geneSnpList = new File("tests/resources/gene_snp_list_non_existing_genotype.txt");
		String[] args = {"-o",outputDir+"deconvolutionTestResults","-c",counts.getAbsolutePath(),"-e",
				expTable.getAbsolutePath(), "-g", dsgTable.getAbsolutePath(), 
				"-sn", geneSnpList.getAbsolutePath(),
		"-t"};

		try {
			Main.main(args);
			fail( "My method didn't throw when I expected it to" );
		} catch (RuntimeException expectedException) {
			assertThat(expectedException.getMessage(), CoreMatchers.containsString("not in genotype file, is your snpsToTest file correct?"));
		}
	}

	/*
	 *  I should refractor the code because the only way I can run theses tests is by running from main method
	 *  But for now rather have tests methods than waiting until I have it refractored. So technically not 
	 *  unittests
	 */
	@Test
	public void runDeconPerGeneSnpPairNotExistingGeneTest() throws Exception {
		File counts = new File("tests/resources/cellcounts.txt");
		File expTable = new File("tests/resources/expression_levels.txt");
		File dsgTable = new File("tests/resources/genotype_dosages.txt");
		File geneSnpList = new File("tests/resources/gene_snp_list_non_existing_gene.txt");
		String[] args = {"-o",outputDir+"deconvolutionTestResults","-c",counts.getAbsolutePath(),"-e",
				expTable.getAbsolutePath(), "-g", dsgTable.getAbsolutePath(), 
				"-sn", geneSnpList.getAbsolutePath(),
		"-t"};

		try {
			Main.main(args);
			fail( "My method didn't throw when I expected it to" );
		} catch (RuntimeException expectedException) {
			assertThat(expectedException.getMessage(), CoreMatchers.containsString("included in gene/snp combinations to test, but not available in the expression file"));
		}
	}
	
	
	/*
	 *  I should refractor the code because the only way I can run theses tests is by running from main method
	 *  But for now rather have tests methods than waiting until I have it refractored. So technically not 
	 *  unittests
	 */
	@Test
	public void runDeconPerGeneSnpPairGeneNotInExpressionFileTest() throws Exception {
		File counts = new File("tests/resources/cellcounts.txt");
		File expTable = new File("tests/resources/expression_levels.txt");
		File dsgTable = new File("tests/resources/genotype_dosages.txt");
		File geneSnpList = new File("tests/resources/gene_snp_list_non_existing_gene.txt");
		String[] args = {"-o",outputDir+"deconvolutionTestResults","-c",counts.getAbsolutePath(),"-e",
				expTable.getAbsolutePath(), "-g", dsgTable.getAbsolutePath(), 
				"-sn", geneSnpList.getAbsolutePath(),
		"-t"};

		try {
			Main.main(args);
			fail( "My method didn't throw when I expected it to" );
		} catch (RuntimeException expectedException) {
			assertThat(expectedException.getMessage(), CoreMatchers.containsString("included in gene/snp combinations to test, but not available in the expression file"));
		}
	}
	
	protected void tearDown() throws Exception {
		deleteDir(new File(outputDir));
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