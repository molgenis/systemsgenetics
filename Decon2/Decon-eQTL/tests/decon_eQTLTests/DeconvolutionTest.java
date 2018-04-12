package decon_eQTLTests;

import static org.junit.Assert.*;

import java.io.File;
import java.util.ArrayList;
import java.util.Arrays;

import org.apache.commons.io.FileUtils;
import org.apache.commons.io.LineIterator;
import org.junit.Before;
import org.junit.Test;
import decon_eQTL.Deconvolution;

public class DeconvolutionTest {
    @Before
    public void init() {
    	// Read in the example data
    }
        
	@Test
	public void mainTest() throws Exception {
		File counts = new File("tests/resources/cellcounts.txt");
		File expTable = new File("tests/resources/expression_levels.txt");
		File dsgTable = new File("tests/resources/genotype_dosages.txt");
		File geneSnpList = new File("tests/resources/gene_snp_list.txt");
		String[] args = {"-o","tests/resources/decovolutionTestResult","-c",counts.getAbsolutePath(),"-e",
						 expTable.getAbsolutePath(), "-g", dsgTable.getAbsolutePath(), "-sn", geneSnpList.getAbsolutePath()};
		Deconvolution.main(args);
		LineIterator deconResults = FileUtils.lineIterator(new File("tests/resources/decovolutionTestResult/deconvolutionResults.csv"), "UTF-8");
		LineIterator deconExpected = FileUtils.lineIterator(new File("tests/resources/deconExpected.txt"), "UTF-8");
		//test if header is same
		assertEquals("File header the same",deconExpected.next(),deconResults.next());
		while (deconResults.hasNext() && deconExpected.hasNext()){
			ArrayList<String> deconResultsStringVector = new ArrayList<String>(Arrays.asList(deconResults.next().split("\t")));
			ArrayList<String> deconExpectedStringVector = new ArrayList<String>(Arrays.asList(deconExpected.next().split("\t")));
			assertEquals("QTL name the same", deconExpectedStringVector.remove(0), deconResultsStringVector.remove(0));
		}
	}
}