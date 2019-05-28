/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package nl.systemsgenetics.depict2;

import java.io.File;
import java.net.URISyntaxException;
import java.util.HashMap;
import java.util.HashSet;
import java.util.List;
import static org.testng.Assert.*;
import org.testng.annotations.Test;
import umcg.genetica.math.matrix2.DoubleMatrixDataset;

/**
 *
 * @author patri
 */
public class PathwayEnrichmentsTest {
	
	public PathwayEnrichmentsTest() {
	}

	@org.testng.annotations.BeforeClass
	public static void setUpClass() throws Exception {
	}

	@org.testng.annotations.BeforeMethod
	public void setUpMethod() throws Exception {
	}


	/**
	 * Test of glsFixedInvCor method, of class PathwayEnrichments.
	 */
	@Test
	public void testGlsFixedInvCor() throws URISyntaxException, Exception {

		File invCorFile = new File(this.getClass().getResource("/invCorMatrix.txt").toURI());
		File pathwayFile = new File(this.getClass().getResource("/pathwayGeneScores.txt").toURI());
		File gwasFile = new File(this.getClass().getResource("/gwasGeneScores.txt").toURI());
		
		DoubleMatrixDataset<String, String> geneZscores = DoubleMatrixDataset.loadDoubleTextData(gwasFile.getAbsolutePath(), '\t');
		DoubleMatrixDataset<String, String> genePathwayZscores = DoubleMatrixDataset.loadDoubleTextData(pathwayFile.getAbsolutePath(), '\t');
		DoubleMatrixDataset<String, String> geneInvCor = DoubleMatrixDataset.loadDoubleTextData(invCorFile.getAbsolutePath(), '\t');
		DoubleMatrixDataset expResult = null;
		DoubleMatrixDataset result = PathwayEnrichments.glsFixedInvCor(geneZscores, genePathwayZscores, geneInvCor);
		
		geneInvCor.printMatrix();
		
		result.printMatrix();
		
		//assertEquals(result, expResult);
	
	}
	
}
