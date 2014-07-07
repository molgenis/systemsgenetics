/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package eqtlmappingpipeline;

import java.io.File;
import java.net.URISyntaxException;
import java.nio.charset.Charset;
import java.text.DateFormat;
import java.text.SimpleDateFormat;
import java.util.Date;
import org.molgenis.genotype.Allele;
import org.molgenis.genotype.Alleles;
import static org.testng.Assert.*;
import org.testng.annotations.BeforeMethod;
import org.testng.annotations.BeforeTest;
import org.testng.annotations.Test;
/**
 *
 * @author MarcJan
 */
public class FullQtlMappingTransMetaTest {
    private File tmpOutputFolder;
	private String fileSep = System.getProperty("file.separator");
	private File testFilesFolder;
	private static final Charset FILE_ENCODING = Charset.forName("UTF-8");
	private static final Alleles MISSING_ALLELES = Alleles.createAlleles(Allele.ZERO, Allele.ZERO);
    
    public FullQtlMappingTransMetaTest() throws URISyntaxException {
		testFilesFolder = new File(this.getClass().getResource("/GeuvadisTestData/").toURI());
	}

	@BeforeTest
	public void setUpMethod() throws Exception {
		File tmpDir = new File(System.getProperty("java.io.tmpdir"));

		DateFormat dateFormat = new SimpleDateFormat("yyyy-MM-dd_HH-mm-ss");
		Date date = new Date();

		tmpOutputFolder = new File(tmpDir, "QTLMappingTransMetaTest_" + dateFormat.format(date));

		Runtime.getRuntime().addShutdownHook(new Thread() {
			@Override
			public void run() {
				System.out.println("Removing tmp dir and files");
				for (File file : tmpOutputFolder.listFiles()) {
					System.out.println(" - Deleting: " + file.getAbsolutePath());
//					file.delete();
				}
				System.out.println(" - Deleting: " + tmpOutputFolder.getAbsolutePath());
//				tmpOutputFolder.delete();
			}
		});

		tmpOutputFolder.mkdir();


		System.out.println("Temp folder with output of this test: " + tmpOutputFolder.getAbsolutePath());
	}
    
    /**
	 * Test of main method, of class GenotypeHarmonizer.
	 */
	@Test
	public void testMain() throws Exception {
		String setingsFile = testFilesFolder + fileSep + "settings.xml";
		System.out.println(setingsFile);

		Main.main("--mode", "metaqtl", "--settings", setingsFile, "--replacetext", "${InputFolder}", "--replacetextwith", testFilesFolder.getAbsolutePath(), "--replacetext2", "${OutputFolder}", "--replacetext2with", tmpOutputFolder.getAbsolutePath());
    }
}
