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
import java.util.Iterator;
import org.molgenis.genotype.Allele;
import org.molgenis.genotype.Alleles;
import static org.testng.Assert.*;
import org.testng.annotations.BeforeMethod;
import org.testng.annotations.BeforeTest;
import org.testng.annotations.Test;
import umcg.genetica.io.trityper.EQTL;
import umcg.genetica.io.trityper.eQTLTextFile;
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
					file.delete();
				}
				System.out.println(" - Deleting: " + tmpOutputFolder.getAbsolutePath());
				tmpOutputFolder.delete();
			}
		});

		tmpOutputFolder.mkdir();


		System.out.println("Temp folder with output of this test: " + tmpOutputFolder.getAbsolutePath());
	}
    
	@Test
	public void testMain() throws Exception {
		String setingsFile = testFilesFolder + fileSep + "settings.xml";
		System.out.println(setingsFile);

		Main.main("--mode", "metaqtl", "--settings", setingsFile, "--replacetext", "${InputFolder}", "--replacetextwith", testFilesFolder.getAbsolutePath(), "--replacetext2", "${OutputFolder}", "--replacetext2with", tmpOutputFolder.getAbsolutePath());
        
        //Read in expected results
        eQTLTextFile eExp = new eQTLTextFile(testFilesFolder+fileSep+"TestOutput"+fileSep+"Trans-Meta-eQTLProbesFDR0.05-ProbeLevel.txt", eQTLTextFile.R);
        eQTLTextFile eActual = new eQTLTextFile(tmpOutputFolder.getAbsolutePath()+fileSep+"eQTLProbesFDR0.05-ProbeLevel.txt", eQTLTextFile.R);

        Iterator<EQTL> eExpIterator = eExp.getEQtlIterator();
        Iterator<EQTL> eActualIterator = eActual.getEQtlIterator();
        
        while(eExpIterator.hasNext() && eActualIterator.hasNext()){
            assertTrue(eActualIterator.next().sameQTL(eExpIterator.next()), "eQTL not identical");
        }
        
        assertFalse(eExpIterator.hasNext(), "not all expected eQTL are found");
        assertFalse(eActualIterator.hasNext(), "found more eQTL than expected");
    }
}
