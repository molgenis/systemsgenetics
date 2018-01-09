/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package eqtlmappingpipeline;

import eqtlmappingpipeline.util.QTLFileSorter;
import java.io.File;
import java.net.URISyntaxException;
import java.text.DateFormat;
import java.text.SimpleDateFormat;
import java.util.Date;
import java.util.Iterator;
import static org.testng.Assert.*;
import org.testng.annotations.BeforeTest;
import org.testng.annotations.Test;
import umcg.genetica.io.trityper.EQTL;
import umcg.genetica.io.trityper.QTLTextFile;

/**
 *
 * @author MarcJan
 */
public class FullQtlMappingCisTest {

    private File tmpOutputFolder;
    private final String fileSep = System.getProperty("file.separator");
    private final File testFilesFolder;

    public FullQtlMappingCisTest() throws URISyntaxException {
        testFilesFolder = new File(this.getClass().getResource("/GeuvadisTestData/").toURI());
    }

    @BeforeTest
    public void setUpMethod() throws Exception {
        File tmpDir = new File(System.getProperty("java.io.tmpdir"));

        DateFormat dateFormat = new SimpleDateFormat("yyyy-MM-dd_HH-mm-ss");
        Date date = new Date();

        tmpOutputFolder = new File(tmpDir, "QTLMappingCisTest_" + dateFormat.format(date));

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
        String inputDir = testFilesFolder + fileSep + "trityper";
        String inputExprs = testFilesFolder + fileSep + "Geuvadis_CEU_YRI_Expr.txt.gz";
        String inputExprsAnnot = testFilesFolder + fileSep + "Geuvadis_CEU_YRI_Annot.txt";
        String inputGte = testFilesFolder + fileSep + "Geuvadis_CEU_gte.txt";

        System.out.println(inputDir);

        Main.main("--mode", "metaqtl", "--in", inputDir, "--out", tmpOutputFolder.getAbsolutePath(), "--cis", "--perm", "10", "--inexp", inputExprs, "--inexpannot", inputExprsAnnot, "--inexpplatform", "Ensembl_v.71", "--gte", inputGte, "--skipqqplot", "--skipdotplot", "--rseed", "0");

        QTLTextFile eExp = new QTLTextFile(testFilesFolder + fileSep + "TestOutput" + fileSep + "Cis-CEU-eQTLsFDR0.05.txt", QTLTextFile.R);
        
        QTLFileSorter r = new QTLFileSorter();
        r.run(tmpOutputFolder.getAbsolutePath() + fileSep + "eQTLsFDR0.05.txt.gz", tmpOutputFolder.getAbsolutePath() + fileSep + "eQTLsFDR0.05_S.txt.gz");
        
        QTLTextFile eActual = new QTLTextFile(tmpOutputFolder.getAbsolutePath() + fileSep + "eQTLsFDR0.05_S.txt.gz", QTLTextFile.R);
        Iterator<EQTL> eExpIterator = eExp.getEQtlIterator();
        Iterator<EQTL> eActualIterator = eActual.getEQtlIterator();
        
        while(eExpIterator.hasNext() && eActualIterator.hasNext()){
            assertTrue(eActualIterator.next().sameQTL(eExpIterator.next()), "eQTL not identical");
        }
        
        assertFalse(eExpIterator.hasNext(), "not all expected eQTL are found");
        assertFalse(eActualIterator.hasNext(), "found more eQTL than expected");
        
        
    }
    
}
