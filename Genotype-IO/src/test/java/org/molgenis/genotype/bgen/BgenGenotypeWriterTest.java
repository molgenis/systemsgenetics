package org.molgenis.genotype.bgen;

import org.molgenis.genotype.GenotypeData;
import org.molgenis.genotype.ResourceTest;
import org.molgenis.genotype.RandomAccessGenotypeDataReaderFormats;
import org.molgenis.genotype.plink.BedBimFamGenotypeData;
import org.molgenis.genotype.util.GenotypeDataCompareTool;
import org.testng.annotations.BeforeClass;
import org.testng.annotations.BeforeMethod;
import org.testng.annotations.Test;

import java.io.File;
import java.io.IOException;
import java.net.URISyntaxException;
import java.text.DateFormat;
import java.text.SimpleDateFormat;
import java.util.Date;

import static org.testng.Assert.*;

public class BgenGenotypeWriterTest extends ResourceTest {

    private GenotypeData genotypeData;
    private BgenGenotypeWriter writer;
    private File tmpOutputFolder;
    private String fileSep = System.getProperty("file.separator");

    @BeforeClass
    public void setUp() throws IOException, URISyntaxException {

        File tmpDir = new File(System.getProperty("java.io.tmpdir"));

        DateFormat dateFormat = new SimpleDateFormat("yyyy-MM-dd_HH-mm-ss");
        Date date = new Date();

        tmpOutputFolder = new File(tmpDir, "GenotypeBgenTest_" + dateFormat.format(date));


        Runtime.getRuntime().addShutdownHook(new Thread() {
            @Override
            public void run() {
                System.out.println("Removing tmp dir and files");
                for (File file : tmpOutputFolder.listFiles()) {
                    System.out.println(" - Deleting: " + file.getAbsolutePath());
                    //file.delete();
                }
                System.out.println(" - Deleting: " + tmpOutputFolder.getAbsolutePath());
                tmpOutputFolder.delete();
            }
        });

        tmpOutputFolder.mkdir();


        genotypeData = new BedBimFamGenotypeData(getTestBed6(), getTestBim6(), getTestFam6(), 2);
        writer = new BgenGenotypeWriter(genotypeData);
    }

    @Test
    public void write() throws IOException, URISyntaxException {

        writer.write(tmpOutputFolder.getAbsolutePath() + fileSep + "test");

        GenotypeData genotypeDataWritten = RandomAccessGenotypeDataReaderFormats.BGEN.createGenotypeData(tmpOutputFolder.getAbsolutePath() + fileSep + "test", 0);

        assertTrue(GenotypeDataCompareTool.same(genotypeData, genotypeDataWritten));
    }
}