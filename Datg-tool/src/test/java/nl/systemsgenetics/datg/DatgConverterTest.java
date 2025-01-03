package nl.systemsgenetics.datg;

import umcg.genetica.math.matrix2.DoubleMatrixDataset;
import umcg.genetica.math.matrix2.DoubleMatrixDatasetRowCompressedReader;

import java.io.File;
import java.io.IOException;
import java.nio.file.FileSystems;
import java.text.DateFormat;
import java.text.SimpleDateFormat;
import java.util.Date;
import java.util.LinkedHashMap;
import java.util.Random;

import static org.testng.Assert.*;

public class DatgConverterTest {

    private File tmpOutputFolder;
    private DoubleMatrixDataset<String, String> testData;
    private final String fs = FileSystems.getDefault().getSeparator();

    @org.testng.annotations.BeforeClass
    public void setUp() {

        File tmpDir = new File(System.getProperty("java.io.tmpdir"));

        DateFormat dateFormat = new SimpleDateFormat("yyyy-MM-dd_HH-mm-ss");
        Date date = new Date();

        tmpOutputFolder = new File(tmpDir, "GenotyperHarmonizerTest_" + dateFormat.format(date));

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




        DatgConverter.TESTNG_MODE = true;

    }

    @org.testng.annotations.AfterMethod
    public void tearDown() {
    }

    @org.testng.annotations.Test
    public void test1() throws Exception {

        final String testName = "test1";

        Random random = new Random(42);
        testData = new DoubleMatrixDataset<>(1000,600);
        for (int r = 0; r < testData.rows() ; ++r){
            for(int c = 0; c < testData.columns() ; ++c){
                testData.setElementQuick(r,c,random.nextDouble());
            }
        }


        testData.save(new File(tmpOutputFolder, testName + ".txt"));

        DatgConverter.main(new String[]{
                "-m", "TXT_2_DATG",
                "-i", tmpOutputFolder.getPath() + fs + testName + ".txt",
                "-o", tmpOutputFolder.getPath() + fs + testName + "_step1",
                "-dn", "So long and thanks for all the fish",
                "-rc", "Fish",
                "-cc", "Dolphins"
        });

        DoubleMatrixDatasetRowCompressedReader datgData = new DoubleMatrixDatasetRowCompressedReader(tmpOutputFolder.getPath() + fs + testName + "_step1");

        assertEquals(datgData.getDatasetName(), "So long and thanks for all the fish");
        assertEquals(datgData.getDataOnRows(), "Fish");
        assertEquals(datgData.getDataOnCols(), "Dolphins");

        compareTwoMatrices(datgData.loadFullDataset(), testData, 0.000001);

        DatgConverter.main(new String[]{
                "-m", "DATG_2_TXT",
                "-i", tmpOutputFolder.getPath() + fs + testName + "_step1.datg",
                "-o", tmpOutputFolder.getPath() + fs + testName + "_step2"
        });

        DoubleMatrixDataset<String, String> txtData = DoubleMatrixDataset.loadDoubleTextData(tmpOutputFolder.getPath() + fs + testName + "_step2.txt", '\t');

        compareTwoMatrices(txtData, testData, 0.000001);

    }

    @org.testng.annotations.Test
    public void test2() throws Exception {

        final String testName = "test2";

        Random random = new Random(42);
        testData = new DoubleMatrixDataset<>(10000,1);
        for (int r = 0; r < testData.rows() ; ++r){
            for(int c = 0; c < testData.columns() ; ++c){
                testData.setElementQuick(r,c,random.nextDouble());
            }
        }


        testData.save(new File(tmpOutputFolder, testName + ".txt.gz"));

        DatgConverter.main(new String[]{
                "-m", "TXT_2_DATG",
                "-i", tmpOutputFolder.getPath() + fs + testName + ".txt.gz",
                "-o", tmpOutputFolder.getPath() + fs + testName + "_step1.datg",
                "-dn", "So sad it had to come to this",
                "-rc", "Mice",
                "-cc", "Humans"
        });

        DoubleMatrixDatasetRowCompressedReader datgData = new DoubleMatrixDatasetRowCompressedReader(tmpOutputFolder.getPath() + fs + testName + "_step1");

        assertEquals(datgData.getDatasetName(), "So sad it had to come to this");
        assertEquals(datgData.getDataOnRows(), "Mice");
        assertEquals(datgData.getDataOnCols(), "Humans");

        compareTwoMatrices(datgData.loadFullDataset(), testData, 0.000001);

        DatgConverter.main(new String[]{
                "-m", "DATG-2-txt",
                "-i", tmpOutputFolder.getPath() + fs + testName + "_step1",
                "-o", tmpOutputFolder.getPath() + fs + testName + "_step2"
        });

        DoubleMatrixDataset<String, String> txtData = DoubleMatrixDataset.loadDoubleTextData(tmpOutputFolder.getPath() + fs + testName + "_step2.txt", '\t');

        compareTwoMatrices(txtData, testData, 0.000001);

    }

    @org.testng.annotations.Test
    public void test3() throws Exception {

        final String testName = "test3";

        Random random = new Random(42);
        testData = new DoubleMatrixDataset<>(10000,60);
        for (int r = 0; r < testData.rows() ; ++r){
            for(int c = 0; c < testData.columns() ; ++c){
                testData.setElementQuick(r,c,random.nextDouble());
            }
        }


        testData.saveBinaryOldFormat(new File(tmpOutputFolder, testName).getPath());

        DatgConverter.main(new String[]{
                "-m", "DAT_2_DATG",
                "-i", tmpOutputFolder.getPath() + fs + testName,
                "-o", tmpOutputFolder.getPath() + fs + testName + "_step1",
                "-dn", "Vogon poetry",
                "-rc", "Please",
                "-cc", "Stop"
        });

        DoubleMatrixDatasetRowCompressedReader datgData = new DoubleMatrixDatasetRowCompressedReader(tmpOutputFolder.getPath() + fs + testName + "_step1");

        assertEquals(datgData.getDatasetName(), "Vogon poetry");
        assertEquals(datgData.getDataOnRows(), "Please");
        assertEquals(datgData.getDataOnCols(), "Stop");

        compareTwoMatrices(datgData.loadFullDataset(), testData, 0.000001);

        DatgConverter.main(new String[]{
                "-m", "DATG_2_TXT",
                "-i", tmpOutputFolder.getPath() + fs + testName + "_step1",
                "-o", tmpOutputFolder.getPath() + fs + testName + "_step2"
        });

        DoubleMatrixDataset<String, String> txtData = DoubleMatrixDataset.loadDoubleTextData(tmpOutputFolder.getPath() + fs + testName + "_step2.txt", '\t');

        compareTwoMatrices(txtData, testData, 0.000001);

    }

    @org.testng.annotations.Test
    public void testInsepect() throws InterruptedException, IOException {

        //Only stdout so does not really do testing expect that it does not crash. This function was used for development

        final String testName = "testInspect";

        LinkedHashMap<String, Integer> hashRows = new LinkedHashMap<String, Integer>();
        LinkedHashMap<String, Integer> hashCols = new LinkedHashMap<String, Integer>();

        for(int r = 0 ; r < 10000 ; ++r){
            if(r == 2){
                hashRows.put("SampleWithVeryLongName", r);
            } else {
                hashRows.put(("R" + r), r);
            }
        }

        for(int c = 0 ; c < 60 ; ++c){
            if(c == 3){
                hashCols.put("ENSG00000136824", c);
            } else if(c == 2){
                hashCols.put("ENSG00000136824.12", c);
            } else if(c == 1){
                hashCols.put("ENSG00000136824.12tooLong", c);
            } else {
                hashCols.put(("C" + c), c);
            }

        }


        Random random = new Random(42);
        testData = new DoubleMatrixDataset<>(hashRows, hashCols);
        for (int r = 0; r < testData.rows() ; ++r){
            for(int c = 0; c < testData.columns() ; ++c){
                testData.setElementQuick(r,c,random.nextDouble());
            }
        }
        testData.setElementQuick(1,2, 0);
        testData.setElementQuick(1,3, 1);
        testData.setElementQuick(3,1, 1000000000001d);
        testData.setElementQuick(3,2, 0.0000002341241234d);
        testData.setElementQuick(4,0, Double.NaN);
        testData.setElementQuick(5,0, -2);
        testData.setElementQuick(5,1, -1.45234234);
        testData.setElementQuick(5,2, -0.0000044);
        testData.setElementQuick(6,0, -9999);
        testData.setElementQuick(6,1, -10000);
        testData.setElementQuick(6,2, -0.001);
        testData.setElementQuick(6,3, -0.0099);
        testData.setElementQuick(7,0, -0.00101);
        testData.setElementQuick(7,1, -10010);
        testData.setElementQuick(7,2, -10001);
        testData.setElementQuick(8,0, 8000);
        testData.setElementQuick(8,1, 13400);
        testData.setElementQuick(8,2, 0.00000123001);

        testData.save(new File(tmpOutputFolder, testName + ".txt.gz"));

        DatgConverter.main(new String[]{
                "-m", "TXT_2_DATG",
                "-i", tmpOutputFolder.getPath() + fs + testName + ".txt.gz",
                "-o", tmpOutputFolder.getPath() + fs + testName + "_step1.datg",
                "-dn", "Inspector gadget",
                "-rc", "Genes",
                "-cc", "Samples"
        });

        DatgConverter.main(new String[]{
                "-m", "INSPECT",
                "-i", tmpOutputFolder.getPath() + fs + testName + "_step1"
        });

    }

    public static void compareTwoMatrices(DoubleMatrixDataset<String, String> m1, DoubleMatrixDataset<String, String> m2, double delta) {

        assertEquals(m1.rows(), m2.rows());
        assertEquals(m1.columns(), m2.columns());

        assertEquals(m1.getRowObjects(), m2.getRowObjects());
        assertEquals(m1.getColObjects(), m2.getColObjects());

        for (int r = 0; r < m1.rows(); ++r) {
            for (int c = 0; c < m1.columns(); ++c) {
                assertEquals(m1.getElementQuick(r, c), m2.getElementQuick(r, c), delta, "Difference at r: " + r + " c: " + c);
            }
        }

    }

}
