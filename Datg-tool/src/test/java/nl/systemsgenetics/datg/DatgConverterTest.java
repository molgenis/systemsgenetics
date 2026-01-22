package nl.systemsgenetics.datg;

import umcg.genetica.math.matrix2.DoubleMatrixDataset;
import umcg.genetica.math.matrix2.DoubleMatrixDatasetRowCompressedReader;

import java.io.BufferedWriter;
import java.io.File;
import java.io.FileWriter;
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

        tmpOutputFolder.mkdirs();


        System.out.println("Temp folder with output of this test: " + tmpOutputFolder.getAbsolutePath());


        DatgConverter.TESTNG_MODE = true;

    }

    @org.testng.annotations.AfterMethod
    public void tearDown() {
    }

    @org.testng.annotations.Test
    public void test1() throws Exception {

        final String testName = "test1";

        final DoubleMatrixDataset<String, String> testData1;
        Random random = new Random(42);
        testData1 = new DoubleMatrixDataset<>(1000, 600);
        for (int r = 0; r < testData1.rows(); ++r) {
            for (int c = 0; c < testData1.columns(); ++c) {
                testData1.setElementQuick(r, c, random.nextDouble());
            }
        }


        testData1.save(new File(tmpOutputFolder, testName + ".txt"));

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

        compareTwoMatrices(datgData.loadFullDataset(), testData1, 0.000001);

        DatgConverter.main(new String[]{
                "-m", "DATG_2_TXT",
                "-i", tmpOutputFolder.getPath() + fs + testName + "_step1.datg",
                "-o", tmpOutputFolder.getPath() + fs + testName + "_step2"
        });

        DoubleMatrixDataset<String, String> txtData = DoubleMatrixDataset.loadDoubleTextData(tmpOutputFolder.getPath() + fs + testName + "_step2.txt", '\t');

        compareTwoMatrices(txtData, testData1, 0.000001);

    }

    @org.testng.annotations.Test
    public void test2() throws Exception {

        final String testName = "test2";

        final DoubleMatrixDataset<String, String> testData1;
        Random random = new Random(42);
        testData1 = new DoubleMatrixDataset<>(10000, 1);
        for (int r = 0; r < testData1.rows(); ++r) {
            for (int c = 0; c < testData1.columns(); ++c) {
                testData1.setElementQuick(r, c, random.nextDouble());
            }
        }


        testData1.save(new File(tmpOutputFolder, testName + ".txt.gz"));

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

        compareTwoMatrices(datgData.loadFullDataset(), testData1, 0.000001);

        DatgConverter.main(new String[]{
                "-m", "DATG-2-txt",
                "-i", tmpOutputFolder.getPath() + fs + testName + "_step1",
                "-o", tmpOutputFolder.getPath() + fs + testName + "_step2"
        });

        DoubleMatrixDataset<String, String> txtData = DoubleMatrixDataset.loadDoubleTextData(tmpOutputFolder.getPath() + fs + testName + "_step2.txt", '\t');

        compareTwoMatrices(txtData, testData1, 0.000001);

    }

    @org.testng.annotations.Test
    public void testOldCompressionMethod() throws Exception {

        final String testName = "testOldCompressionMethod";

        final DoubleMatrixDataset<String, String> testData1;
        Random random = new Random(42);
        testData1 = new DoubleMatrixDataset<>(10000, 1);
        for (int r = 0; r < testData1.rows(); ++r) {
            for (int c = 0; c < testData1.columns(); ++c) {
                testData1.setElementQuick(r, c, random.nextDouble());
            }
        }

        DoubleMatrixDatasetRowCompressedWriterOldCompressionMethod.saveDataset(new File(tmpOutputFolder, testName + "_step1").getPath(), testData1, "So sad it had to come to this", "Mice", "Humans");

        DoubleMatrixDatasetRowCompressedReader datgData = new DoubleMatrixDatasetRowCompressedReader(tmpOutputFolder.getPath() + fs + testName + "_step1");

        assertEquals(datgData.getDatasetName(), "So sad it had to come to this");
        assertEquals(datgData.getDataOnRows(), "Mice");
        assertEquals(datgData.getDataOnCols(), "Humans");

        compareTwoMatrices(datgData.loadFullDataset(), testData1, 0.000001);

        DatgConverter.main(new String[]{
                "-m", "DATG-2-txt",
                "-i", tmpOutputFolder.getPath() + fs + testName + "_step1",
                "-o", tmpOutputFolder.getPath() + fs + testName + "_step2"
        });

        DoubleMatrixDataset<String, String> txtData = DoubleMatrixDataset.loadDoubleTextData(tmpOutputFolder.getPath() + fs + testName + "_step2.txt", '\t');

        compareTwoMatrices(txtData, testData1, 0.000001);

    }

    @org.testng.annotations.Test
    public void testUpgrade() throws Exception {

        final String testName = "testUpgrade";

        final DoubleMatrixDataset<String, String> testData1;
        Random random = new Random(42);
        testData1 = new DoubleMatrixDataset<>(100, 5000);
        for (int r = 0; r < testData1.rows(); ++r) {
            for (int c = 0; c < testData1.columns(); ++c) {
                testData1.setElementQuick(r, c, random.nextDouble());
            }
        }

        DoubleMatrixDatasetRowCompressedWriterOldCompressionMethod.saveDataset(new File(tmpOutputFolder, testName).getPath(), testData1, "So sad it had to come to this", "Mice", "Humans");

//        DoubleMatrixDatasetRowCompressedReader datgData = new DoubleMatrixDatasetRowCompressedReader(tmpOutputFolder.getPath() + fs + testName);
//        assertFalse(datgData.isLz4FrameCompression());
//        datgData.close();

        DatgConverter.main(new String[]{
                "-m", "UPGRADE",
                "-i", tmpOutputFolder.getPath() + fs + testName
        });

        DoubleMatrixDatasetRowCompressedReader datgData = new DoubleMatrixDatasetRowCompressedReader(tmpOutputFolder.getPath() + fs + testName);

        assertTrue(datgData.isLz4FrameCompression());
        assertEquals(datgData.getDatasetName(), "So sad it had to come to this");
        assertEquals(datgData.getDataOnRows(), "Mice");
        assertEquals(datgData.getDataOnCols(), "Humans");

        compareTwoMatrices(datgData.loadFullDataset(), testData1, 0.000001);

    }

    @org.testng.annotations.Test
    public void test3() throws Exception {

        final String testName = "test3";

        final DoubleMatrixDataset<String, String> testData1;
        Random random = new Random(42);
        testData1 = new DoubleMatrixDataset<>(10000, 60);
        for (int r = 0; r < testData1.rows(); ++r) {
            for (int c = 0; c < testData1.columns(); ++c) {
                testData1.setElementQuick(r, c, random.nextDouble());
            }
        }


        testData1.saveBinaryOldFormat(new File(tmpOutputFolder, testName).getPath());

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

        compareTwoMatrices(datgData.loadFullDataset(), testData1, 0.000001);

        DatgConverter.main(new String[]{
                "-m", "DATG_2_TXT",
                "-i", tmpOutputFolder.getPath() + fs + testName + "_step1",
                "-o", tmpOutputFolder.getPath() + fs + testName + "_step2"
        });

        DoubleMatrixDataset<String, String> txtData = DoubleMatrixDataset.loadDoubleTextData(tmpOutputFolder.getPath() + fs + testName + "_step2.txt", '\t');

        compareTwoMatrices(txtData, testData1, 0.000001);

    }

    @org.testng.annotations.Test
    public void testRowConcat() throws InterruptedException, IOException {

        final String testName = "testRowConcat";

        LinkedHashMap<String, Integer> hashRows = new LinkedHashMap<String, Integer>();
        LinkedHashMap<String, Integer> hashCols = new LinkedHashMap<String, Integer>();

        for (int r = 0; r < 100; ++r) {
            hashRows.put(("R" + r), r);
        }

        for (int c = 0; c < 60; ++c) {
            hashCols.put(("C" + c), c);
        }

        LinkedHashMap<String, Integer> hashRowsValidation = new LinkedHashMap<String, Integer>();

        int validationR = 0;
        for (String s : new String[]{"ENSG001", "ENSG002", "ENSG003"}) {
            for (int r = 0; r < 100; ++r) {
                hashRowsValidation.put(("R" + r + "-" + s), validationR++);
            }
        }

        final DoubleMatrixDataset<String, String> validationData = new DoubleMatrixDataset<>(hashRowsValidation, hashCols);
        validationR = 0;


        Random random = new Random(42);
        final DoubleMatrixDataset<String, String> testData1;
        testData1 = new DoubleMatrixDataset<>(hashRows, hashCols);
        for (int r = 0; r < testData1.rows(); ++r) {
            for (int c = 0; c < testData1.columns(); ++c) {
                double d = random.nextDouble();
                testData1.setElementQuick(r, c, d);
                validationData.setElementQuick(validationR, c, d);
            }
            validationR++;
        }
        final DoubleMatrixDataset<String, String> testData2 = new DoubleMatrixDataset<>(hashRows, hashCols);
        for (int r = 0; r < testData2.rows(); ++r) {
            for (int c = 0; c < testData2.columns(); ++c) {
                double d = random.nextDouble();
                testData2.setElementQuick(r, c, d);
                validationData.setElementQuick(validationR, c, d);
            }
            validationR++;
        }
        final DoubleMatrixDataset<String, String> testData3 = new DoubleMatrixDataset<>(hashRows, hashCols);
        for (int r = 0; r < testData3.rows(); ++r) {
            for (int c = 0; c < testData3.columns(); ++c) {
                double d = random.nextDouble();
                testData3.setElementQuick(r, c, d);
                validationData.setElementQuick(validationR, c, d);
            }
            validationR++;
        }

        testData1.save(new File(tmpOutputFolder, testName + "zscores_ENSG001" + ".txt.gz"));
        testData2.save(new File(tmpOutputFolder, testName + "zscores_ENSG002" + ".txt.gz"));
        testData3.save(new File(tmpOutputFolder, testName + "zscores_ENSG003" + ".txt.gz"));

        DatgConverter.main(new String[]{
                "-m", "ROW_CONCAT",
                "-i", tmpOutputFolder.getPath(),
                "-o", tmpOutputFolder.getPath() + fs + testName + "_concat1.datg",
                "-fp", "zscores_(.*)\\.txt",
                "-dn", "Row row row",
                "-rc", "My boat",
                "-cc", "Gently down te stream"
        });

        DoubleMatrixDataset<String, String> result = DoubleMatrixDataset.loadDoubleBinaryData(tmpOutputFolder.getPath() + fs + testName + "_concat1.datg");

        compareTwoMatrices(result, validationData, 0.000001);

    }

    @org.testng.annotations.Test
    public void testRowConcat2() throws InterruptedException, IOException {

        final String testName = "testRowConcat2";

        LinkedHashMap<String, Integer> hashRowsValidation = new LinkedHashMap<String, Integer>();


        LinkedHashMap<String, Integer> hashRows1 = new LinkedHashMap<String, Integer>();
        LinkedHashMap<String, Integer> hashCols = new LinkedHashMap<String, Integer>();

        int validationR = 0;
        for (int r = 0; r < 100; ++r) {
            hashRows1.put(("R" + r), r);
            hashRowsValidation.put(("R" + r), validationR++);
        }

        LinkedHashMap<String, Integer> hashRows2 = new LinkedHashMap<String, Integer>();

        for (int r = 0; r < 100; ++r) {
            hashRows2.put(("Arrr" + r), r);
            hashRowsValidation.put(("Arrr" + r), validationR++);
        }

        LinkedHashMap<String, Integer> hashRows3 = new LinkedHashMap<String, Integer>();

        for (int r = 0; r < 100; ++r) {
            hashRows3.put(("Brrr" + r), r);
            hashRowsValidation.put(("Brrr" + r), validationR++);
        }

        for (int c = 0; c < 60; ++c) {
            hashCols.put(("C" + c), c);
        }

        final DoubleMatrixDataset<String, String> validationData = new DoubleMatrixDataset<>(hashRowsValidation, hashCols);
        validationR = 0;


        Random random = new Random(42);
        final DoubleMatrixDataset<String, String> testData1;
        testData1 = new DoubleMatrixDataset<>(hashRows1, hashCols);
        for (int r = 0; r < testData1.rows(); ++r) {
            for (int c = 0; c < testData1.columns(); ++c) {
                double d = random.nextDouble();
                testData1.setElementQuick(r, c, d);
                validationData.setElementQuick(validationR, c, d);
            }
            validationR++;
        }
        final DoubleMatrixDataset<String, String> testData2 = new DoubleMatrixDataset<>(hashRows2, hashCols);
        for (int r = 0; r < testData2.rows(); ++r) {
            for (int c = 0; c < testData2.columns(); ++c) {
                double d = random.nextDouble();
                testData2.setElementQuick(r, c, d);
                validationData.setElementQuick(validationR, c, d);
            }
            validationR++;
        }
        final DoubleMatrixDataset<String, String> testData3 = new DoubleMatrixDataset<>(hashRows3, hashCols);
        for (int r = 0; r < testData3.rows(); ++r) {
            for (int c = 0; c < testData3.columns(); ++c) {
                double d = random.nextDouble();
                testData3.setElementQuick(r, c, d);
                validationData.setElementQuick(validationR, c, d);
            }
            validationR++;
        }

        testData1.save(new File(tmpOutputFolder, testName + "xscores_ENSG001" + ".txt.gz"));
        testData2.save(new File(tmpOutputFolder, testName + "xscores_ENSG002" + ".txt.gz"));
        testData3.save(new File(tmpOutputFolder, testName + "xscores_ENSG003" + ".txt.gz"));

        DatgConverter.main(new String[]{
                "-m", "ROW_CONCAT",
                "-i", tmpOutputFolder.getPath(),
                "-o", tmpOutputFolder.getPath() + fs + testName + "_concat1.datg",
                "-fp", "xscores_.*\\.txt",
                "-dn", "Row row row",
                "-rc", "My boat",
                "-cc", "Gently down te stream"
        });

        DoubleMatrixDataset<String, String> result = DoubleMatrixDataset.loadDoubleBinaryData(tmpOutputFolder.getPath() + fs + testName + "_concat1.datg");

        compareTwoMatrices(result, validationData, 0.000001);

    }

    @org.testng.annotations.Test
    public void testInsepect() throws InterruptedException, IOException {

        //Only stdout so does not really do testing expect that it does not crash. This function was used for development

        final String testName = "testInspect";

        LinkedHashMap<String, Integer> hashRows = new LinkedHashMap<String, Integer>();
        LinkedHashMap<String, Integer> hashCols = new LinkedHashMap<String, Integer>();

        for (int r = 0; r < 10000; ++r) {
            if (r == 2) {
                hashRows.put("SampleWithVeryLongName", r);
            } else {
                hashRows.put(("R" + r), r);
            }
        }

        for (int c = 0; c < 60; ++c) {
            if (c == 3) {
                hashCols.put("ENSG00000136824", c);
            } else if (c == 2) {
                hashCols.put("ENSG00000136824.12", c);
            } else if (c == 1) {
                hashCols.put("ENSG00000136824.12tooLong", c);
            } else {
                hashCols.put(("C" + c), c);
            }

        }

        final DoubleMatrixDataset<String, String> testData1;
        Random random = new Random(42);
        testData1 = new DoubleMatrixDataset<>(hashRows, hashCols);
        for (int r = 0; r < testData1.rows(); ++r) {
            for (int c = 0; c < testData1.columns(); ++c) {
                testData1.setElementQuick(r, c, random.nextDouble());
            }
        }
        testData1.setElementQuick(1, 2, 0);
        testData1.setElementQuick(1, 3, 1);
        testData1.setElementQuick(3, 1, 1000000000001d);
        testData1.setElementQuick(3, 2, 0.0000002341241234d);
        testData1.setElementQuick(4, 0, Double.NaN);
        testData1.setElementQuick(5, 0, -2);
        testData1.setElementQuick(5, 1, -1.45234234);
        testData1.setElementQuick(5, 2, -0.0000044);
        testData1.setElementQuick(6, 0, -9999);
        testData1.setElementQuick(6, 1, -10000);
        testData1.setElementQuick(6, 2, -0.001);
        testData1.setElementQuick(6, 3, -0.0099);
        testData1.setElementQuick(7, 0, -0.00101);
        testData1.setElementQuick(7, 1, -10010);
        testData1.setElementQuick(7, 2, -10001);
        testData1.setElementQuick(8, 0, 8000);
        testData1.setElementQuick(8, 1, 13400);
        testData1.setElementQuick(8, 2, 0.00000123001);

        testData1.save(new File(tmpOutputFolder, testName + ".txt.gz"));

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

    @org.testng.annotations.Test
    public void testGrep() throws Exception {

        final String testName = "testGrep1";

        final DoubleMatrixDataset<String, String> testData1;
        Random random = new Random(42);
        testData1 = new DoubleMatrixDataset<>(1000, 600);
        for (int r = 0; r < testData1.rows(); ++r) {
            for (int c = 0; c < testData1.columns(); ++c) {
                testData1.setElementQuick(r, c, random.nextDouble());
            }
        }

        testData1.save(new File(tmpOutputFolder, testName + ".txt.gz"));

        DatgConverter.main(new String[]{
                "-m", "TXT_2_DATG",
                "-i", tmpOutputFolder.getPath() + fs + testName + ".txt.gz",
                "-o", tmpOutputFolder.getPath() + fs + testName + "_step1",
                "-dn", "So long and thanks for all the fish",
                "-rc", "Fish",
                "-cc", "Dolphins"
        });

        DoubleMatrixDatasetRowCompressedReader datgData = new DoubleMatrixDatasetRowCompressedReader(tmpOutputFolder.getPath() + fs + testName + "_step1");

        assertEquals(datgData.getDatasetName(), "So long and thanks for all the fish");
        assertEquals(datgData.getDataOnRows(), "Fish");
        assertEquals(datgData.getDataOnCols(), "Dolphins");

        compareTwoMatrices(datgData.loadFullDataset(), testData1, 0.000001);


        BufferedWriter grepWriter = new BufferedWriter(new FileWriter(new File(tmpOutputFolder, testName + "grep1.txt")));
        grepWriter.write("R");
        grepWriter.newLine();
        grepWriter.close();

        DatgConverter.main(new String[]{
                "-m", "DATG_2_TXT",
                "-i", tmpOutputFolder.getPath() + fs + testName + "_step1.datg",
                "-o", tmpOutputFolder.getPath() + fs + testName + "_step2",
                "-rg", tmpOutputFolder.getPath() + fs + testName + "grep1.txt"
        });

        DoubleMatrixDataset<String, String> txtData = DoubleMatrixDataset.loadDoubleTextData(tmpOutputFolder.getPath() + fs + testName + "_step2.txt", '\t');

        compareTwoMatrices(txtData, testData1, 0.000001);

    }

    @org.testng.annotations.Test
    public void testGrep2() throws Exception {

        final String testName = "testGrep2";

        final DoubleMatrixDataset<String, String> testData1;
        Random random = new Random(42);
        testData1 = new DoubleMatrixDataset<>(1000, 600);
        for (int r = 0; r < testData1.rows(); ++r) {
            for (int c = 0; c < testData1.columns(); ++c) {
                testData1.setElementQuick(r, c, random.nextDouble());
            }
        }

        testData1.save(new File(tmpOutputFolder, testName + ".txt.gz"));

        DatgConverter.main(new String[]{
                "-m", "TXT_2_DATG",
                "-i", tmpOutputFolder.getPath() + fs + testName + ".txt.gz",
                "-o", tmpOutputFolder.getPath() + fs + testName + "_step1",
                "-dn", "So long and thanks for all the fish",
                "-rc", "Fish",
                "-cc", "Dolphins"
        });

        DoubleMatrixDatasetRowCompressedReader datgData = new DoubleMatrixDatasetRowCompressedReader(tmpOutputFolder.getPath() + fs + testName + "_step1");

        assertEquals(datgData.getDatasetName(), "So long and thanks for all the fish");
        assertEquals(datgData.getDataOnRows(), "Fish");
        assertEquals(datgData.getDataOnCols(), "Dolphins");

        compareTwoMatrices(datgData.loadFullDataset(), testData1, 0.000001);


        BufferedWriter grepWriter = new BufferedWriter(new FileWriter(new File(tmpOutputFolder, testName + "grep1.txt")));
        grepWriter.write("X");
        grepWriter.newLine();
        grepWriter.close();

        DatgConverter.main(new String[]{
                "-m", "DATG_2_TXT",
                "-i", tmpOutputFolder.getPath() + fs + testName + "_step1.datg",
                "-o", tmpOutputFolder.getPath() + fs + testName + "_step2",
                "-rg", tmpOutputFolder.getPath() + fs + testName + "grep1.txt"
        });

        DoubleMatrixDataset<String, String> txtData = DoubleMatrixDataset.loadDoubleTextData(tmpOutputFolder.getPath() + fs + testName + "_step2.txt", '\t');

        assertEquals(txtData.rows(), 0);


        DatgConverter.main(new String[]{
                "-m", "DATG_2_TXT",
                "-i", tmpOutputFolder.getPath() + fs + testName + "_step1.datg",
                "-o", tmpOutputFolder.getPath() + fs + testName + "_step2b",
                "-rg", "X"
        });

        DoubleMatrixDataset<String, String> txtData2 = DoubleMatrixDataset.loadDoubleTextData(tmpOutputFolder.getPath() + fs + testName + "_step2b.txt", '\t');

        assertEquals(txtData2.rows(), 0);


    }

    @org.testng.annotations.Test
    public void testGrep3() throws Exception {

        final String testName = "testGrep3";

        final DoubleMatrixDataset<String, String> testData1;
        Random random = new Random(42);
        testData1 = new DoubleMatrixDataset<>(500, 600);
        for (int r = 0; r < testData1.rows(); ++r) {
            for (int c = 0; c < testData1.columns(); ++c) {
                testData1.setElementQuick(r, c, random.nextDouble());
            }
        }

        testData1.save(new File(tmpOutputFolder, testName + ".txt.gz"));

        DatgConverter.main(new String[]{
                "-m", "TXT_2_DATG",
                "-i", tmpOutputFolder.getPath() + fs + testName + ".txt.gz",
                "-o", tmpOutputFolder.getPath() + fs + testName + "_step1",
                "-dn", "So long and thanks for all the fish",
                "-rc", "Fish",
                "-cc", "Dolphins"
        });

        DoubleMatrixDatasetRowCompressedReader datgData = new DoubleMatrixDatasetRowCompressedReader(tmpOutputFolder.getPath() + fs + testName + "_step1");

        assertEquals(datgData.getDatasetName(), "So long and thanks for all the fish");
        assertEquals(datgData.getDataOnRows(), "Fish");
        assertEquals(datgData.getDataOnCols(), "Dolphins");

        compareTwoMatrices(datgData.loadFullDataset(), testData1, 0.000001);


        BufferedWriter grepWriter = new BufferedWriter(new FileWriter(new File(tmpOutputFolder, testName + "grep1.txt")));
        grepWriter.write("R45");
        grepWriter.newLine();
        grepWriter.write("X");
        grepWriter.newLine();
        grepWriter.write("R101");
        grepWriter.newLine();
        grepWriter.close();

        DatgConverter.main(new String[]{
                "-m", "DATG_2_TXT",
                "-i", tmpOutputFolder.getPath() + fs + testName + "_step1.datg",
                "-o", tmpOutputFolder.getPath() + fs + testName + "_step2",
                "-rg", tmpOutputFolder.getPath() + fs + testName + "grep1.txt"
        });

        DoubleMatrixDataset<String, String> txtData = DoubleMatrixDataset.loadDoubleTextData(tmpOutputFolder.getPath() + fs + testName + "_step2.txt", '\t');

        DoubleMatrixDataset<String, String> testData1b = testData1.viewRowSelection(new String[]{"R45", "R101", "R450", "R451", "R452", "R453", "R454", "R455", "R456", "R457", "R458", "R459"});

        compareTwoMatrices(txtData, testData1b, 0.000001);


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
