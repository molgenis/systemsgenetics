/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package umcg.genetica.methylation;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileInputStream;
import java.io.IOException;
import java.io.InputStreamReader;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.logging.Level;
import java.util.logging.Logger;
import java.util.regex.Pattern;
import umcg.genetica.io.text.TextFile;
import umcg.genetica.math.matrix.DoubleMatrixDataset;

/**
 *
 * @author Marc Jan and Lude
 */
public class ParseTcgaFile {

    private static Pattern SPLIT_ON_TAB = Pattern.compile("\\t");
    protected static final String ENCODING = "ISO-8859-1";
    
    /**
     * Read Level 1 data from TGCA a single DoubleMatrixDatasets containing beta
     * values. Beta values are calculated as : M/U+M (Non-bead studio / genome
     * studio).
     *
     * @param fileInputFolder
     * @param printToFile
     * @param fileOutputFolder
     * @param TcgaMethod
     * @return all data out of the TCGA files Methylated / un-methylated and
     * beta-values
     */
    public static DoubleMatrixDataset<String, String> parseTCGAData_lvl1(String fileInputFolder, boolean printToFile, String fileOutputFolder, boolean TcgaMethod) {

        File file = new File(fileInputFolder);
        File[] files = file.listFiles();
        ArrayList<File> vecFiles = new ArrayList<File>();
        for (int f = 0; f < files.length; f++) {
            if (files[f].getAbsolutePath().endsWith(".txt")) {
                vecFiles.add(files[f]);
            }
        }

        System.out.println("Files to parse:\t" + vecFiles.size());
        int nrSamples = vecFiles.size();

        int nrProbes = 0;
        ArrayList<String> vecProbes = new ArrayList<String>();
        try {
            BufferedReader in = new BufferedReader(new InputStreamReader(new FileInputStream(vecFiles.get(0)), ENCODING), 8096);
            String str;
            in.readLine();
            in.readLine();
            while ((str = in.readLine()) != null) {
                String[] data = SPLIT_ON_TAB.split(str);
                vecProbes.add(data[0]);
                nrProbes++;
            }
            in.close();
        } catch (IOException e) {
            System.out.println(e.getMessage());
            System.exit(-1);
        }

        System.out.println(nrProbes);

        DoubleMatrixDataset<String, String> dataset3 = new DoubleMatrixDataset<String, String>(nrProbes, nrSamples);
        //ExpressionDataset dataset3 = new ExpressionDataset(nrProbes * 2, nrSamples);
        for (int p = 0; p < vecProbes.size(); p++) {
            dataset3.rowObjects.set(p, vecProbes.get(p));
            //dataset3.probeNames[p * 2] = "M-" + (String) vecProbes.get(p);
            //dataset3.probeNames[p * 2 + 1] = "U-" + (String) vecProbes.get(p);
        }
        for (int f = 0; f < nrSamples; f++) {
            File currentFile = vecFiles.get(f);
            //System.out.println("Processing:\t" + f + "\t" + currentFile.getAbsolutePath());
            try {
                BufferedReader in = new BufferedReader(new InputStreamReader(new FileInputStream(currentFile), ENCODING), 8096);
                String str = in.readLine();
                String[] data = SPLIT_ON_TAB.split(str);
                dataset3.colObjects.set(f, data[1]);
                str = in.readLine();
                data = SPLIT_ON_TAB.split(str);
                int columnM = -1;
                int columnU = -1;
                for (int d = 0; d < data.length; d++) {
                    if (data[d].toLowerCase().trim().replace(" ", "_").equals("methylated_signal_intensity_(m)")) {
                        columnM = d;
                    }
                    if (data[d].toLowerCase().trim().replace(" ", "_").equals("un-methylated_signal_intensity_(u)")) {
                        columnU = d;
                    }
                }
                //System.out.println(columnM + "\t" + columnU);
                int p = 0;
                while ((str = in.readLine()) != null) {
                    data = SPLIT_ON_TAB.split(str);
                    if (data[columnM].equals("NA") || data[columnM].equals("NaN")) {
                        data[columnM] = "-999";
                    }
                    if (data[columnU].equals("NA") || data[columnU].equals("NaN")) {
                        data[columnU] = "-999";
                    }
                    double methylatedSignal = Double.parseDouble(data[columnM]);
                    double unmethylatedSignal = Double.parseDouble(data[columnU]);

                    //Non bead-studio / genome-studio beta calculations. TCGA way.
                    if (methylatedSignal == -999 || unmethylatedSignal == -999) {
                        dataset3.rawData[p][f] = -999;
                    } else if(methylatedSignal == 0 && unmethylatedSignal == 0){
                        dataset3.rawData[p][f] = 0;
                    } else if (methylatedSignal <= 0 || unmethylatedSignal <= 0) {
                        dataset3.rawData[p][f] = -999;
                    } else if (TcgaMethod) {
                        dataset3.rawData[p][f] = methylatedSignal / (methylatedSignal + unmethylatedSignal);
                    } else {
                        dataset3.rawData[p][f] = methylatedSignal / ((methylatedSignal + unmethylatedSignal) + 100);
                    }

                    String probe = (String) vecProbes.get(p);
                    if (!data[0].equals(probe)) {
                        System.out.println("Error!:\t" + f + "\t" + data[0] + "\t" + probe);
                    }
                    p++;
                }
                in.close();
            } catch (IOException e) {
                System.out.println(e.getMessage());
                System.exit(-1);
            }

        }

        dataset3.recalculateHashMaps();

        if (printToFile) {

            try {
                dataset3.save(fileOutputFolder + "/TCGADataBeta.txt");
            } catch (IOException ex) {
                Logger.getLogger(ParseTcgaFile.class.getName()).log(Level.SEVERE, null, ex);
            }
        }
        ArrayList<DoubleMatrixDataset<String, String>> tmp = new ArrayList<DoubleMatrixDataset<String, String>>(3);

        return (dataset3);
    }

    /**
     * Read Level 1 450K data from TGCA a single DoubleMatrixDatasets containing beta
     * values. Beta values are calculated as : M/U+M (Non-bead studio / genome
     * studio).
     *
     * @param fileInputFolder
     * @param printToFile
     * @param fileOutputFolder
     * @param TcgaMethod
     * @return all data out of the TCGA files Methylated / un-methylated and
     * beta-values
     */
    public static DoubleMatrixDataset<String, String> parseTCGAData450As27K_lvl1(String fileInputFolder, boolean printToFile, String fileOutputFolder, boolean TcgaMethod, HashMap<String, Boolean> probeList) {
        HashMap<String, Integer> probeIndex = new HashMap<String, Integer>();
        File file = new File(fileInputFolder);
        File[] files = file.listFiles();
        ArrayList<File> vecFiles = new ArrayList<File>();
        for (int f = 0; f < files.length; f++) {
            if (files[f].getAbsolutePath().endsWith(".txt")) {
                vecFiles.add(files[f]);
            }
        }

        System.out.println("Files to parse:\t" + vecFiles.size());
        int nrSamples = vecFiles.size();

        int nrProbes = 0;
        ArrayList<String> vecProbes = new ArrayList<String>();
        try {
            BufferedReader in = new BufferedReader(new InputStreamReader(new FileInputStream(vecFiles.get(0)), ENCODING), 8096);
            String str;
            in.readLine();
            in.readLine();
            while ((str = in.readLine()) != null) {
                String[] data = SPLIT_ON_TAB.split(str);
                vecProbes.add(data[0]);
                nrProbes++;
            }
            in.close();
        } catch (IOException e) {
            System.out.println(e.getMessage());
            System.exit(-1);
        }

        System.out.println(nrProbes + "\t" + probeList.size());

        DoubleMatrixDataset<String, String> dataset3 = new DoubleMatrixDataset<String, String>(probeList.size(), nrSamples);
        //ExpressionDataset dataset3 = new ExpressionDataset(nrProbes * 2, nrSamples);
        int position = 0;
        for (int p = 0; p < vecProbes.size(); p++) {
            if (probeList.containsKey(vecProbes.get(p))) {
                dataset3.rowObjects.set(position, vecProbes.get(p));
                probeIndex.put(vecProbes.get(p), position);
                position++;
            }
        }

        for (int f = 0; f < nrSamples; f++) {
            String currentFile = vecFiles.get(f).getAbsolutePath();
            System.out.println("Processing:\t" + f + "\t" + currentFile);

            int columnM = 1;
            int columnU = 2;

            try {
                TextFile in = new TextFile(currentFile, TextFile.R);

                String str = in.readLine();
                String[] data = SPLIT_ON_TAB.split(str);
                dataset3.colObjects.set(f, data[1]);
                in.readLine();

                //System.out.println(columnM + "\t" + columnU);
                int nrSet = 0;
                while ((str = in.readLine()) != null && nrSet != probeIndex.size()) {
                    data = SPLIT_ON_TAB.split(str);

                    if (probeIndex.containsKey(data[0])) {
                        nrSet++;
                        int p = probeIndex.get(data[0]);
                        double methylatedSignal;
                        double unmethylatedSignal;
                        if (data[columnM].equals("NA") || data[columnM].equals("NaN")) {
                            methylatedSignal = -999d;
                        } else {
                            methylatedSignal = Double.parseDouble(data[columnM]);
                        }
                        if (data[columnU].equals("NA") || data[columnU].equals("NaN")) {
                            unmethylatedSignal = -999d;
                        } else {
                            unmethylatedSignal = Double.parseDouble(data[columnU]);
                        }

                        //Non bead-studio / genome-studio beta calculations. TCGA way.
                        if (methylatedSignal == -999 || unmethylatedSignal == -999) {
                            dataset3.rawData[p][f] = -999;
                        } else if(methylatedSignal == 0 && unmethylatedSignal == 0){
                            dataset3.rawData[p][f] = 0;
                        } else if (methylatedSignal <= 0 || unmethylatedSignal <= 0) {
                            dataset3.rawData[p][f] = -999;
                        } else if (TcgaMethod) {
                            dataset3.rawData[p][f] = methylatedSignal / (methylatedSignal + unmethylatedSignal);
                        } else {
                            dataset3.rawData[p][f] = methylatedSignal / ((methylatedSignal + unmethylatedSignal) + 100);
                        }
                        
//                        if (!data[0].equals(dataset3.rowObjects.get(p))) {
//                            System.out.println("Error!:\t" + f + "\t" + data[0] + "\t" + dataset3.rowObjects.get(p));
//                        }
                    }
                }
                in.close();
            } catch (IOException e) {
                System.out.println(e.getMessage());
                System.exit(-1);
            }

        }

        dataset3.recalculateHashMaps();

        if (printToFile) {

            try {
                dataset3.save(fileOutputFolder + "TCGA_450K-27K_DataBeta.txt");
            } catch (IOException ex) {
                Logger.getLogger(ParseTcgaFile.class.getName()).log(Level.SEVERE, null, ex);
            }
        }
        ArrayList<DoubleMatrixDataset<String, String>> tmp = new ArrayList<DoubleMatrixDataset<String, String>>(3);

        return (dataset3);
    }

    /**
     * Read Level 1 data from TGCA Return 3 DoubleMatrixDatasets, methylated
     * probe values, un-methylated probe values and beta values. Beta values are
     * calculated as : M/U+M (Non-bead studio / genome studio).
     *
     * @param fileInputFolder
     * @param printToFile
     * @param fileOutputFolder
     * @param TcgaMethod
     * @return all data out of the TCGA files Methylated / un-methylated and
     * beta-values
     */
    public static ArrayList<DoubleMatrixDataset<String, String>> parseTCGAData_lvl1_all_matrices(String fileInputFolder, boolean printToFile, String fileOutputFolder, boolean TcgaMethod) {

        File file = new File(fileInputFolder);
        File[] files = file.listFiles();
        ArrayList<File> vecFiles = new ArrayList<File>();
        for (int f = 0; f < files.length; f++) {
            if (files[f].getAbsolutePath().endsWith(".txt")) {
                vecFiles.add(files[f]);
            }
        }

        System.out.println("Files to parse:\t" + vecFiles.size());
        int nrSamples = vecFiles.size();

        int nrProbes = 0;
        ArrayList<String> vecProbes = new ArrayList<String>();
        try {
            BufferedReader in = new BufferedReader(new InputStreamReader(new FileInputStream(vecFiles.get(0)), ENCODING), 8096);
            String str;
            in.readLine();
            in.readLine();
            while ((str = in.readLine()) != null) {
                String[] data = SPLIT_ON_TAB.split(str);
                vecProbes.add(data[0]);
                nrProbes++;
            }
            in.close();
        } catch (IOException e) {
            System.out.println(e.getMessage());
            System.exit(-1);
        }

        System.out.println(nrProbes);

        DoubleMatrixDataset<String, String> dataset1 = new DoubleMatrixDataset<String, String>(nrProbes, nrSamples);
        DoubleMatrixDataset<String, String> dataset2 = new DoubleMatrixDataset<String, String>(nrProbes, nrSamples);
        DoubleMatrixDataset<String, String> dataset3 = new DoubleMatrixDataset<String, String>(nrProbes, nrSamples);
        //ExpressionDataset dataset3 = new ExpressionDataset(nrProbes * 2, nrSamples);
        for (int p = 0; p < vecProbes.size(); p++) {
            dataset1.rowObjects.set(p, vecProbes.get(p));
            dataset2.rowObjects.set(p, vecProbes.get(p));
            dataset3.rowObjects.set(p, vecProbes.get(p));
            //dataset3.probeNames[p * 2] = "M-" + (String) vecProbes.get(p);
            //dataset3.probeNames[p * 2 + 1] = "U-" + (String) vecProbes.get(p);
        }
        for (int f = 0; f < nrSamples; f++) {
            File currentFile = vecFiles.get(f);
            System.out.println("Processing:\t" + f + "\t" + currentFile.getAbsolutePath());
            try {
                BufferedReader in = new BufferedReader(new InputStreamReader(new FileInputStream(currentFile), ENCODING), 8096);
                String str = in.readLine();
                String[] data = SPLIT_ON_TAB.split(str);
                dataset1.colObjects.set(f, data[1]);
                dataset2.colObjects.set(f, data[1]);
                dataset3.colObjects.set(f, data[1]);
                str = in.readLine();
                data = SPLIT_ON_TAB.split(str);
                int columnM = -1;
                int columnU = -1;
                for (int d = 0; d < data.length; d++) {
                    if (data[d].toLowerCase().trim().replace(" ", "_").equals("methylated_signal_intensity_(m)")) {
                        columnM = d;
                    }
                    if (data[d].toLowerCase().trim().replace(" ", "_").equals("un-methylated_signal_intensity_(u)")) {
                        columnU = d;
                    }
                }
                //System.out.println(columnM + "\t" + columnU);
                int p = 0;
                while ((str = in.readLine()) != null) {
                    data = SPLIT_ON_TAB.split(str);
                    if (data[columnM].equals("NA") || data[columnM].equals("NaN")) {
                        data[columnM] = "-999";
                    }
                    if (data[columnU].equals("NA") || data[columnU].equals("NaN")) {
                        data[columnU] = "-999";
                    }
                    dataset1.rawData[p][f] = Double.parseDouble(data[columnM]);
                    dataset2.rawData[p][f] = Double.parseDouble(data[columnU]);

                    //Non bead-studio / genome-studio beta calculations. TCGA way.
                    if (dataset1.rawData[p][f] == -999 || dataset2.rawData[p][f] == -999) {
                        dataset3.rawData[p][f] = -999;
                    } else if(dataset1.rawData[p][f] == 0 && dataset2.rawData[p][f] == 0){
                        dataset3.rawData[p][f] = 0;
                    } else if (dataset1.rawData[p][f] <= 0 || dataset2.rawData[p][f] <= 0) {
                        dataset3.rawData[p][f] = -999;
                    } else if (TcgaMethod) {
                        dataset3.rawData[p][f] = dataset1.rawData[p][f] / (dataset1.rawData[p][f] + dataset2.rawData[p][f]);
                    } else {
                        dataset3.rawData[p][f] = dataset1.rawData[p][f] / ((dataset1.rawData[p][f] + dataset2.rawData[p][f]) + 100);
                    }
                    //dataset3.rawData[p * 2][f] = Double.parseDouble(data[columnM]);
                    //dataset3.rawData[p * 2 + 1][f] = Double.parseDouble(data[columnU]);

                    String probe = (String) vecProbes.get(p);
                    if (!data[0].equals(probe)) {
                        System.out.println("Error!:\t" + f + "\t" + data[0] + "\t" + probe);
                    }
                    p++;
                }
                in.close();
            } catch (IOException e) {
                System.out.println(e.getMessage());
                System.exit(-1);
            }

        }
        dataset1.recalculateHashMaps();
        dataset2.recalculateHashMaps();
        dataset3.recalculateHashMaps();

        if (printToFile) {
            try {
                dataset1.save(fileOutputFolder + "/TCGADataM.txt");
            } catch (IOException ex) {
                Logger.getLogger(ParseTcgaFile.class.getName()).log(Level.SEVERE, null, ex);
            }
            try {
                dataset2.save(fileOutputFolder + "/TCGADataU.txt");
            } catch (IOException ex) {
                Logger.getLogger(ParseTcgaFile.class.getName()).log(Level.SEVERE, null, ex);
            }
            try {
                dataset3.save(fileOutputFolder + "/TCGADataBeta.txt");
            } catch (IOException ex) {
                Logger.getLogger(ParseTcgaFile.class.getName()).log(Level.SEVERE, null, ex);
            }
        }
        ArrayList<DoubleMatrixDataset<String, String>> tmp = new ArrayList<DoubleMatrixDataset<String, String>>(3);

        tmp.add(dataset1);
        tmp.add(dataset2);
        tmp.add(dataset3);

        return (tmp);
    }

    /**
     * Read level 3 data from TCGA Return DoubleMatrixDataset beta-values
     *
     * @param fileInputFolder
     * @param nrProbes
     * @param printToFile
     * @param fileOut
     * @return
     */
    public static DoubleMatrixDataset<String, String> parseTCGAData_lvl3(String fileInputFolder, int nrProbes, boolean printToFile, String fileOut) {

        File file = new File(fileInputFolder);
        File[] files = file.listFiles();

        ArrayList<File> vecFiles = new ArrayList<File>();
        for (int f = 0; f < files.length; f++) {
            if (files[f].getAbsolutePath().endsWith(".txt")) {
                vecFiles.add(files[f]);
            }
        }


        System.out.println("Files to parse:\t" + vecFiles.size());
        int nrSamples = vecFiles.size();

        int nrP = 0;
        ArrayList<String> vecProbes = new ArrayList<String>();
        try {
            BufferedReader in = new BufferedReader(new InputStreamReader(new FileInputStream(vecFiles.get(0)), ENCODING), 8096);
            String str;
            in.readLine();
            in.readLine();
            while ((str = in.readLine()) != null) {
                String[] data = str.split("\t");
                vecProbes.add(data[0]);
                nrP++;
            }
            in.close();
        } catch (IOException e) {
            System.out.println(e.getMessage());
            System.exit(-1);
        }

        //
        DoubleMatrixDataset<String, String> dataset1 = new DoubleMatrixDataset<String, String>(nrProbes, nrSamples);

        //ExpressionDataset dataset3 = new ExpressionDataset(nrProbes * 2, nrSamples);
        for (int p = 0; p < vecProbes.size(); p++) {
            dataset1.rowObjects.set(p, (String) vecProbes.get(p));
            //dataset3.probeNames[p * 2] = "M-" + (String) vecProbes.get(p);
            //dataset3.probeNames[p * 2 + 1] = "U-" + (String) vecProbes.get(p);
        }
        for (int f = 0; f < nrSamples; f++) {
            File currentFile = vecFiles.get(f);
            System.out.println("Processing:\t" + f + "\t" + currentFile.getAbsolutePath());
            try {
                BufferedReader in = new BufferedReader(new InputStreamReader(new FileInputStream(currentFile), ENCODING), 8096);
                String str = in.readLine();
                String[] data = str.split("\t");
                dataset1.colObjects.set(f, data[1]);
                in.readLine();

                int p = 0;
                while ((str = in.readLine()) != null) {
                    data = str.split("\t");
                    if (f == 0) {
                        vecProbes.add(data[0]);
                    }

                    if (data[1].equals("NA")) {
                        data[1] = "-999";
                    }
                    dataset1.rawData[p][f] = Double.parseDouble(data[1]);

                    String probe = (String) vecProbes.get(p);
                    if (!data[0].equals(probe)) {
                        System.out.println("Error!:\t" + f + "\t" + data[0] + "\t" + probe);
                    }
                    p++;
                }
                in.close();
            } catch (IOException e) {
                System.out.println(e.getMessage());
                System.exit(-1);
            }

        }

        dataset1.recalculateHashMaps();

        if (printToFile) {
            try {
                dataset1.save(fileOut);
            } catch (IOException ex) {
                Logger.getLogger(ParseTcgaFile.class.getName()).log(Level.SEVERE, null, ex);
            }
        }
        System.out.println(dataset1.colObjects.toString());
        return (dataset1);
    }
}
