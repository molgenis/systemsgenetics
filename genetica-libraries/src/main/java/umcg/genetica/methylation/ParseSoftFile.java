/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package umcg.genetica.methylation;

import java.io.*;
import java.util.*;
import java.util.regex.Pattern;
import umcg.genetica.math.matrix.DoubleMatrixDataset;
import static umcg.genetica.methylation.ParseTcgaFile.ENCODING;

/**
 *
 * @author Lude & Marc Jan
 */
public class ParseSoftFile {

    private static Pattern SPLIT_ON_TAB = Pattern.compile("\\t");
    private static Pattern SPLIT_ON_EQUALS = Pattern.compile(" = ");

    public static HashMap<String, SoftfileAnnotation> importAnnotationFromSOFTFile(String fileLocation) throws Exception {

        HashMap<String, SoftfileAnnotation> dataAnnotation = new HashMap<String, SoftfileAnnotation>();

        if (fileLocation.equals("") || !fileLocation.endsWith(".soft")) {
            throw new Exception("No (correct) file specified");
        }

        try {
            BufferedReader in = new BufferedReader(new InputStreamReader(new FileInputStream(new File(fileLocation)), ENCODING), 8096);
            String str = "";
            while ((str = in.readLine()) != null) {
                if (str.startsWith("^SAMPLE")) {
                    String sampleName = SPLIT_ON_EQUALS.split(str)[1];
                    SoftfileAnnotation tmp = new SoftfileAnnotation();

                    while ((str = in.readLine()) != null) {

                        if (str.startsWith("^SAMPLE")) {
                            tmp.putAnnotationInformation("Error_In_SoftFile", "True");
                            dataAnnotation.put(sampleName, tmp);
                            sampleName = SPLIT_ON_EQUALS.split(str)[1];
                            tmp = new SoftfileAnnotation();
                        }

                        if (str.startsWith("!Sample_title = ")) {
                            tmp.setTitle(SPLIT_ON_EQUALS.split(str)[1]);
                        } else if (str.startsWith("!Sample_geo_accession = ")) {
                            tmp.setAccession(SPLIT_ON_EQUALS.split(str)[1]);
                        }

                        if (str.startsWith("!Sample_") && str.contains(" = ")) {
                            String[] stringSplit = SPLIT_ON_EQUALS.split(str);
                            if (stringSplit.length == 2) {
                                if (tmp.getAnnotationInformation().containsKey(SPLIT_ON_EQUALS.split(str)[0])) {
                                    tmp.putAnnotationInformation(SPLIT_ON_EQUALS.split(str)[0], tmp.getAnnotationInformation().get(SPLIT_ON_EQUALS.split(str)[0]) + " // " + SPLIT_ON_EQUALS.split(str)[1]);
                                } else {
                                    tmp.putAnnotationInformation(SPLIT_ON_EQUALS.split(str)[0], SPLIT_ON_EQUALS.split(str)[1]);
                                }
                            } else {
                                //System.out.println(str);
                            }

                        } else if (str.startsWith("#") && str.contains(" = ")) {
                            String[] stringSplit = SPLIT_ON_EQUALS.split(str);
                            if (stringSplit.length == 2) {
                                if (tmp.getAnnotationInformation().containsKey(SPLIT_ON_EQUALS.split(str)[0])) {
                                    tmp.putAnnotationInformation(SPLIT_ON_EQUALS.split(str)[0], tmp.getAnnotationInformation().get(SPLIT_ON_EQUALS.split(str)[0]) + " " + SPLIT_ON_EQUALS.split(str)[1]);
                                } else {
                                    tmp.putAnnotationInformation(SPLIT_ON_EQUALS.split(str)[0], SPLIT_ON_EQUALS.split(str)[1]);
                                }
                            } else {
                                //System.out.println(str);
                            }
                        }

                        if (str.startsWith("!sample_table_begin")) {
                            dataAnnotation.put(sampleName, tmp);
                            break;
                        }
                    }
                }
            }
            in.close();
        } catch (IOException e) {
            System.out.println(e.getMessage());
            System.exit(-1);
        }

        return (dataAnnotation);
    }

    public static DoubleMatrixDataset<String, String> importSOFTFile(String fileLocation, int numberOfSamples, int numberOfProbes) throws Exception {

        if (fileLocation.equals("") || !fileLocation.endsWith(".soft")) {
            throw new Exception("No (correct) file specified");
        }

        boolean debug = false;

        String probeHeader = "ID_REF";
        String valueHeader = "VALUE";
        String intensityHeader = "Intensity";

        DoubleMatrixDataset<String, String> dataset = new DoubleMatrixDataset<String, String>(numberOfProbes, numberOfSamples);

        HashMap<String, Integer> hashUniqueProbes = new HashMap<String, Integer>();
        ArrayList<String> uniqueProbes = new ArrayList<String>();
        ArrayList<String> uniqueSamples = new ArrayList<String>();

        int sampleID = 0;

        try {
            BufferedReader in = new BufferedReader(new InputStreamReader(new FileInputStream(new File(fileLocation)), ENCODING), 8096);
            String str = "";
            while ((str = in.readLine()) != null) {
                if (str.startsWith("^")) {

                    if (debug) {
                        System.out.println(str);
                    }

                    if (str.startsWith("^SAMPLE")) {
                        String sampleName = SPLIT_ON_EQUALS.split(str)[1];
                        uniqueSamples.add(str);
                        while ((str = in.readLine()) != null) {
                            //System.out.println(str);
                            int nrProbesThisSample = 0;
                            if (str.startsWith("!Sample_supplementary_file") && debug) {
                                System.out.println(str);
                            }
                            if (str.startsWith("!Sample_characteristics_ch1") && debug) {
                                if (str.toLowerCase().contains("male") || str.toLowerCase().contains("female")) {
                                    System.out.println(sampleName + "\t" + str);
                                }
                            }
                            if (str.startsWith("!sample_table_begin")) {
                                str = in.readLine();

                                if (debug) {
                                    System.out.println(str);
                                }

                                String[] data = SPLIT_ON_TAB.split(str);

                                int probeHeaderColumn = -1;
                                int valueHeaderColumn = -1;
                                int intensityHeaderColumn = -1;

                                double[] valsValue = new double[numberOfProbes];

                                for (int d = 0; d < data.length; ++d) {
                                    if (data[d].equals(probeHeader)) {
                                        probeHeaderColumn = d;
                                    }
                                    if (data[d].equals(valueHeader)) {
                                        valueHeaderColumn = d;
                                    }
                                    if (data[d].equals(intensityHeader)) {
                                        intensityHeaderColumn = d;
                                    }
                                }

                                if (intensityHeaderColumn != -1) {
                                    //valueHeaderColumn = intensityHeaderColumn;
                                }

                                int nrMissingValues = 0;

                                while ((str = in.readLine()) != null) {

                                    if (str.startsWith("!sample_table_end")) {
                                        break;
                                    }

                                    if (valueHeaderColumn != -1) {
                                        data = SPLIT_ON_TAB.split(str);
                                        double value = 0;
                                        if (data.length <= valueHeaderColumn || data[valueHeaderColumn] == null || data[valueHeaderColumn].length() == 0 || data[valueHeaderColumn].equalsIgnoreCase("null")) {
                                            value = -999;
                                            nrMissingValues++;
                                        } else {
                                            value = Double.parseDouble(data[valueHeaderColumn]);
                                            nrProbesThisSample++;
                                        }

                                        if (!hashUniqueProbes.containsKey(data[probeHeaderColumn])) {
                                            int probeID = hashUniqueProbes.size();
                                            hashUniqueProbes.put(data[probeHeaderColumn], probeID);
                                            uniqueProbes.add(data[probeHeaderColumn]);
                                            valsValue[probeID] = value;
                                        } else {
                                            int probeID = ((Integer) hashUniqueProbes.get(data[probeHeaderColumn])).intValue();
                                            valsValue[probeID] = value;
                                        }
                                    }
                                }
                                if (probeHeaderColumn != -1) {
                                    for (int p = 0; p < valsValue.length; p++) {
                                        dataset.rawData[p][sampleID] = valsValue[p];
                                    }
                                    dataset.colObjects.set(sampleID, sampleName);
                                    if (debug) {
                                        System.out.println(sampleName + "\t" + sampleID + "\tNrProbesThisSample:\t" + nrProbesThisSample + "\tNrMissingProbeValues:\t" + nrMissingValues + "\t" + hashUniqueProbes.size() + "\t" + uniqueSamples.size());
                                    }
                                    sampleID++;
                                }
                                break;
                            }
                        }
                    }
                }
            }
            in.close();
        } catch (IOException e) {
            System.out.println(e.getMessage());
            System.exit(-1);
        }

        if (debug) {
            System.out.println("Total number of samples:\t" + sampleID);
        }

        for (int r = 0; r < dataset.nrRows; r++) {
            dataset.rowObjects.set(r, (String) uniqueProbes.get(r));
        }

        dataset.recalculateHashMaps();
        return (dataset);
    }

    public static DoubleMatrixDataset<String, String> importSOFTFileSelection(String fileLocation, int numberOfSamples, int numberOfProbes, int samplesPerDataset) throws Exception {

        if (fileLocation.equals("") || !fileLocation.endsWith(".soft")) {
            throw new Exception("No (correct) file specified");
        }

        boolean debug = false;

        String probeHeader = "ID_REF";
        String valueHeader = "VALUE";
        String intensityHeader = "Intensity";

        HashMap<String, Integer> SelectionBasedOnSeriesId = new HashMap<String, Integer>();

        DoubleMatrixDataset<String, String> dataset = new DoubleMatrixDataset<String, String>(numberOfProbes, numberOfSamples);

        HashMap<String, Integer> hashUniqueProbes = new HashMap<String, Integer>();
        ArrayList<String> uniqueProbes = new ArrayList<String>();
        ArrayList<String> uniqueSamples = new ArrayList<String>();

        int sampleID = 0;

        try {
            BufferedReader in = new BufferedReader(new InputStreamReader(new FileInputStream(new File(fileLocation)), ENCODING), 8096);
            String str = "";
            while ((str = in.readLine()) != null) {
                if (str.startsWith("^")) {

                    if (debug) {
                        System.out.println(str);
                    }

                    if (str.startsWith("^SAMPLE")) {
                        String sampleName = SPLIT_ON_EQUALS.split(str)[1];
                        String sampleSeriesId = "";

                        while ((str = in.readLine()) != null) {
                            //System.out.println(str);
                            int nrProbesThisSample = 0;

                            if (str.startsWith("!Sample_series_id")) {
                                if (debug) {
                                    System.out.println(str);
                                }
                                sampleSeriesId = SPLIT_ON_EQUALS.split(str)[1];
                            }

                            if (str.startsWith("!Sample_supplementary_file") && debug) {
                                System.out.println(str);
                            }

                            if (str.startsWith("!Sample_characteristics_ch1") && debug) {
                                if (str.toLowerCase().contains("male") || str.toLowerCase().contains("female")) {
                                    System.out.println(sampleName + "\t" + str);
                                }
                            }

                            if (str.startsWith("!sample_table_begin")) {
                                str = in.readLine();

                                if (debug) {
                                    System.out.println(str);
                                }

                                String[] data = SPLIT_ON_TAB.split(str);

                                int probeHeaderColumn = -1;
                                int valueHeaderColumn = -1;
                                int intensityHeaderColumn = -1;

                                double[] valsValue = new double[numberOfProbes];

                                for (int d = 0; d < data.length; ++d) {
                                    if (data[d].equals(probeHeader)) {
                                        probeHeaderColumn = d;
                                    }
                                    if (data[d].equals(valueHeader)) {
                                        valueHeaderColumn = d;
                                    }
                                    if (data[d].equals(intensityHeader)) {
                                        intensityHeaderColumn = d;
                                    }
                                }

                                if (intensityHeaderColumn != -1) {
                                    //valueHeaderColumn = intensityHeaderColumn;
                                }

                                int nrMissingValues = 0;

                                while ((str = in.readLine()) != null) {

                                    if (str.startsWith("!sample_table_end")) {
                                        break;
                                    }

                                    if (valueHeaderColumn != -1) {
                                        data = SPLIT_ON_TAB.split(str);
                                        double value = 0;
                                        if (data.length <= valueHeaderColumn || data[valueHeaderColumn] == null || data[valueHeaderColumn].length() == 0 || data[valueHeaderColumn].equalsIgnoreCase("null")) {
                                            value = -999;
                                            nrMissingValues++;
                                        } else {
                                            value = Double.parseDouble(data[valueHeaderColumn]);
                                            nrProbesThisSample++;
                                        }

                                        if (!hashUniqueProbes.containsKey(data[probeHeaderColumn])) {
                                            int probeID = hashUniqueProbes.size();
                                            hashUniqueProbes.put(data[probeHeaderColumn], probeID);
                                            uniqueProbes.add(data[probeHeaderColumn]);
                                            valsValue[probeID] = value;
                                        } else {
                                            int probeID = ((Integer) hashUniqueProbes.get(data[probeHeaderColumn])).intValue();
                                            valsValue[probeID] = value;
                                        }
                                    }
                                }

                                if (SelectionBasedOnSeriesId.containsKey(sampleSeriesId)) {
                                    SelectionBasedOnSeriesId.put(sampleSeriesId, SelectionBasedOnSeriesId.get(sampleSeriesId) + 1);
                                } else {
                                    SelectionBasedOnSeriesId.put(sampleSeriesId, 1);
                                }

                                if (probeHeaderColumn != -1 && !(SelectionBasedOnSeriesId.get(sampleSeriesId) > samplesPerDataset)) {
                                    if (debug) {
                                        System.out.println(sampleSeriesId + "\t" + SelectionBasedOnSeriesId.get(sampleSeriesId));
                                    }

                                    uniqueSamples.add(sampleName);

                                    for (int p = 0; p < valsValue.length; p++) {
                                        dataset.rawData[p][sampleID] = valsValue[p];
                                    }

                                    dataset.colObjects.set(sampleID, sampleName);
                                    if (debug) {
                                        System.out.println(sampleName + "\t" + sampleID + "\tNrProbesThisSample:\t" + nrProbesThisSample + "\tNrMissingProbeValues:\t" + nrMissingValues + "\t" + hashUniqueProbes.size() + "\t" + uniqueSamples.size());
                                    }
                                    sampleID++;
                                }
                                break;
                            }
                        }
                    }
                }
            }
            in.close();
        } catch (IOException e) {
            System.out.println(e.getMessage());
            System.exit(-1);
        }

        if (debug) {
            System.out.println("Total number of samples:\t" + sampleID);
        }

        for (int r = 0; r < dataset.nrRows; r++) {
            dataset.rowObjects.set(r, (String) uniqueProbes.get(r));
        }


        dataset.recalculateHashMaps();
        return (dataset);
    }
}
