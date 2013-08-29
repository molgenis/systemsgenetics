///*
// * To change this template, choose Tools | Templates
// * and open the template in the editor.
// */
//
//package umcg.genetica.legacy;
//
//import java.io.*;
//import java.util.*;
//import java.awt.image.*;
//import java.awt.*;
//import java.rmi.Remote;
//import java.util.regex.Pattern;
//import umcg.genetica.math.stats.*;
//import umcg.genetica.methylation.DeepCopy;
//
///**
// *
// * @author lude
// */
//public class ParseSOFTFile_lude {
//	
//	
//	private static Pattern SPLIT_ON_TAB = Pattern.compile("\\t");
//	private static Pattern SPLIT_ON_EQUALS = Pattern.compile(" = ");
//
//    public ParseSOFTFile_lude(){
//		
//	}
//	
//	public ParseSOFTFile_lude(String fileLocation, int numberOfSamples, int numberOfProbes) throws Exception {
//		
//		importSOFTFile(fileLocation, numberOfSamples, numberOfProbes);
//		
//		//performQC();
//
//        //performPCA();
//
//        //convertPCsToGeneLevel();
//
//        //parseDataPeterBram();
//
//        /*
//        performCoexpressionOnPCScoresFromPCAOverSamples();
//        
//        if (1==2) {
//            ExpressionDataset dataset = new ExpressionDataset("/Users/lude/Downloads/PCATCGABeta/PrincipalcomponentsTCGA+GPL8490Combined.txt");
//            dataset.transposeDataset();
//            for (int p=0; p<100; p++) {
//                for (int q=100; q<200; q++) {
//                    double corr = JSci.maths.ArrayMath.correlation(dataset.rawData[p], dataset.rawData[q]);
//                    double r2  = corr * corr;
//                    if (r2 > 0.35) {
//                        System.out.println(p + "\t" + q + "\t" + corr + "\t" + r2);
//                        for (int  s=0; s<dataset.nrSamples; s++) {
//                            System.out.println(dataset.rawData[p][s] + "\t" + dataset.rawData[q][s]);
//                        }
//                    }
//                }
//            }
//            System.exit(0);
//        }
//         * 
//         */
//
//        //comparePCScoresWithGPL570TCs();
//        
//        //compareSexPhenotype();
//        
//        //parseTCGAData();
//
//        //compareTCGAPCsWithFeatures();
//        
//        //compareTCGAEVsWithBatches();
//        
//        //compareTCGAWithGPL8490();
//        
//        System.exit(0);
//
//    }
//	
//	public ExpressionDataset_lude importSOFTFile(String fileLocation, int numberOfSamples, int numberOfProbes) throws Exception {
//		
//		if(fileLocation.equals("") || !fileLocation.endsWith(".soft")){
//			throw new Exception("No (correct) file specified");
//		}
//		
//		boolean debug=false;
//		
//        String probeHeader = "ID_REF";
//        String valueHeader = "VALUE";
//        String intensityHeader = "Intensity";
//
//        ExpressionDataset_lude dataset = new ExpressionDataset_lude(numberOfProbes, numberOfSamples);
//
//        HashMap<String, Integer> hashUniqueProbes = new HashMap<String, Integer>();
//        ArrayList<String> vecUniqueProbes = new ArrayList<String>();
//        ArrayList<String> vecUniqueSamples = new ArrayList<String>();
//		
//        int sampleID = 0;
//		
//        try {
//			java.io.BufferedReader in = new java.io.BufferedReader(new java.io.FileReader(new File(fileLocation)));
//            String str = "";
//            while ((str = in.readLine()) != null) {
//                if (str.startsWith("^")) {
//					
//                    if(debug) System.out.println(str);
//					
//                    if (str.startsWith("^SAMPLE")) {
//                        String sampleName = SPLIT_ON_EQUALS.split(str)[1];
//                        vecUniqueSamples.add(str);
//                        while ((str = in.readLine()) != null) {
//                            //System.out.println(str);
//                            int nrProbesThisSample = 0;
//                            if (str.startsWith("!Sample_supplementary_file") && debug) {
//                                System.out.println(str);
//                            }
//                            if (str.startsWith("!Sample_characteristics_ch1") && debug) {
//                                if (str.toLowerCase().contains("male") || str.toLowerCase().contains("female") ) {
//                                    System.out.println(sampleName + "\t" + str);
//                                }
//                            }
//                            if (str.startsWith("!sample_table_begin")) {
//                                str = in.readLine();
//								
//                                if(debug) System.out.println(str);
//								
//                                String[] data = SPLIT_ON_TAB.split(str);
//								
//                                int probeHeaderColumn = -1;
//                                int valueHeaderColumn = -1;
//                                int intensityHeaderColumn = -1;
//								
//                                double[] valsValue = new double[numberOfProbes];
//								
//                                for (int d=0; d<data.length; ++d) {
//                                    if (data[d].equals(probeHeader)) {
//                                        probeHeaderColumn = d;
//                                    }
//                                    if (data[d].equals(valueHeader)) {
//                                        valueHeaderColumn = d;
//                                    }
//                                    if (data[d].equals(intensityHeader)) {
//                                        intensityHeaderColumn = d;
//                                    }
//                                }
//								
//                                if (intensityHeaderColumn!=-1) {
//                                    //valueHeaderColumn = intensityHeaderColumn;
//                                }
//								
//                                int nrMissingValues = 0;
//								
//                                while ((str = in.readLine()) != null) {
//									
//                                    if (str.startsWith("!sample_table_end")){
//										break;
//									}
//									
//                                    if (valueHeaderColumn!=-1){
//                                        data = SPLIT_ON_TAB.split(str);
//                                        double value = 0;
//                                        if (data.length <= valueHeaderColumn || data[valueHeaderColumn]==null || data[valueHeaderColumn].length()==0 || data[valueHeaderColumn].equalsIgnoreCase("null")) {
//                                            value = -999;
//                                            nrMissingValues++;
//                                        } else {
//                                            value = Double.parseDouble(data[valueHeaderColumn]);
//                                            nrProbesThisSample++;
//                                        }
//										
//                                        if (!hashUniqueProbes.containsKey(data[probeHeaderColumn])) {
//                                            int probeID = hashUniqueProbes.size();
//                                            hashUniqueProbes.put(data[probeHeaderColumn], probeID);
//                                            vecUniqueProbes.add(data[probeHeaderColumn]);
//                                            valsValue[probeID] = value;
//                                        } else {
//                                            int probeID = ((Integer) hashUniqueProbes.get(data[probeHeaderColumn])).intValue();
//                                            valsValue[probeID] = value;
//                                        }
//                                    }
//                                }
//                                if (probeHeaderColumn!=-1) {
//									
//                                    for (int p=0; p<valsValue.length; p++) {
//                                        dataset.rawData[p][sampleID] = valsValue[p];
//                                    }
//									
//                                    dataset.sampleNames[sampleID] = sampleName;
//                                    if(debug) System.out.println(sampleName + "\t" + sampleID + "\tNrProbesThisSample:\t" + nrProbesThisSample + "\tNrMissingProbeValues:\t" + nrMissingValues + "\t" + hashUniqueProbes.size() + "\t" + vecUniqueSamples.size());
//                                    sampleID++;
//                                }
//                                break;
//                            }
//                        }
//                    }
//                }
//            }
//            in.close();
//        } catch (IOException e) {
//            e.printStackTrace();
//            System.out.println(e.getMessage());
//            System.exit(-1);
//        }
//		
//        if(debug) System.out.println("Total number of samples:\t" + sampleID);
//		
//        for (int p=0; p<dataset.nrProbes; p++) {
//            dataset.probeNames[p] = (String) vecUniqueProbes.get(p);
//        }
//		
//		dataset.recalculateHashMaps();
//		return(dataset);
//    }
//
//    public ExpressionDataset_lude importSOFTFileSelection(String fileLocation, int numberOfSamples, int numberOfProbes, int samplesPerDataset) throws Exception {
//		
//		if(fileLocation.equals("") || !fileLocation.endsWith(".soft")){
//			throw new Exception("No (correct) file specified");
//		}
//		
//		boolean debug = false;
//		
//        String probeHeader = "ID_REF";
//        String valueHeader = "VALUE";
//        String intensityHeader = "Intensity";
//		
//		HashMap<String, Integer> SelectionBasedOnSeriesId = new HashMap<String, Integer>();
//
//        ExpressionDataset_lude dataset = new ExpressionDataset_lude(numberOfProbes, numberOfSamples);
//
//        HashMap<String, Integer> hashUniqueProbes = new HashMap<String, Integer>();
//        ArrayList<String> uniqueProbes = new ArrayList<String>();
//        ArrayList<String> uniqueSamples = new ArrayList<String>();
//		
//        int sampleID = 0;
//		
//        try {
//			java.io.BufferedReader in = new java.io.BufferedReader(new java.io.FileReader(new File(fileLocation)));
//            String str = "";
//            while ((str = in.readLine()) != null) {
//                if (str.startsWith("^")) {
//					
//                    if(debug) System.out.println(str);
//					
//                    if (str.startsWith("^SAMPLE")) {
//                        String sampleName = SPLIT_ON_EQUALS.split(str)[1];
//						String sampleSeriesId = "";
//                        
//                        while ((str = in.readLine()) != null) {
//                            //System.out.println(str);
//                            int nrProbesThisSample = 0;
//							
//							if (str.startsWith("!Sample_series_id")) {
//								if(debug)System.out.println(str);
//								sampleSeriesId = SPLIT_ON_EQUALS.split(str)[1];
//                            }
//							
//                            if (str.startsWith("!Sample_supplementary_file") && debug) {
//                                System.out.println(str);
//                            }
//							
//                            if (str.startsWith("!Sample_characteristics_ch1") && debug) {
//                                if (str.toLowerCase().contains("male") || str.toLowerCase().contains("female") ) {
//                                    System.out.println(sampleName + "\t" + str);
//                                }
//                            }
//							
//                            if (str.startsWith("!sample_table_begin")) {
//                                str = in.readLine();
//								
//                                if(debug) System.out.println(str);
//								
//                                String[] data = SPLIT_ON_TAB.split(str);
//								
//                                int probeHeaderColumn = -1;
//                                int valueHeaderColumn = -1;
//                                int intensityHeaderColumn = -1;
//								
//                                double[] valsValue = new double[numberOfProbes];
//								
//                                for (int d=0; d<data.length; ++d) {
//                                    if (data[d].equals(probeHeader)) {
//                                        probeHeaderColumn = d;
//                                    }
//                                    if (data[d].equals(valueHeader)) {
//                                        valueHeaderColumn = d;
//                                    }
//                                    if (data[d].equals(intensityHeader)) {
//                                        intensityHeaderColumn = d;
//                                    }
//                                }
//								
//                                if (intensityHeaderColumn!=-1) {
//                                    //valueHeaderColumn = intensityHeaderColumn;
//                                }
//								
//                                int nrMissingValues = 0;
//								
//                                while ((str = in.readLine()) != null) {
//									
//                                    if (str.startsWith("!sample_table_end")){
//										break;
//									}
//									
//                                    if (valueHeaderColumn!=-1){
//                                        data = SPLIT_ON_TAB.split(str);
//                                        double value = 0;
//                                        if (data.length <= valueHeaderColumn || data[valueHeaderColumn]==null || data[valueHeaderColumn].length()==0 || data[valueHeaderColumn].equalsIgnoreCase("null")) {
//                                            value = -999;
//                                            nrMissingValues++;
//                                        } else {
//                                            value = Double.parseDouble(data[valueHeaderColumn]);
//                                            nrProbesThisSample++;
//                                        }
//										
//                                        if (!hashUniqueProbes.containsKey(data[probeHeaderColumn])) {
//                                            int probeID = hashUniqueProbes.size();
//                                            hashUniqueProbes.put(data[probeHeaderColumn], probeID);
//                                            uniqueProbes.add(data[probeHeaderColumn]);
//                                            valsValue[probeID] = value;
//                                        } else {
//                                            int probeID = ((Integer) hashUniqueProbes.get(data[probeHeaderColumn])).intValue();
//                                            valsValue[probeID] = value;
//                                        }
//                                    }
//								}								
//								
//								if(SelectionBasedOnSeriesId.containsKey(sampleSeriesId)){
//									SelectionBasedOnSeriesId.put(sampleSeriesId, SelectionBasedOnSeriesId.get(sampleSeriesId)+1);
//								} else {
//									SelectionBasedOnSeriesId.put(sampleSeriesId, 1);
//								}
//								
//                                if (probeHeaderColumn!=-1 && !(SelectionBasedOnSeriesId.get(sampleSeriesId)>5)) {
//									if(debug) System.out.println(sampleSeriesId+"\t"+SelectionBasedOnSeriesId.get(sampleSeriesId));
//
//									uniqueSamples.add(sampleName);
//
//									for (int p=0; p<valsValue.length; p++) {
//										dataset.rawData[p][sampleID] = valsValue[p];
//									}
//
//									dataset.sampleNames[sampleID] = sampleName;
//									if(debug) System.out.println(sampleName + "\t" + sampleID + "\tNrProbesThisSample:\t" + nrProbesThisSample + "\tNrMissingProbeValues:\t" + nrMissingValues + "\t" + hashUniqueProbes.size() + "\t" + uniqueSamples.size());
//									sampleID++;
//								}
//                                break;
//                            }
//                        }
//                    }
//                }
//            }
//            in.close();
//        } catch (IOException e) {
//            e.printStackTrace();
//            System.out.println(e.getMessage());
//            System.exit(-1);
//        }
//		
//        if(debug) System.out.println("Total number of samples:\t" + sampleID);
//		
//        for (int p=0; p<dataset.nrProbes; p++) {
//            dataset.probeNames[p] = (String) uniqueProbes.get(p);
//        }
//		
//		dataset.recalculateHashMaps();
//		return(dataset);
//    }
//	
//	public ExpressionDataset_lude performQC(ExpressionDataset_lude dataset, int maxMissingSamplesMissingAProbe, int maxMissingProbesMissingPerSample) {
//		
//		boolean debug = true;
//		int originalNumberSamples = dataset.nrSamples;
//		int originalNumberProbes = dataset.nrProbes;
//        
//        String[] samplesToExplicitlyExlcude = {""};
//        HashMap hashSamplesToExplicitlyExlcude = new HashMap();
//		
//		if(samplesToExplicitlyExlcude.length>0){		
//			for (int s=0; s<samplesToExplicitlyExlcude.length; s++) {
//				hashSamplesToExplicitlyExlcude.put(samplesToExplicitlyExlcude[s], null);
//				hashSamplesToExplicitlyExlcude.put(samplesToExplicitlyExlcude[s].replace("GSM", ""), null);
//			}
//		}
//
//        HashMap<String, Boolean> hashSamplesPassingQC = new HashMap<String, Boolean>();
//		HashMap<String, Boolean> hashProbesPassingQC = new HashMap<String, Boolean>();
//		
//        for (int s=0; s<dataset.nrSamples; ++s) {
//            int nrMissing = 0;
//            for (int p=0; p<dataset.nrProbes; ++p) {
//                if (dataset.rawData[p][s]==-999){
//					nrMissing++;
//				}
//            }
//            if (nrMissing>=maxMissingProbesMissingPerSample) {
//                if(debug) System.out.println(s + "\t" + dataset.sampleNames[s] + "\t" + nrMissing);
//            } else {
//                if (!hashSamplesToExplicitlyExlcude.containsKey(dataset.sampleNames[s])) {
//                    hashSamplesPassingQC.put(dataset.sampleNames[s], true);
//                }
//            }
//        }
//		
//        dataset = dataset.RemoveProbesOrSamples(dataset, hashProbesPassingQC, hashSamplesPassingQC);
//		System.out.println("Number of samples removed: "+(originalNumberSamples-dataset.rawData[1].length));
//		
//        for (int p=0; p<dataset.nrProbes; ++p) {
//            int nrMissing = 0;
//            for (int s=0; s<dataset.nrSamples; ++s) {
//                if (dataset.rawData[p][s]==-999) nrMissing++;
//            }
//            if (nrMissing>=maxMissingSamplesMissingAProbe) {
//				if(debug) System.out.println("Excluding probe:\t" + p + "\t" + dataset.probeNames[p] + "\t" + nrMissing);
//            } else {
//                hashProbesPassingQC.put(dataset.probeNames[p], null);
//            }
//        }
//		
//        dataset = dataset.RemoveProbesOrSamples(dataset, hashProbesPassingQC, hashSamplesPassingQC);
//		System.out.println("Number of probes removed: "+(originalNumberProbes-dataset.rawData.length));
//		
//		QuantileNormalization.QuantileNormAdressingNaValuesBeforeQN(dataset.rawData, true, false);
//		
//		
//		//OLD:
//        int[] nanPerProbe = new int[dataset.nrProbes];
//        int[] nanPerSample = new int[dataset.nrSamples];
//        for (int p=0; p<dataset.nrProbes; p++) {
//            for (int s=0; s<dataset.nrSamples; s++) {
//                if (Double.isNaN(dataset.rawData[p][s]) || Double.isInfinite(dataset.rawData[p][s]) || dataset.rawData[p][s]==-999) {
//                    nanPerProbe[p]++;
//                    nanPerSample[s]++;
//                }
//            }
//        }
//
//		ExpressionDataset_lude datasetSorted = new ExpressionDataset_lude(dataset.nrSamples, dataset.nrProbes);
//        //Quantile normalisation, allowing for missing values:
//		
//        for (int s=0; s<dataset.nrSamples; s++) {
//            double[] vals = new double[dataset.nrProbes];
//            for (int p=0; p<dataset.nrProbes; p++) {
//                if (dataset.rawData[p][s]==-999) {
//					//Shouldn't this be the median?
//                    vals[p] = 0;
//                } else {
//                    vals[p] = dataset.rawData[p][s];
//                }
//            }
//            Arrays.sort(vals);
//            datasetSorted.rawData[s] = vals;
//        }
//        
//		datasetSorted.transposeDataset();
//		
//        double[] dist = new double[dataset.nrProbes];
//        for (int p=0; p<dataset.nrProbes; p++) {
//            dist[p] = JSci.maths.ArrayMath.median(datasetSorted.rawData[p]);
//            if(debug) System.out.println(p + "\t" + dist[p] + "\t" + JSci.maths.ArrayMath.mean(datasetSorted.rawData[p]));
//        }
//		datasetSorted = null;
//		
//        Arrays.sort(dist);
//
//        //Perform quantile normalization
//		//Original Lude
//		
//        for (int s=0; s<dataset.nrSamples; s++) {
//            ArrayList<Double> vec1 = new ArrayList<Double>();
//            for (int p=0; p<dataset.nrProbes; p++) {
//                if (dataset.rawData[p][s]!=-999) {
//                    vec1.add(dataset.rawData[p][s]);
//                }
//            }
//            double[] vals1 = new double[vec1.size()];
//            for (int v = 0; v < vals1.length; v++) {
//                vals1[v] = ((Double) vec1.get(v)).doubleValue();
//            }
//            umcg.genetica.util.Rank rank = new umcg.genetica.util.Rank(vals1, 0d);
//            double[] valsRanked = rank.getRanks();
//            for (int v=0; v<vals1.length; v++){
//				valsRanked[v]--;
//			}
//            for (int v=0; v<vals1.length; v++) {
//                double quantile = ((double) valsRanked[v]) / ((double) vals1.length + 1d);
//                int distIndex = (int) (quantile * (double) dataset.nrProbes);
//                vals1[v] = dist[distIndex];
//            }
//            int itr = 0;
//            for (int p=0; p<dataset.nrProbes; p++) {
//                if (dataset.rawData[p][s]!=-999) {
//                    dataset.rawData[p][s] = vals1[itr];
//                    itr++;
//                }
//            }
//        }
//		
//        //Replace missing values:
//        for (int p=0; p<dataset.nrProbes; p++) {
//            double valSum = 0;
//            int nr = 0;
//            for (int s=0; s<dataset.nrSamples; s++) {
//                if (dataset.rawData[p][s]!=-999) {
//                    valSum += dataset.rawData[p][s];
//                    nr++;
//                }
//            }
//            double mean = valSum / nr;
//            for (int s=0; s<dataset.nrSamples; s++) {
//                if (dataset.rawData[p][s]==-999) {
//                    dataset.rawData[p][s] = mean;
//                }
//            }
//        }
//		
//		return(dataset);
//    }
//
//	public void compareTCGAWithGPL8490() {
//
//        int nrTCs = 100;
//        HashMap hashTCsToIncludeHuman = new HashMap();
//        for (int tc=0; tc<nrTCs; tc++) {
//            hashTCsToIncludeHuman.put("Comp" + String.valueOf(tc + 1), null);
//        }
//        
//        ExpressionDataset_lude dataset1 = new ExpressionDataset_lude("/Users/lude/Documents/DMG/Data/GPL8490/PCAOver4432SamplesProbesCenteredScaled/principalcomponents.txt", "\t", null, hashTCsToIncludeHuman);
//        ExpressionDataset_lude dataset2 = new ExpressionDataset_lude("/Users/lude/Downloads/PCATCGABeta/principalcomponents.txt", "\t", dataset1.hashProbes, hashTCsToIncludeHuman);
//        dataset1 = new ExpressionDataset_lude(dataset1.fileName, "\t", dataset2.hashProbes, hashTCsToIncludeHuman);
//        ExpressionDataset_lude dataset3 = new ExpressionDataset_lude(dataset2.nrProbes, dataset2.nrSamples);
//        for (int tc=0; tc<dataset2.nrSamples; tc++) {
//            dataset3.sampleNames[tc] = dataset2.sampleNames[tc];
//        }
//        for (int p=0; p<dataset2.nrProbes; p++) {
//            int probeIndex = ((Integer) dataset1.hashProbes.get(dataset2.probeNames[p])).intValue();
//            for (int tc=0; tc<dataset2.nrSamples; tc++) {
//                dataset3.rawData[probeIndex][tc] = dataset2.rawData[p][tc];
//            }
//            dataset3.probeNames[probeIndex] = dataset1.probeNames[p];
//        }
//        dataset3.save(dataset2.fileName + ".ProbeOrderingAsGPL8490.txt");
//        dataset1.transposeDataset();
//        dataset2 = null;
//        dataset3.transposeDataset();
//        for (int tc1=0; tc1<dataset1.nrProbes; tc1++) {
//            double r2Sum = 0;
//            for (int tc2=0; tc2<dataset3.nrProbes; tc2++) {
//                double correlation = JSci.maths.ArrayMath.correlation(dataset1.rawData[tc1], dataset3.rawData[tc2]);
//                double r2 = correlation * correlation; 
//                if (r2 > 0.01) {
//                    System.out.println(tc1 + "\t" + tc2 + "\t" + r2);
//                }
//                r2Sum+=r2;
//            }
//            System.out.println(tc1 + "\t" + r2Sum);
//        }
//        
//        System.exit(0);
//    }
//    
//    public void parseTCGAData() {
//        
//        String fileDirTCGAFiles = "/Users/lude/Downloads/774bde41-d002-447f-805b-0da8d94b30fa/DNA_Methylation/JHU_USC__HumanMethylation27/Level_1/";
//        File file = new File(fileDirTCGAFiles);
//        File[] files = file.listFiles();
//        Vector vecFiles = new Vector();
//        for (int f=0; f<files.length; f++) {
//            if (files[f].getAbsolutePath().endsWith(".txt")) {
//                vecFiles.add(files[f]);
//            }
//        }
//        System.out.println("Files to parse:\t" + vecFiles.size());
//        int nrSamples = vecFiles.size();
//        
//        int nrProbes = 0;
//        Vector vecProbes = new Vector();
//        try {
//            java.io.BufferedReader in = new java.io.BufferedReader(new java.io.FileReader((File) vecFiles.get(0)));
//            String str = in.readLine();
//            str = in.readLine();
//            while ((str = in.readLine()) != null) {
//                String[] data = str.split("\t");
//                vecProbes.add(data[0]);
//                nrProbes++;
//            }
//            in.close();
//        } catch (IOException e) {
//            e.printStackTrace();
//            System.out.println(e.getMessage());
//            System.exit(-1);
//        }
//        
//        System.out.println (nrProbes);
//        
//        ExpressionDataset_lude dataset1 = new ExpressionDataset_lude(nrProbes, nrSamples);
//        ExpressionDataset_lude dataset2 = new ExpressionDataset_lude(nrProbes, nrSamples);
//        ExpressionDataset_lude dataset3 = new ExpressionDataset_lude(nrProbes, nrSamples);
//        //ExpressionDataset dataset3 = new ExpressionDataset(nrProbes * 2, nrSamples);
//        for (int p=0; p<vecProbes.size(); p++) {
//            dataset1.probeNames[p] = (String) vecProbes.get(p);
//            dataset2.probeNames[p] = (String) vecProbes.get(p);
//            dataset3.probeNames[p] = (String) vecProbes.get(p);
//            //dataset3.probeNames[p * 2] = "M-" + (String) vecProbes.get(p);
//            //dataset3.probeNames[p * 2 + 1] = "U-" + (String) vecProbes.get(p);
//        }
//        for (int f=0; f<nrSamples; f++) {
//            File currentFile = (File) vecFiles.get(f);
//            System.out.println("Processing:\t" + f + "\t" + currentFile.getAbsolutePath());
//            try {
//                java.io.BufferedReader in = new java.io.BufferedReader(new java.io.FileReader(currentFile));
//                String str = in.readLine();
//                String[] data = str.split("\t");
//                dataset1.sampleNames[f] = data[1];
//                dataset2.sampleNames[f] = data[1];
//                dataset3.sampleNames[f] = data[1];
//                str = in.readLine();
//                data = str.split("\t");
//                int columnM = -1;
//                int columnU = -1;
//                for (int d=0; d<data.length; d++) {
//                    if (data[d].toLowerCase().trim().replace(" ", "_").equals("methylated_signal_intensity_(m)")) {
//                        columnM = d;
//                    }
//                    if (data[d].toLowerCase().trim().replace(" ", "_").equals("un-methylated_signal_intensity_(u)")) {
//                        columnU = d;
//                    }
//                }
//                System.out.println(columnM + "\t" + columnU);
//                int p = 0;
//                while ((str = in.readLine()) != null) {
//                    data = str.split("\t");
//                    if (data[columnM].equals("NA")) {
//                        data[columnM] = "0";
//                    }
//                    if (data[columnU].equals("NA")) {
//                        data[columnU] = "0";
//                    }
//                    dataset1.rawData[p][f] = Double.parseDouble(data[columnM]);
//                    dataset2.rawData[p][f] = Double.parseDouble(data[columnU]);
//                    
//                    dataset3.rawData[p][f] = dataset1.rawData[p][f] / (dataset1.rawData[p][f] + dataset2.rawData[p][f]);
//                    //dataset3.rawData[p * 2][f] = Double.parseDouble(data[columnM]);
//                    //dataset3.rawData[p * 2 + 1][f] = Double.parseDouble(data[columnU]);
//                    
//                    String probe = (String) vecProbes.get(p);
//                    if (!data[0].equals(probe)) {
//                        System.out.println("Error!:\t" + f + "\t" + data[0] + "\t" + probe);
//                    }
//                    p++;
//                }
//                in.close();
//            } catch (IOException e) {
//                e.printStackTrace();
//                System.out.println(e.getMessage());
//                System.exit(-1);
//            }
//            
//        }
//        dataset1.recalculateHashMaps();
//        dataset2.recalculateHashMaps();
//        dataset3.recalculateHashMaps();
//        
//        dataset1.save("/Users/lude/Downloads/TCGADataM.txt");
//        dataset2.save("/Users/lude/Downloads/TCGADataU.txt");
//        dataset3.save("/Users/lude/Downloads/TCGADataBeta.txt");
//        
//        System.exit(0);
//        
//    }
//    
//    public void compareTCGAEVsWithBatches() {
//        
//        HashMap hashSampleDisease = new HashMap();
//        String fileDir = "/Users/lude/Downloads/774bde41-d002-447f-805b-0da8d94b30fa/METADATA/JHU_USC__HumanMethylation27/";
//        File file = new File(fileDir);
//        File[] files = file.listFiles();
//        for (int f=0; f<files.length; f++) {
//            try {
//                if (files[f].getAbsolutePath().endsWith(".sdrf.txt")) {
//                    java.io.BufferedReader in = new java.io.BufferedReader(new java.io.FileReader(files[f]));
//                    String str = in.readLine();
//                    while ((str = in.readLine()) != null) {
//                        String[] data = str.split("\t");
//                        String disease = data[7].replace("jhu-usc.edu_", "").replace(".HumanMethylation27", "").replace(".adf.txt", "");
//                        disease = disease.split("\\.")[0];
//                        hashSampleDisease.put(data[0].trim(), disease);
//                    }
//                    in.close();
//                }
//            } catch (IOException e) {
//                e.printStackTrace();
//                System.out.println(e.getMessage());
//                System.exit(-1);
//            }        
//        }
//        
//        
//        ExpressionDataset_lude dataset = new ExpressionDataset_lude("/Users/lude/Downloads/PCATCGADataM+U/eigenvectors.txt");
//        
//        HashMap hashBatches = new HashMap();
//        HashMap hashDiseases = new HashMap();
//        for (int s=0; s<dataset.nrProbes; s++) {
//            String sample = dataset.probeNames[s];
//            String[] data = sample.split("-");
//            String batch  = data[5];
//            if (!hashBatches.containsKey(batch)) {
//                hashBatches.put(batch, hashBatches.size());
//            }
//            int batchID = ((Integer) hashBatches.get(batch)).intValue();
//            String disease = (String) hashSampleDisease.get(sample);
//            if (!hashDiseases.containsKey(disease)) {
//                hashDiseases.put(disease, hashDiseases.size());
//            }
//            int diseaseID = ((Integer) hashDiseases.get(disease)).intValue();
//            System.out.println(sample + "\t" + batch + "\t" + batchID + "\t" + disease + "\t" + diseaseID);
//        }
//        for (int s=0; s<dataset.nrProbes; s++) {
//            
//        }
//        System.out.println(hashBatches.size());
//        
//        System.exit(0);
//    }
//    
//    public void compareTCGAPCsWithFeatures() {
//        
//        HashMap hashProbesToInclude = new HashMap();
//        
//        String fileProbeAnnotation = "/Users/lude/Documents/DMG/Data/GPL8490/GPL8490_HumanMethylation27_270596_v.1.2.csv";
//        try {
//            java.io.BufferedReader in = new java.io.BufferedReader(new java.io.FileReader(new File(fileProbeAnnotation)));
//            String str = "";
//            while ((str = in.readLine()) != null) {
//                String[] data = str.split(",");
//                if (data.length>10) {
//                    hashProbesToInclude.put("M-" + data[0], null);
//                }
//            }
//            in.close();
//        } catch (IOException e) {
//            e.printStackTrace();
//            System.out.println(e.getMessage());
//            System.exit(-1);
//        }        
//        
//        //ExpressionDataset dataset = new ExpressionDataset("/Users/lude/Downloads/PCATCGADataUv2/principalcomponents.txt");
//        ExpressionDataset_lude dataset = new ExpressionDataset_lude("/Users/lude/Downloads/PCATCGADataM+U/principalcomponents.txt", "\t", hashProbesToInclude, null);
//
//        
//        double[] gcContent1 = new double[dataset.nrProbes];
//        double[] gcContent2 = new double[dataset.nrProbes];
//        double[] gcContent3 = new double[dataset.nrProbes];
//        double[] gcContent4 = new double[dataset.nrProbes];
//        try {
//            java.io.BufferedReader in = new java.io.BufferedReader(new java.io.FileReader(new File(fileProbeAnnotation)));
//            String str = "";
//            while ((str = in.readLine()) != null) {
//                String[] data = str.split(",");
//                if (data.length>10) {
//                    int probeIndex = ((Integer) dataset.hashProbes.get("M-" + data[0])).intValue();
//                    if (1==1) {
//                        byte[] bytes = data[4].getBytes();
//                        for (int b=0; b<bytes.length; b++) {
//                            if (bytes[b]==67 || bytes[b]==71) {
//                                gcContent1[probeIndex]++;
//                            }
//                        }
//                    }
//                    if (1==1) {
//                        byte[] bytes = data[6].getBytes();
//                        for (int b=0; b<bytes.length; b++) {
//                            if (bytes[b]==67 || bytes[b]==71) {
//                                gcContent2[probeIndex]++;
//                            }
//                        }
//                    }
//                    if (1==1) {
//                        byte[] bytes = data[15].getBytes();
//                        for (int b=0; b<bytes.length; b++) {
//                            if (bytes[b]==67 || bytes[b]==71) {
//                                gcContent3[probeIndex]++;
//                            }
//                        }
//                    }
//                    if (1==1) {
//                        byte[] bytes = data[16].getBytes();
//                        for (int b=0; b<bytes.length; b++) {
//                            if (bytes[b]==67 || bytes[b]==71) {
//                                gcContent4[probeIndex]++;
//                            }
//                        }
//                    }
//                    double base = 0;
//                    if (data[17].equals("A")) base = 0;
//                    if (data[17].equals("C")) base = 1;
//                    if (data[17].equals("G")) base = 2;
//                    if (data[17].equals("T")) base = 3;
//                    double grnRed = 0; 
//                    if (data[18].equals("Grn")) grnRed = 0;
//                    if (data[18].equals("Red")) grnRed = 1;
//                    System.out.println(probeIndex + "\t" + dataset.probeNames[probeIndex] + "\t" + gcContent1[probeIndex] + "\t" + gcContent2[probeIndex] + "\t" + gcContent3[probeIndex] + "\t" + gcContent4[probeIndex] + "\t" + base + "\t" + grnRed );
//                }
//            }
//            in.close();
//        } catch (IOException e) {
//            e.printStackTrace();
//            System.out.println(e.getMessage());
//            System.exit(-1);
//        }
//        
//        ExpressionDataset_lude datasetProbeProperties = new ExpressionDataset_lude("/Users/lude/Documents/DMG/Data/GPL8490/GPL8490ProbeSequencePropertiesPCA/principalcomponents.txt");
//        datasetProbeProperties.transposeDataset();
//        
//        
//        dataset.transposeDataset();;
//        for (int tc=0; tc<dataset.nrProbes; tc++) {
//            String output = tc + "\t" + JSci.maths.ArrayMath.correlation(gcContent1, dataset.rawData[tc]) + "\t" + JSci.maths.ArrayMath.correlation(gcContent2, dataset.rawData[tc]) + "\t" + JSci.maths.ArrayMath.correlation(gcContent3, dataset.rawData[tc]) + "\t" + JSci.maths.ArrayMath.correlation(gcContent4, dataset.rawData[tc]);
//            for (int n=0; n<datasetProbeProperties.nrProbes; n++) {
//                output+="\t" + JSci.maths.ArrayMath.correlation(datasetProbeProperties.rawData[n], dataset.rawData[tc]);
//            }
//            System.out.println(output);
//        }
//        
//        
//        System.exit(0);
//    }
//    
//    private double mean(double[] v) {
//        double sum = 0;
//        for (int k = 0; k < v.length; k++) sum += v[k];
//        return (sum / (double) v.length);
//    }
//
//    private double variance(double[] v, double mean) {
//        double ans = 0.0;
//        for (int i = 0; i < v.length; i++) ans += (v[i] - mean) * (v[i] - mean);
//        return ans / (v.length - 1);
//    }
//
//    public void parseDataPeterBram() {
//
//
//        ExpressionDataset_lude dataset = new ExpressionDataset_lude("/Users/lude/Downloads/ForLude061011/expressionphenotype_240811_lastcolumnremoved.txt");
//        dataset.transposeDataset();
//
//        ExpressionDataset_lude datasetSorted = new ExpressionDataset_lude(dataset.nrSamples, dataset.nrProbes);
//
//        //Perform QN, allowing for missing values:
//        for (int s=0; s<dataset.nrSamples; s++) {
//            double[] vals = new double[dataset.nrProbes];
//            for (int p=0; p<dataset.nrProbes; p++) {
//                vals[p] = dataset.rawData[p][s];
//            }
//            Arrays.sort(vals);
//            datasetSorted.rawData[s] = vals;
//        }
//        datasetSorted.transposeDataset();
//        double[] dist = new double[dataset.nrProbes];
//        for (int p=0; p<dataset.nrProbes; p++) {
//            dist[p] = JSci.maths.ArrayMath.median(datasetSorted.rawData[p]);
//            //System.out.println(p + "\t" + dist[p] + "\t" + JSci.maths.ArrayMath.mean(datasetSorted.rawData[p]));
//        }
//        Arrays.sort(dist);
//        double min = JSci.maths.ArrayMath.min(dist);
//        double max = JSci.maths.ArrayMath.max(dist);
//        int nrBins = 25;
//        int[] expDist = new int[nrBins];
//        for (int p=0; p<dataset.nrProbes; p++) {
//            int bin = (int) Math.round( (dist[p] - min) / (max - min) * (double) (nrBins - 1));
//            expDist[bin]++;
//        }
//        System.out.println(min + "\t" + max);
//        for (int bin=0; bin<nrBins; bin++) {
//            double exp = min + ((double) bin / (double) nrBins) * (max - min);
//            System.out.println(bin + "\t" + exp + "\t" + expDist[bin]);
//        }
//
//        dataset.transposeDataset();
//        dataset.save(dataset.fileName + ".Transposed.txt");
//        dataset.transposeDataset();
//
//        //Perform quantile normalization:
//        System.out.println(dataset.nrSamples);
//        for (int s=0; s<dataset.nrSamples; s++) {
//            Vector vec1 = new Vector();
//            for (int p=0; p<dataset.nrProbes; p++) {
//                vec1.add(dataset.rawData[p][s]);
//            }
//            double[] vals1 = new double[vec1.size()];
//            for (int v = 0; v < vals1.length; v++) {
//                vals1[v] = ((Double) vec1.get(v)).doubleValue();
//            }
//            umcg.genetica.util.Rank rank = new umcg.genetica.util.Rank(vals1, 0d);
//            double[] valsRanked = rank.getRanks();
//            for (int v=0; v<vals1.length; v++) valsRanked[v]--;
//            for (int v=0; v<vals1.length; v++) {
//                int distIndex = (int) valsRanked[v];
//                vals1[v] = dist[distIndex];
//            }
//            int itr = 0;
//            for (int p=0; p<dataset.nrProbes; p++) {
//                dataset.rawData[p][s] = vals1[itr];
//                itr++;
//            }
//            System.out.println(s + "\t" + dataset.sampleNames[s] + "\t" + JSci.maths.ArrayMath.median(vals1) + "\t" + JSci.maths.ArrayMath.mean(vals1) + "\t" + JSci.maths.ArrayMath.standardDeviation(vals1) + "\t" + JSci.maths.ArrayMath.min(vals1) + "\t" + JSci.maths.ArrayMath.max(vals1));
//        }
//
//        dataset.transposeDataset();
//        dataset.save(dataset.fileName + ".Transposed.QuantileNormalized.txt");
//        dataset.transposeDataset();
//
//        for (int p=0; p<dataset.nrProbes; p++) {
//            double mean = JSci.maths.ArrayMath.mean(dataset.rawData[p]);
//            if (mean > 100) {
//                System.out.println(p + "\t" + dataset.probeNames[p] + "\t" + mean);
//            }
//        }
//
//        dataset.standardNormalizeData();
//        dataset.transposeDataset();
//        dataset.standardNormalizeData();
//        //dataset.transposeDataset();
//        
//        File fileDatasetEVs = new File(dataset.fileName + ".PCAOverSamplesEigenvectorsAfterQNProbesCenteredScaled.txt");
//        if (!fileDatasetEVs.exists()) {
//            int sampleCountMinusOne = dataset.nrSamples - 1;
//            double[][] corr = new double[dataset.nrProbes][dataset.nrProbes];
//            for (int f = 0; f < dataset.nrProbes; f++) {
//                for (int g = f + 1; g < dataset.nrProbes; g++) {
//                    double covarianceInterim = 0;
//                    for (int s = 0; s < dataset.nrSamples; s++) {
//                        covarianceInterim += dataset.rawData[f][s] * dataset.rawData[g][s];
//                    }
//                    corr[f][g] = covarianceInterim / (double) (sampleCountMinusOne);
//                    corr[g][f] = corr[f][g];
//                }
//                corr[f][f] = 1;
//                System.out.println(f);
//            }
//
//            System.out.println("EVD decomposition ongoing:");
//            Jama.EigenvalueDecomposition eig = eigenValueDecomposition(corr);
//            double[] eigenValues = eig.getRealEigenvalues();
//
//            double[][] eigenVectors = new double[corr.length][corr.length];
//            ExpressionDataset_lude datasetEVs = new ExpressionDataset_lude(corr.length, corr.length);
//            for (int f = 0; f < dataset.nrProbes; f++) {
//                datasetEVs.sampleNames[f] = dataset.probeNames[f];
//            }
//            for (int pca=0; pca<corr.length; pca++) {
//                datasetEVs.probeNames[pca] = "Comp" + (pca + 1);
//                eigenVectors[pca] = getEigenVector(eig, pca);
//                datasetEVs.rawData[pca] = getEigenVector(eig, pca);
//                System.out.println("PCA" + (pca+1) + "\tExplainedVariance:\t" + getEigenValueVar(eigenValues, pca) + "\tEigenvalue:\t" + eigenValues[eigenValues.length - 1 - pca]);
//            }
//            datasetEVs.save(dataset.fileName + ".PCAOverSamplesEigenvectorsAfterQNProbesCenteredScaled.txt");
//        }
//
//
//
//
//        System.exit(0);
//    }
//
//    public void convertPCsToGeneLevel () {
//
//        ExpressionDataset_lude datasetPCs = new ExpressionDataset_lude("/Users/lude/Downloads/PCATCGABeta/PCAOverTCGA+GPL8490Combined/principalcomponents.txt");
//        
//        int nrPCs = 100;
//
//         HashMap hashAnnotation = new HashMap();
//         Vector vecUniqueGenes = new Vector();
//         HashMap hashUniqueGenes = new HashMap();
//         try {
//            java.io.BufferedReader in = new java.io.BufferedReader(new java.io.FileReader(new File("/Users/lude/Documents/DMG/Data/GPL8490/ProbeAnnotationNCBIBuild36.txt")));
//            String str = "";
//            while ((str = in.readLine()) != null) {
//                String[] data = str.split("\t");
//                if (datasetPCs.hashProbes.containsKey(data[0])) {
//                    String annotation = "-";
//                    if (data.length > 3) {
//                        annotation = data[3];
//                        if (!hashUniqueGenes.containsKey(annotation)) {
//                            vecUniqueGenes.add(annotation);
//                            hashUniqueGenes.put(annotation, data[0]);
//                        } else {
//                            hashUniqueGenes.put(annotation, (String) hashUniqueGenes.get(annotation) + "," + data[0]);
//                        }
//                    }
//                    hashAnnotation.put(data[0], annotation);
//                }
//            }
//        } catch (IOException e) {
//            e.printStackTrace();
//            System.out.println(e.getMessage());
//            System.exit(-1);
//        }
//        System.out.println("Number of unique genes:\t" + vecUniqueGenes.size());
//
//        
//        /*
//        ExpressionDataset dataset = new ExpressionDataset("/Users/lude/Documents/DMG/Data/GPL8490/GPL8490_family.soft.intensities.txt.binary.qced.binary");
//        dataset.transposeDataset();
//        dataset.standardNormalizeData();
//
//        
//        File fileDatasetEVs = new File(dataset.fileName + ".PCAOverSamplesEigenvectors.txt");
//        if (!fileDatasetEVs.exists()) {
//            int sampleCountMinusOne = dataset.nrSamples - 1;
//            double[][] corr = new double[dataset.nrProbes][dataset.nrProbes];
//            for (int f = 0; f < dataset.nrProbes; f++) {
//                for (int g = f + 1; g < dataset.nrProbes; g++) {
//                    double covarianceInterim = 0;
//                    for (int s = 0; s < dataset.nrSamples; s++) {
//                        covarianceInterim += dataset.rawData[f][s] * dataset.rawData[g][s];
//                    }
//                    corr[f][g] = covarianceInterim / (double) (sampleCountMinusOne);
//                    corr[g][f] = corr[f][g];
//                }
//                corr[f][f] = 1;
//                System.out.println(f);
//            }
//
//            System.out.println("EVD decomposition ongoing:");
//            Jama.EigenvalueDecomposition eig = eigenValueDecomposition(corr);
//            double[] eigenValues = eig.getRealEigenvalues();
//
//            double[][] eigenVectors = new double[corr.length][corr.length];
//            ExpressionDataset datasetEVs = new ExpressionDataset(corr.length, corr.length);
//            for (int f = 0; f < dataset.nrProbes; f++) {
//                datasetEVs.sampleNames[f] = dataset.probeNames[f];
//            }
//            for (int pca=0; pca<corr.length; pca++) {
//                datasetEVs.probeNames[pca] = "Comp" + (pca + 1);
//                eigenVectors[pca] = getEigenVector(eig, pca);
//                datasetEVs.rawData[pca] = getEigenVector(eig, pca);
//                System.out.println("PCA" + (pca+1) + "\tExplainedVariance:\t" + getEigenValueVar(eigenValues, pca) + "\tEigenvalue:\t" + eigenValues[eigenValues.length - 1 - pca]);
//            }
//            datasetEVs.save(dataset.fileName + ".PCAOverSamplesEigenvectors.txt");
//        }
//
//        ExpressionDataset datasetEVs = new ExpressionDataset(dataset.fileName + ".PCAOverSamplesEigenvectors.txt");
//
//
//        //Calculate PC Scores:
//        ExpressionDataset datasetPCs = new ExpressionDataset(dataset.nrSamples, nrPCs);
//        for (int p=0; p<dataset.nrSamples; p++) {
//            datasetPCs.probeNames[p] = dataset.sampleNames[p];
//        }
//        for (int pca=0; pca<nrPCs; pca++) {
//            datasetPCs.sampleNames[pca] = "Comp" + (pca + 1);
//            for (int p=0; p<dataset.nrSamples; p++) {
//                double val = 0;
//                for (int s=0; s<dataset.nrProbes; s++) {
//                    val+=datasetEVs.rawData[pca][s] * dataset.rawData[s][p];
//                }
//                datasetPCs.rawData[p][pca] = val;
//                //System.out.println(pca + "\t" + p + "\t" + dataset.sampleNames[p] + "\t" + val + "\t" + (String) hashAnnotation.get(dataset.sampleNames[p]));
//            }
//            if (pca == 0) {
//                double[] valsX = new double[dataset.nrSamples];
//                dataset.transposeDataset();
//                for (int p=0; p<dataset.nrProbes; p++) {
//                    valsX[p] = JSci.maths.ArrayMath.mean(dataset.rawData[p]);
//                }
//                dataset.transposeDataset();
//                double[] valsY = new double[dataset.nrSamples];
//                for (int p=0; p<dataset.nrSamples; p++) {
//                    valsY[p] = datasetPCs.rawData[p][0];
//                }
//                System.out.println("Correlation PC0 and average:\t" + JSci.maths.ArrayMath.correlation(valsX, valsY));
//            }
//            System.out.println(pca);
//        }
//        datasetPCs.recalculateHashMaps();
//        datasetPCs.save(dataset.fileName + ".PCAOverSamplesPrincipalComponents.txt");
//         * 
//         */
//        
//
//        ExpressionDataset_lude datasetPCsGeneLevel = new ExpressionDataset_lude(vecUniqueGenes.size(), nrPCs);
//        for (int pca=0; pca<nrPCs; pca++) {
//            datasetPCsGeneLevel.sampleNames[pca] = "Comp" + (pca + 1);
//        }
//        for (int v=0; v<datasetPCsGeneLevel.nrProbes; v++) {
//            String gene = (String) vecUniqueGenes.get(v);
//            datasetPCsGeneLevel.probeNames[v] = gene;
//            String[] probes = ((String) hashUniqueGenes.get(gene)).split(",");
//            for (int p=0; p<probes.length; p++) {
//                String probe = probes[p];
//                int probeID = ((Integer) datasetPCs.hashProbes.get(probe)).intValue();
//                for (int pca=0; pca<nrPCs; pca++) {
//                    datasetPCsGeneLevel.rawData[v][pca]+=datasetPCs.rawData[probeID][pca];
//                }
//            }
//            for (int pca=0; pca<nrPCs; pca++) {
//                datasetPCsGeneLevel.rawData[v][pca]/=(double) probes.length;
//            }
//        }
//        datasetPCsGeneLevel.save(datasetPCs.fileName + ".GeneLevel.txt");
//
//        System.exit(0);
//    }
//
//    public void performPCA() {
//        ExpressionDataset_lude dataset = new ExpressionDataset_lude("/Users/lude/Documents/DMG/Data/GPL8490/GPL8490_family.soft.intensities.txt.binary.qced.binary");
//        dataset.save("/Users/lude/Documents/DMG/Data/GPL8490/GPL8490_family.soft.intensities.txt.binary.qced.txt");
//
//        for (int p=0; p<dataset.nrProbes; p++) {
//
//            umcg.genetica.util.Rank rank = new umcg.genetica.util.Rank(dataset.rawData[p], 0d);
//            dataset.rawData[p] = rank.getRanks();
//
//        }
//
//         HashMap hashAnnotation = new HashMap();
//         try {
//            java.io.BufferedReader in = new java.io.BufferedReader(new java.io.FileReader(new File("/Users/lude/Documents/DMG/Data/GPL8490/ProbeAnnotationNCBIBuild36.txt")));
//            String str = "";
//            while ((str = in.readLine()) != null) {
//                String[] data = str.split("\t");
//                String annotation = "-";
//                if (data.length > 3) {
//                    annotation = data[3];
//                }
//                hashAnnotation.put(data[0], annotation);
//            }
//        } catch (IOException e) {
//            e.printStackTrace();
//            System.out.println(e.getMessage());
//            System.exit(-1);
//        }
//
//
//        String outputFile = "/Users/lude/Documents/DMG/Data/GPL8490/Plots/";
//        dataset.standardNormalizeData();
//        int sampleCountMinusOne = dataset.nrSamples - 1;
//        for (int f = 0; f < dataset.nrProbes; f++) {
//
//            String annotation = (String) hashAnnotation.get(dataset.probeNames[f]);
//            //if (annotation.contains("F13A1")) {
//                for (int g = f + 1; g < dataset.nrProbes; g++) {
//                //for (int g = 0; g < dataset.nrProbes; g++) {
//                    //Calculate correlation:
//                    double covarianceInterim = 0;
//                    for (int s = 0; s < dataset.nrSamples; s++) {
//                        covarianceInterim += dataset.rawData[f][s] * dataset.rawData[g][s];
//                    }
//                    double covariance = covarianceInterim / (double) (sampleCountMinusOne);
//                    double correlation = covariance;
//                    if (correlation > 0.90) {
//                        System.out.println(dataset.probeNames[f] + "\t" + (String) hashAnnotation.get(dataset.probeNames[f]) + "\t" + dataset.probeNames[g] + "\t" + (String) hashAnnotation.get(dataset.probeNames[g]) + "\t" + correlation);
//                    }
//                    if (correlation > 0.90) {
//                        //System.out.println(dataset.probeNames[f] + "\t" + (String) hashAnnotation.get(dataset.probeNames[f]) + "\t" + dataset.probeNames[g] + "\t" + (String) hashAnnotation.get(dataset.probeNames[g]) + "\t" + correlation);
//
//                        int width = 500 + 200;
//                        int height = 500 + 200;
//                        int marginLeft = 100; int marginRight = 100; int marginTop = 100; int marginBottom = 100;
//                        double innerWidth = width - marginLeft - marginRight;
//                        double innerHeight = height - marginBottom - marginTop;
//                        double x0 = marginLeft;
//                        double x1 = x0 + innerWidth;
//                        double y0 = marginTop;
//                        double y1 = y0 + innerHeight;
//                        double centerX = (x1 + x0) / 2;
//                        double centerY = (y1 + y0) / 2;
//
//                        BufferedImage bimage = new BufferedImage(width, height, BufferedImage.TYPE_INT_RGB);
//                        Graphics2D g2d = bimage.createGraphics();
//                        g2d.setRenderingHint(RenderingHints.KEY_ANTIALIASING, RenderingHints.VALUE_ANTIALIAS_ON);
//                        g2d.setColor(new Color(255, 255, 255));
//                        g2d.fillRect(0,0, width, height);
//
//                        g2d.setColor(new Color(0, 0, 0));
//                        g2d.setComposite(java.awt.AlphaComposite.getInstance(java.awt.AlphaComposite.SRC_OVER, 1.0f));
//                        double minF = JSci.maths.ArrayMath.min(dataset.rawData[f]);
//                        double maxF = JSci.maths.ArrayMath.max(dataset.rawData[f]);
//                        double minG = JSci.maths.ArrayMath.min(dataset.rawData[g]);
//                        double maxG = JSci.maths.ArrayMath.max(dataset.rawData[g]);
//                        for (int s=0; s<dataset.nrSamples; s++) {
//                            int posX = (int) Math.round(x0 + innerWidth * (dataset.rawData[f][s] - minF) / (maxF - minF));
//                            int posY = (int) Math.round(y1 - innerHeight * (dataset.rawData[g][s] - minG) / (maxG - minG));
//
//                            float hue = (float) s / (float) dataset.nrSamples;
//                            g2d.setColor(g2d.getColor().getHSBColor(hue, 1.0f, 1.0f));
//
//                            g2d.fillOval(posX - 3, posY - 3, 7, 7);
//                        }
//                        String correlationString = (new java.text.DecimalFormat("0.000", new java.text.DecimalFormatSymbols(java.util.Locale.US))).format(correlation);
//                        try {
//                            javax.imageio.ImageIO.write(bimage, "png", new File(outputFile + "/" + correlationString + "-Coexpression-" + dataset.probeNames[f] + "-" + dataset.probeNames[g] + ".png"));
//                        } catch (IOException e) {
//                            System.out.println(e.getMessage());
//                            e.printStackTrace();
//                        }
//
//                    }
//                }
//            //}
//            if (f % 100 == 99) {
//                System.out.println((f + 1) + " Probes processed");
//            }
//        }
//
//        System.exit(0);
//    }
//
//    public void performCoexpressionOnPCScoresFromPCAOverSamples() {
//
//        int nrTCs = 100;
//        HashMap hashTCsToIncludeHuman = new HashMap();
//        for (int tc=0; tc<nrTCs; tc++) {
//            hashTCsToIncludeHuman.put("Comp" + String.valueOf(tc + 1), null);
//        }
//
//        ExpressionDataset_lude datasetMapping = new ExpressionDataset_lude("/Users/lude/Documents/DMG/Data/GPL8490/ProbeAnnotationNCBIBuild36Sorted.txt");
//        
//        HashMap hashAutosomalProbes = new HashMap();
//        for (int p=0; p<datasetMapping.nrProbes; p++) {
//            //if (datasetMapping.rawData[p][0]<23) {
//                hashAutosomalProbes.put(datasetMapping.probeNames[p], null);
//            //}
//        }
//        datasetMapping = new ExpressionDataset_lude("/Users/lude/Documents/DMG/Data/GPL8490/ProbeAnnotationNCBIBuild36Sorted.txt", "\t", hashAutosomalProbes, null);
//        
//        //ExpressionDataset datasetInitial = new ExpressionDataset("/Users/lude/Documents/DMG/Data/GPL8490/PCAOver4432SamplesProbesCenteredScaled/principalcomponents.txt", "\t", datasetMapping.hashProbes, hashTCsToIncludeHuman);
//        //ExpressionDataset datasetInitial = new ExpressionDataset("/Users/lude/Downloads/PCATCGADataUv2/principalcomponents.txt", "\t", datasetMapping.hashProbes, hashTCsToIncludeHuman);
//        ExpressionDataset_lude datasetInitial = new ExpressionDataset_lude("/Users/lude/Downloads/PCATCGABeta/PCAOverTCGA+GPL8490Combined/principalcomponents.txt", "\t", datasetMapping.hashProbes, hashTCsToIncludeHuman);
//         
//        ExpressionDataset_lude dataset = datasetInitial;
//        
//        datasetMapping = new ExpressionDataset_lude("/Users/lude/Documents/DMG/Data/GPL8490/ProbeAnnotationNCBIBuild36Sorted.txt", "\t", datasetInitial.hashProbes, null);
//        
//        if (1==2) {
//            dataset = new ExpressionDataset_lude(datasetInitial.nrProbes, datasetInitial.nrSamples);
//            for (int tc=0; tc<nrTCs; tc++) {
//                dataset.sampleNames[tc] = datasetInitial.sampleNames[tc];
//            }
//            for (int p=0; p<dataset.nrProbes; p++) {
//                dataset.probeNames[p] = datasetMapping.probeNames[p];
//                int probeIndex = ((Integer) datasetInitial.hashProbes.get(datasetMapping.probeNames[p])).intValue();
//                for (int tc=0; tc<nrTCs; tc++) {
//                    dataset.rawData[p][tc] = datasetInitial.rawData[probeIndex][tc];
//                }
//            }
//            dataset.recalculateHashMaps();
//            dataset.save(datasetInitial.fileName + ".sorted.txt");
//
//            for (int tc=0; tc<nrTCs; tc++) {
//                int lag = 2;
//                double[] vals1 = new double[dataset.nrProbes - lag];
//                double[] vals2 = new double[dataset.nrProbes - lag];
//                for (int p=0; p<dataset.nrProbes - lag; p++) {
//                    vals1[p] =dataset.rawData[p][tc];
//                    vals2[p] =dataset.rawData[p + lag][tc];                
//                    if (tc==5) {
//                        //System.out.println(tc + "\t" + p + "\t" + vals1[p]);
//                    }
//                }
//                System.out.println(tc + "\t" + JSci.maths.ArrayMath.correlation(vals1, vals2));
//            }
//
//            if (1==1) System.exit(0);
//        }
//
//         HashMap hashAnnotation = new HashMap();
//         try {
//            java.io.BufferedReader in = new java.io.BufferedReader(new java.io.FileReader(new File("/Users/lude/Documents/DMG/Data/GPL8490/ProbeAnnotationNCBIBuild36.txt")));
//            String str = "";
//            while ((str = in.readLine()) != null) {
//                String[] data = str.split("\t");
//                String annotation = "-";
//                if (data.length > 3) {
//                    annotation = data[3];
//                }
//                hashAnnotation.put(data[0], annotation);
//            }
//        } catch (IOException e) {
//            e.printStackTrace();
//            System.out.println(e.getMessage());
//            System.exit(-1);
//        }
//
//
//        String outputFile = "/Users/lude/Documents/DMG/Data/GPL8490/PCAOver4432SamplesProbesCenteredScaled/Plots/";
//        dataset.transposeDataset();
//        dataset.standardNormalizeData();
//        dataset.transposeDataset();
//        dataset.standardNormalizeData();
//        int sampleCountMinusOne = dataset.nrSamples - 1;
//
//        int[] count = new int[dataset.nrProbes];
//        
//        
//        if (1==2) {
//            int[] degree = new int[dataset.nrProbes];
//
//            for (int f = 0; f < dataset.nrProbes; f++) {
//                for (int g = f + 1; g < dataset.nrProbes; g++) {
//                    if (f!=g) {
//                        //Calculate correlation:
//                        double covarianceInterim = 0;
//                        for (int s = 0; s < dataset.nrSamples; s++) {
//                            covarianceInterim += dataset.rawData[f][s] * dataset.rawData[g][s];
//                        }
//                        double covariance = covarianceInterim / (double) (sampleCountMinusOne);
//                        double correlation = covariance;
//                        if (correlation > 0.4) {
//                            degree[f]++;
//                            degree[g]++;
//                        }
//                    }
//                }
//                if (f%100==0) System.out.println(f);
//            }
//
//            int[] degreeDist = new int[dataset.nrProbes];
//            for (int f = 0; f < dataset.nrProbes; f++) {
//                degreeDist[degree[f]]++;
//            }
//            System.out.println("Degree distribution:");
//            for (int f = 0; f < dataset.nrProbes; f++) {
//                if (degreeDist[f]>0) {
//                    System.out.println(f + "\t" + degreeDist[f]);
//                }
//            }
//            System.out.println("");
//        }        
//        
//        //String[] alsGenes = {"SOD1","ALS2","SETX","FUS","VAPB","ANG","TARDBP","OPTN","VCP","APEX1","ATXN2","CHMP2B","NEFH","SMN1","SMN2","PRPH","VEGFA","UNC13A","C9orf72","NIPA1","PFN1"};
//        //String[] alsGenes = {"RGS1","REL","AHSA2","IL18RAP","IL18R1","IL1RL1","IL1RL2","ITGA4","UBE2E3","CTLA4","ICOS","CD28","CCR1","CCR2","CCRL2","CCR3","CCR5","CCR9","IL12A","LPP","IL2","IL21","HLA-DQA1","HLA-DQB1","TNFAIP3","TAGAP","SH2B3","PTPN2","TNFRSF14","MMEL1","RUNX3","PLEK","CCR4","CD80","KTELC1","BACH2","MAP3K7","PTPRK","THEMIS","ZMIZ1","ETS1","CIITA","SOCS1","CLEC16A","ICOSLG","PARK7","TNFRSF9","NFIA","CD247","FASLG","TNFSF18","TNFSF4","FRMD4B","IRF4","ELMO1","ZFP36L1","UBE2L3","YDJC","TLR7","TLR8"};
//        //String[] alsGenes = {"RGS1","REL","AHSA2","IL18RAP","IL18R1","IL1RL1","IL1RL2","ITGA4","UBE2E3","CTLA4","ICOS","CD28","CCR1","CCR2","CCRL2","CCR3","CCR5","CCR9","IL12A","LPP","IL2","IL21","HLA-DQA1","HLA-DQB1","TNFAIP3","TAGAP","SH2B3","PTPN2","TNFRSF14","MMEL1","RUNX3","PLEK","CCR4","CD80","KTELC1","BACH2","MAP3K7","PTPRK","THEMIS","ZMIZ1","ETS1","CIITA","SOCS1","CLEC16A","ICOSLG","PARK7","TNFRSF9","NFIA","CD247","FASLG","TNFSF18","TNFSF4","FRMD4B","IRF4","ELMO1","ZFP36L1","UBE2L3","YDJC","TLR7","TLR8"};
//        //String[] alsGenes = {"FTO", "TMEM18", "MC4R", "GNPDA2", "BDNF", "NEGR1", "SH2B1", "ETV5", "MTCH2", "KCTD15", "identified", "SEC16B", "TFAP2B", "FAIM2", "NRXN3", "identified", "RBJ", "GPRC5B", "MAP2K5", "QPCTL", "TNNI3K", "SLC39A8", "FLJ35779", "LRRN6C", "TMEM160", "FANCL", "CADM2", "PRKD1", "LRP1B", "PTBP2", "MTIF3", "ZNF608", "RPL27A", "NUDT3", "APOB48R", "SULT1A2", "AC138894.2", "ATXN2L", "TUFM", "NDUFS3", "CUGBP1", "SEC16B", "TFAP2B", "FAIM2", "NRXN3", "RBJ", "GPRC5B", "MAP2K5", "QPCTL", "TNNI3K", "SLC39A8", "FLJ35779", "LRRN6C", "TMEM160", "FANCL", "CADM2", "PRKD1", "LRP1B", "PTBP2", "MTIF3", "ZNF608", "RPL27A", "NUDT3", "ADCY3", "POMC", "IQCK", "LBXCOR1", "GIPR", "HMGCR", "ZC3H4", "GTF3A", "TUB", "HMGA1"};
//        String[] alsGenes = {"LACTB", "TP53"};
//        //String[] alsGenes = {"OAS1", "TAGAP"};
//        
//        
//        for (int f = 0; f < dataset.nrProbes; f++) {
//
//            String annotation = (String) hashAnnotation.get(dataset.probeNames[f]);
//            boolean include = false;
//            //if (annotation.toUpperCase().startsWith("TARDBP".toUpperCase())) {
//            //    include = true;
//            //}
//            for (int a=0; a<alsGenes.length; a++) {
//                if (annotation.toUpperCase().equals(alsGenes[a].toUpperCase())) {
//                    include = true;
//                }
//            }
//            if (include) {
//                System.out.println(f + "\t" + annotation);
//                //for (int g = f + 1; g < dataset.nrProbes; g++) {
//                for (int g = 0; g < dataset.nrProbes; g++) {
//                    if (f!=g) {
//                        //Calculate correlation:
//                        double covarianceInterim = 0;
//                        for (int s = 0; s < dataset.nrSamples; s++) {
//                            covarianceInterim += dataset.rawData[f][s] * dataset.rawData[g][s];
//                        }
//                        double covariance = covarianceInterim / (double) (sampleCountMinusOne);
//                        double correlation = covariance;
//                        if (Math.abs(correlation) > 0.3) {
//                            System.out.println(dataset.probeNames[f] + "\t" + (String) hashAnnotation.get(dataset.probeNames[f]) + "\t" + dataset.probeNames[g] + "\t" + (String) hashAnnotation.get(dataset.probeNames[g]) + "\t" + correlation);
//                            count[g]++;
//                        }
//                        if (Math.abs(correlation) > 10.4) {
//                            //System.out.println(dataset.probeNames[f] + "\t" + (String) hashAnnotation.get(dataset.probeNames[f]) + "\t" + dataset.probeNames[g] + "\t" + (String) hashAnnotation.get(dataset.probeNames[g]) + "\t" + correlation);
//
//                            int width = 500 + 200;
//                            int height = 500 + 200;
//                            int marginLeft = 100; int marginRight = 100; int marginTop = 100; int marginBottom = 100;
//                            double innerWidth = width - marginLeft - marginRight;
//                            double innerHeight = height - marginBottom - marginTop;
//                            double x0 = marginLeft;
//                            double x1 = x0 + innerWidth;
//                            double y0 = marginTop;
//                            double y1 = y0 + innerHeight;
//                            double centerX = (x1 + x0) / 2;
//                            double centerY = (y1 + y0) / 2;
//
//                            BufferedImage bimage = new BufferedImage(width, height, BufferedImage.TYPE_INT_RGB);
//                            Graphics2D g2d = bimage.createGraphics();
//                            g2d.setRenderingHint(RenderingHints.KEY_ANTIALIASING, RenderingHints.VALUE_ANTIALIAS_ON);
//                            g2d.setColor(new Color(255, 255, 255));
//                            g2d.fillRect(0,0, width, height);
//
//                            g2d.setColor(new Color(0, 0, 0));
//                            g2d.setComposite(java.awt.AlphaComposite.getInstance(java.awt.AlphaComposite.SRC_OVER, 1.0f));
//                            double minF = JSci.maths.ArrayMath.min(dataset.rawData[f]);
//                            double maxF = JSci.maths.ArrayMath.max(dataset.rawData[f]);
//                            double minG = JSci.maths.ArrayMath.min(dataset.rawData[g]);
//                            double maxG = JSci.maths.ArrayMath.max(dataset.rawData[g]);
//                            for (int s=0; s<dataset.nrSamples; s++) {
//                                int posX = (int) Math.round(x0 + innerWidth * (dataset.rawData[f][s] - minF) / (maxF - minF));
//                                int posY = (int) Math.round(y1 - innerHeight * (dataset.rawData[g][s] - minG) / (maxG - minG));
//
//                                float hue = (float) s / (float) dataset.nrSamples;
//                                g2d.setColor(g2d.getColor().getHSBColor(hue, 1.0f, 1.0f));
//
//                                g2d.fillOval(posX - 3, posY - 3, 7, 7);
//                            }
//                            if (1==2)  {
//                                try {
//                                    javax.imageio.ImageIO.write(bimage, "png", new File(outputFile + "/" + correlation + "-Coexpression-" + dataset.probeNames[f] + "-" + dataset.probeNames[g] + ".png"));
//                                } catch (IOException e) {
//                                    System.out.println(e.getMessage());
//                                    e.printStackTrace();
//                                }
//                            }
//
//                        }
//                    }
//                }
//            }
//            if (f % 100 == 99) {
//                //System.out.println((f + 1) + " Probes processed");
//            }
//        }
//
//        for (int f = 0; f < dataset.nrProbes; f++) {
//            if (count[f]>1) {
//                System.out.println(dataset.probeNames[f] + "\t" + (String) hashAnnotation.get(dataset.probeNames[f]) + "\t" + count[f]);
//            }
//        }
//
//        System.exit(0);
//    }
//
//    public void comparePCScoresWithGPL570TCs() {
//
//        String sortedEigenvectorFile         = "/Users/lude/Documents/DMG/Data/HumanUCSC/PCAOverAllProbesGPL570/eigenvectors.txt.binary";
//        String probeAnnotationFile = "/Users/lude/Documents/DMG/Data/AffymetrixAnnotation2010-08-23/HG-U133_Plus_2/HG-U133_Plus_2.na31.annot.csv";
//
//        HashMap hashProbeToChr = new HashMap();
//        HashMap hashProbeToHGNC = new HashMap();
//        HashMap hashHGNCToProbe = new HashMap();
//        try {
//            java.io.BufferedReader in = new java.io.BufferedReader(new java.io.FileReader(new File(probeAnnotationFile)));
//            String str = "";
//            while ((str = in.readLine()) != null) {
//                if (str.startsWith("\"Probe Set ID")) break;
//            }
//            while ((str = in.readLine()) != null) {
//                String[] data = str.trim().split("\",\"");
//                String mapping = data[12].replace("//", "\t");
//                String[] mappings = mapping.split("\t");
//                if (mappings.length==3) {
//                    String[] mappingsSplit = mappings[0].split(":");
//                    String chr = mappingsSplit[0];
//                    int chrInt = -1;
//                    try {
//                        chrInt = Integer.parseInt(chr.replace("chr", ""));
//                    } catch (Exception e) {
//                        if (chr.equals("chrX")) chrInt = 23;
//                        if (chr.equals("chrY")) chrInt = 24;
//                    }
//                    if (chrInt!=-1) {
//                        hashProbeToChr.put(data[0].substring(1), chrInt);
//                    }
//                }
//                String hgnc = data[14].trim(); if (hgnc.equals("---")) hgnc = "";
//                if (hgnc.length()>0) {
//                    hashProbeToHGNC.put(data[0].substring(1), hgnc.trim().replace(" ", ""));
//                    if (!hashHGNCToProbe.containsKey(hgnc)) {
//                        hashHGNCToProbe.put(hgnc, data[0].substring(1));
//                    } else {
//                        hashHGNCToProbe.put(hgnc, (String) hashHGNCToProbe.get(hgnc) + "," + data[0].substring(1));
//                    }
//                }
//            }
//        } catch (Exception e) {
//            System.out.println("Error:\t" + e.getMessage());
//            e.printStackTrace();
//        }
//
//        int nrTCs = 100;
//        HashMap hashTCsToIncludeHuman = new HashMap();
//        for (int tc=0; tc<nrTCs; tc++) {
//            hashTCsToIncludeHuman.put("Comp" + String.valueOf(tc + 1), null);
//        }
//
//        //ExpressionDataset dataset = new ExpressionDataset("/Users/lude/Documents/DMG/Data/GPL8490/PCAOVerSamplesGPL8490/principalcomponents.txt.binary", "\t", null, hashTCsToIncludeHuman);
//        //ExpressionDataset dataset = new ExpressionDataset("/Users/lude/Documents/DMG/Data/GPL8490/PCAOverSamplesProbesCenteredScaled/principalcomponents.txt", "\t", null, hashTCsToIncludeHuman);
//        //ExpressionDataset dataset = new ExpressionDataset("/Users/lude/Documents/DMG/Data/GPL8490/PCAOver4432SamplesProbesForcedNormalDistribution/principalcomponents.txt", "\t", null, hashTCsToIncludeHuman);
//        //ExpressionDataset dataset = new ExpressionDataset("/Users/lude/Downloads/PCATCGABeta/principalcomponents.txt", "\t", null, hashTCsToIncludeHuman);
//        ExpressionDataset_lude dataset = new ExpressionDataset_lude("/Users/lude/Downloads/PCATCGABeta/PCAOverTCGA+GPL8490Combined/principalcomponents.txt", "\t", null, hashTCsToIncludeHuman);
//
//        if (1==2) {
//            ExpressionDataset_lude datasetPCs = new ExpressionDataset_lude("/Users/lude/Documents/DMG/Data/GPL8490/PCAOVerSamplesGPL8490/principalcomponents.txt.binary", "\t", null, hashTCsToIncludeHuman);
//            datasetPCs.transposeDataset();
//            dataset.transposeDataset();
//
//            for (int s=0; s<dataset.nrProbes; s++) {
//                double r2Sum = 0;
//                for (int t=0; t<datasetPCs.nrProbes; t++) {
//                    double corr = JSci.maths.ArrayMath.correlation(datasetPCs.rawData[t], dataset.rawData[s]);
//                    double r2 = corr * corr;
//                    r2Sum+=r2;
//                }
//                System.out.println(s + "\t" + r2Sum);
//            }
//
//            System.exit(0);
//        }
//
//        HashMap hashAnnotation = new HashMap();
//         try {
//            java.io.BufferedReader in = new java.io.BufferedReader(new java.io.FileReader(new File("/Users/lude/Documents/DMG/Data/GPL8490/ProbeAnnotationNCBIBuild36.txt")));
//            String str = "";
//            while ((str = in.readLine()) != null) {
//                String[] data = str.split("\t");
//                String annotation = "-";
//                if (data.length > 3) {
//                    annotation = data[3];
//                }
//                hashAnnotation.put(data[0], annotation);
//            }
//        } catch (IOException e) {
//            e.printStackTrace();
//            System.out.println(e.getMessage());
//            System.exit(-1);
//        }
//
//
//        HashMap hash27KToAffy = new HashMap();
//        HashMap hashAffyTo27K = new HashMap();
//        HashMap hash27KChr = new HashMap();
//        for (int p=0; p<dataset.nrProbes; p++) {
//            String hgnc = (String) hashAnnotation.get(dataset.probeNames[p]);
//            if (hashHGNCToProbe.containsKey(hgnc)) {
//                String[] probes = ((String) hashHGNCToProbe.get(hgnc)).split(",");
//                if (!hashAffyTo27K.containsKey(probes[0])) {
//                    if (hashProbeToChr.containsKey(probes[0])) {
//                        int chr = ((Integer) hashProbeToChr.get(probes[0])).intValue();
//                        if (chr < 23) {
//                            hash27KToAffy.put(dataset.probeNames[p], probes[0]);
//                            hash27KChr.put(dataset.probeNames[p], chr);
//                            hashAffyTo27K.put(probes[0], dataset.probeNames[p]);
//                        }
//                    }
//                }
//            }
//        }
//        System.out.println(hash27KToAffy.size() + "\t" + hashAffyTo27K.size());
//
//        //dataset = new ExpressionDataset("/Users/lude/Documents/DMG/Data/GPL8490/PCAOverSamplesProbesCenteredScaled/principalcomponents.txt", "\t", hash27KToAffy, hashTCsToIncludeHuman);
//        //dataset = new ExpressionDataset("/Users/lude/Documents/DMG/Data/GPL8490/PCAOver4432SamplesProbesForcedNormalDistribution/principalcomponents.txt", "\t", hash27KToAffy, hashTCsToIncludeHuman);
//        //dataset = new ExpressionDataset("/Users/lude/Downloads/PCATCGABeta/principalcomponents.txt", "\t", hash27KToAffy, hashTCsToIncludeHuman);
//        dataset = new ExpressionDataset_lude("/Users/lude/Downloads/PCATCGABeta/PCAOverTCGA+GPL8490Combined/principalcomponents.txt", "\t", hash27KToAffy, hashTCsToIncludeHuman);
//        int[] probeChr = new int[dataset.nrProbes];
//        for (int p=0; p<dataset.nrProbes; p++) {
//            probeChr[p] = ((Integer) hash27KChr.get(dataset.probeNames[p])).intValue();
//            //System.out.println(p + "\t" + dataset.probeNames[p] + "\t" + probeChr[p]);
//        }
//
//        int nrTCsGPL570 = 250; //779;
//        HashMap hashGPL570TCsToInclude = new HashMap();
//        for (int tc=0; tc<nrTCsGPL570; tc++) {
//            hashGPL570TCsToInclude.put("Comp" + String.valueOf(tc + 1), null);
//        }
//
//        ExpressionDataset_lude datasetEVsInitial = new ExpressionDataset_lude(sortedEigenvectorFile, "\t", hashAffyTo27K, hashGPL570TCsToInclude);
//        ExpressionDataset_lude datasetEVs = new ExpressionDataset_lude(sortedEigenvectorFile, "\t", hashAffyTo27K, hashGPL570TCsToInclude);
//        for (int p=0; p<datasetEVs.nrProbes; p++) {
//            String probe27K = (String) hashAffyTo27K.get(datasetEVsInitial.probeNames[p]);
//            int probeID = ((Integer) dataset.hashProbes.get(probe27K)).intValue();
//            datasetEVs.probeNames[probeID] = datasetEVsInitial.probeNames[p];
//            for (int s=0; s<datasetEVs.nrSamples; s++) {
//                datasetEVs.rawData[probeID][s] = datasetEVsInitial.rawData[p][s];
//            }
//        }
//        datasetEVs.recalculateHashMaps();
//
//        dataset.transposeDataset();
//        dataset.standardNormalizeData();
//        datasetEVs.transposeDataset();
//
//        /*
//        for (int p=0; p<dataset.nrProbes; p++) {
//
//            jsc.util.Rank rank = new jsc.util.Rank(dataset.rawData[p], 0d);
//            dataset.rawData[p] = rank.getRanks();
//
//        }
//        for (int p=0; p<datasetEVs.nrProbes; p++) {
//
//            jsc.util.Rank rank = new jsc.util.Rank(datasetEVs.rawData[p], 0d);
//            datasetEVs.rawData[p] = rank.getRanks();
//
//        }
//         *
//         */
//        //dataset.standardNormalizeData();
//        //datasetEVs.standardNormalizeData();
//
//        System.out.println("Expression and methylation PCs that correlate with each other:");
//        for (int tc2=0; tc2<nrTCsGPL570; tc2++) {
//            double r2Sum = 0;
//            for (int tc=0; tc<nrTCs; tc++) {
//                double correlation = JSci.maths.ArrayMath.correlation(dataset.rawData[tc], datasetEVs.rawData[tc2]);
//                double r2 = correlation * correlation;
//                if (r2 > 0.01) {
//                    //System.out.println(tc + "\t" + tc2 + "\t" + correlation + "\t" + r2);
//                }
//                r2Sum+=r2;
//            }
//            System.out.println(tc2 + "\t" + r2Sum);
//        }
//
//        System.out.println("");
//
//        dataset.transposeDataset();
//        datasetEVs.transposeDataset();
//
//        dataset.standardNormalizeData();
//        datasetEVs.standardNormalizeData();
//
//
//        int width = 400 + 200 + 500;
//        int height = 400 + 200;
//        int marginLeft = 100; int marginRight = 600; int marginTop = 100; int marginBottom = 100;
//        double innerWidth = width - marginLeft - marginRight;
//        double innerHeight = height - marginBottom - marginTop;
//        double x0 = marginLeft;
//        double x1 = x0 + innerWidth;
//        double y0 = marginTop;
//        double y1 = y0 + innerHeight;
//        double centerX = (x1 + x0) / 2;
//        double centerY = (y1 + y0) / 2;
//
//        BufferedImage bimage = new BufferedImage(width, height, BufferedImage.TYPE_INT_RGB);
//        Graphics2D g2d = bimage.createGraphics();
//        g2d.setRenderingHint(RenderingHints.KEY_ANTIALIASING, RenderingHints.VALUE_ANTIALIAS_ON);
//        g2d.setColor(new Color(255, 255, 255));
//        g2d.fillRect(0,0, width, height);
//
//        g2d.setColor(new Color(0, 0, 0));
//        g2d.setComposite(java.awt.AlphaComposite.getInstance(java.awt.AlphaComposite.SRC_OVER, 1.0f));
//
//        long[][] corrDist = new long[201][201];
//        int nrSamples = dataset.nrSamples;
//        int nrSamplesGPL570 = datasetEVs.nrSamples;
//        double sampleCountMinusOne = dataset.nrSamples - 1;
//        double sampleCountMinusOneGPL570 = datasetEVs.nrSamples - 1;
//
//        for (int p=0; p<dataset.nrProbes; p++) {
//            for (int q=p+1; q<dataset.nrProbes; q++) {
//                if (probeChr[p]!=probeChr[q]) {
//                    double covarianceInterim = 0;
//                    for (int s = 0; s < nrSamples; s++) {
//                        covarianceInterim += dataset.rawData[p][s] * dataset.rawData[q][s];
//                    }
//                    double correlation = covarianceInterim / sampleCountMinusOne;
//
//                    double covarianceInterimGPL570 = 0;
//                    for (int s = 0; s < nrSamplesGPL570; s++) {
//                        covarianceInterimGPL570 += datasetEVs.rawData[p][s] * datasetEVs.rawData[q][s];
//                    }
//                    double correlationGPL570 = covarianceInterimGPL570 / sampleCountMinusOneGPL570;
//
//                    int correlationInt = (int) Math.round((correlation + 1d) * 100d);
//                    int correlationGPL570Int = (int) Math.round((correlationGPL570 + 1d) * 100d);
//                    //if (correlationGPL570Int<0) correlationGPL570Int = 0;
//                    //if (correlationGPL570Int>200) correlationGPL570Int = 200;
//                    corrDist[correlationInt][correlationGPL570Int]++;
//                }
//            }
//            if (p%100==0) System.out.println(p);
//        }
//
//        long[] dist = new long[201];
//        long[] distGPL570 = new long[201];
//        long sum = 0;
//        for (int d=0; d<201; d++) {
//            for (int e=0; e<201; e++) {
//                long value = corrDist[d][e];
//                dist[d]+=value;
//                distGPL570[e]+=value;
//                sum+=value;
//                if (value > 0) {
//                    double scale = Math.log10((double) value + 1d);
//                    if (scale > 5) {
//                        scale = 5;
//                    }
//                    int colorInt = 255 - (int) Math.round(scale * 51d);
//                    //System.out.println(value + "\t" + scale + "\t" + colorInt);
//                    g2d.setColor(new Color(colorInt, colorInt, colorInt));
//                    g2d.fillRect((int) x0 + d * 2, (int) y1 - e * 2, 2, 2);
//                }
//            }
//        }
//        
//        long[][] corrDistExp = new long[201][201];
//        boolean[][] corrDistSign = new boolean[201][201];
//        for (int d=0; d<201; d++) {
//            for (int e=0; e<201; e++) {
//                corrDistExp[d][e] = (long) ((double) dist[d] * (double) (distGPL570[e] / (double) sum));
//                long value = corrDist[d][e];
//
//                if (value > 0) {
//                    double scale = Math.log10((double) value);
//                    if (scale > 5) {
//                        scale = 5;
//                    }
//                    int colorInt = 255 - (int) Math.round(scale * 51d);
//                    g2d.setColor(new Color(colorInt, colorInt, colorInt));
//                    g2d.setComposite(java.awt.AlphaComposite.getInstance(java.awt.AlphaComposite.SRC_OVER, 1.0f));
//                    g2d.fillRect((int) x0 + d * 2 + 500, (int) y1 - e * 2, 2, 2);
//                }
//
//                double obsDivExp = (double) corrDist[d][e] / (double) corrDistExp[d][e];
//                if (Double.isInfinite(obsDivExp)) obsDivExp = 15;
//
//                double colorRange = 0;
//                if (obsDivExp > 1) {
//                    colorRange = Math.log(obsDivExp) / 2;
//                    if (colorRange > 1d) colorRange = 1d;
//                    if (obsDivExp > 2) {
//                        corrDistSign[d][e] = true;
//                        System.out.println(d + "\t" + e + "\t" + corrDist[d][e] + "\t" + corrDistExp[d][e] + "\t" + obsDivExp);
//                    }
//                    float hue = (float) colorRange;
//                    g2d.setColor(g2d.getColor().getHSBColor(hue, 1.0f, 1.0f));
//                    g2d.setComposite(java.awt.AlphaComposite.getInstance(java.awt.AlphaComposite.SRC_OVER, (float) colorRange));
//                    g2d.fillRect((int) x0 + d * 2 + 500, (int) y1 - e * 2, 2, 2);
//                }
//
//                
//            }
//        }
//
//        g2d.setColor(new Color(0, 0, 0));
//        g2d.setComposite(java.awt.AlphaComposite.getInstance(java.awt.AlphaComposite.SRC_OVER, 0.25f));
//        g2d.drawLine((int) x0, (int) centerY, (int) x1, (int) centerY);
//        g2d.drawLine((int) centerX, (int) y0, (int) centerX, (int) y1);
//        g2d.drawLine((int) x0 + 500, (int) centerY, (int) x1 + 500, (int) centerY);
//        g2d.drawLine((int) centerX + 500, (int) y0, (int) centerX + 500, (int) y1);
//
//        try {
//            //javax.imageio.ImageIO.write(bimage, "png", new File("/Users/lude/Documents/DMG/Data/GPL8490/PCAOverSamplesProbesCenteredScaled/Plots/CorrelationGPL8490-" + nrTCs + "TCs-GPL570-" + nrTCsGPL570 + "TCs.png"));
//            //javax.imageio.ImageIO.write(bimage, "png", new File("/Users/lude/Documents/DMG/Data/GPL8490/PCAOver4432SamplesProbesForcedNormalDistribution/Plots/CorrelationGPL8490-" + nrTCs + "TCs-GPL570-" + nrTCsGPL570 + "TCs.png"));
//            //javax.imageio.ImageIO.write(bimage, "png", new File("/Users/lude/Downloads/PCATCGABeta/CorrelationTCGA-" + nrTCs + "TCs-GPL570-" + nrTCsGPL570 + "TCs.png"));
//            javax.imageio.ImageIO.write(bimage, "png", new File("/Users/lude/Downloads/PCATCGABeta/PCAOverTCGA+GPL8490Combined/CorrelationTCGA+GPL8490-" + nrTCs + "TCs-GPL570-" + nrTCsGPL570 + "TCs.png"));
//        } catch (IOException e) {
//            System.out.println(e.getMessage());
//            e.printStackTrace();
//        }
//        
//
//        for (int p=0; p<dataset.nrProbes; p++) {
//            //if (dataset.probeNames[p].equals("cg19789466")) { //OAS1
//            //if (dataset.probeNames[p].equals("cg19248557")) { //HLA-DRA
//            //if (dataset.probeNames[p].equals("cg14451276")) { //AOAH
//            if (dataset.probeNames[p].equals("cg11519508")) { //TP53
//                for (int q=0; q<dataset.nrProbes; q++) {
//                    if (probeChr[p]!=probeChr[q] && p!=q) {
//                        double covarianceInterim = 0;
//                        for (int s = 0; s < nrSamples; s++) {
//                            covarianceInterim += dataset.rawData[p][s] * dataset.rawData[q][s];
//                        }
//                        double correlation = covarianceInterim / sampleCountMinusOne;
//
//                        double covarianceInterimGPL570 = 0;
//                        for (int s = 0; s < nrSamplesGPL570; s++) {
//                            covarianceInterimGPL570 += datasetEVs.rawData[p][s] * datasetEVs.rawData[q][s];
//                        }
//                        double correlationGPL570 = covarianceInterimGPL570 / sampleCountMinusOneGPL570;
//
//                        int correlationInt = (int) Math.round((correlation + 1d) * 100d);
//                        int correlationGPL570Int = (int) Math.round((correlationGPL570 + 1d) * 100d);
//                        //if (corrDistSign[correlationInt][correlationGPL570Int]) {
//                            //if (correlation < 0 && correlationGPL570 < 0) {
//                        if (correlation > 0.25) {
//                                System.out.println(p + "\t" + dataset.probeNames[p] + "\t" + datasetEVs.probeNames[p] + "\t" + (String) hashProbeToHGNC.get(datasetEVs.probeNames[p]) + "\t" + q + "\t" + dataset.probeNames[q] + "\t" + datasetEVs.probeNames[q] + "\t" + (String) hashProbeToHGNC.get(datasetEVs.probeNames[q]) + "\t" + correlation + "\t" + correlationGPL570);
//                            //}
//                        }
//                    }
//                }
//            }
//            if (p%100==0) System.out.println(p);
//        }
//
//
//
//        System.exit(0);
//    }
//
//    private Jama.EigenvalueDecomposition eigenValueDecomposition(double[][] data) {
//        Jama.Matrix m = new Jama.Matrix(data);
//        Jama.EigenvalueDecomposition eig = m.eig();
//        return eig;
//    }
//
//    private double[] getEigenVector(Jama.EigenvalueDecomposition eig, double[] eigenValues, int pca) {
//        Jama.Matrix eigenValueMatrix = eig.getV();
//        double[][] eigenValueMat = eigenValueMatrix.getArray();
//        double[] eigenVector = new double[eigenValueMat.length];
//        for (int i = 0; i < eigenValueMat.length; i++) {
//            eigenVector[i] = eigenValueMat[i][eigenValueMat.length - 1 - pca]; // * Math.sqrt(eigenValues[eigenValues.length - 1 - pca]);
//        }
//        return eigenVector;
//    }
//
//    private double[] getEigenVector(Jama.EigenvalueDecomposition eig, int pca) {
//        Jama.Matrix eigenValueMatrix = eig.getV();
//        double[][] eigenValueMat = eigenValueMatrix.getArray();
//        double[] eigenVector = new double[eigenValueMat.length];
//        for (int i = 0; i < eigenValueMat.length; i++) {
//            eigenVector[i] = eigenValueMat[i][eigenValueMat.length - 1 - pca]; // * Math.sqrt(eigenValues[eigenValues.length - 1 - pca]);
//        }
//        return eigenVector;
//    }
//
//    private double getEigenValueVar(double[] eigenValues, int pca) {
//        double sumEigenvalues = 0.0;
//        for (Double d : eigenValues) {
//            sumEigenvalues += Math.abs(d);
//        }
//        double result = eigenValues[eigenValues.length - 1 - pca] / sumEigenvalues;
//        return result;
//    }
//
//    private double[] getEigenVectorSVD(Jama.SingularValueDecomposition svd, double[] singularValues, int pca) {
//        Jama.Matrix eigenValueMatrix = svd.getV();
//        double[][] eigenValueMat = eigenValueMatrix.getArray();
//        double[] eigenVector = new double[eigenValueMat.length];
//        for (int i = 0; i < eigenValueMat.length; i++) {
//            eigenVector[i] = eigenValueMat[i][pca] * Math.sqrt(singularValues[pca]);
//        }
//        return eigenVector;
//    }
//
//	//    public void compareSexPhenotype() {
////        
////        
////        int nrTCs = 100;
////        HashMap hashTCsToIncludeHuman = new HashMap();
////        for (int tc=0; tc<nrTCs; tc++) {
////            hashTCsToIncludeHuman.put("Comp" + String.valueOf(tc + 1), null);
////        }
////        
////        ExpressionDataset datasetEVs = new ExpressionDataset("/Users/lude/Documents/DMG/Data/GPL8490/PCAOver4432SamplesProbesCenteredScaled/eigenvectors.txt", "\t", null, hashTCsToIncludeHuman);
////        
////        ExpressionDataset datasetSex = new ExpressionDataset("/Users/lude/Documents/DMG/Data/GPL8490/Gender.txt", "\t", datasetEVs.hashProbes, null);
////        
////        datasetEVs = new ExpressionDataset("/Users/lude/Documents/DMG/Data/GPL8490/PCAOver4432SamplesProbesCenteredScaled/eigenvectors.txt", "\t", datasetSex.hashProbes, hashTCsToIncludeHuman);
////        
////        for (int tc=0; tc<nrTCs; tc++) {
////            Vector vec1 = new Vector();
////            Vector vec2 = new Vector();
////            for (int s=0; s<datasetEVs.nrProbes; s++) { 
////                int sampleIndex = ((Integer) datasetSex.hashProbes.get(datasetEVs.probeNames[s])).intValue();
////                if (datasetSex.rawData[sampleIndex][0]==1) {
////                    vec1.add(datasetEVs.rawData[s][tc]);
////                } else {
////                    vec2.add(datasetEVs.rawData[s][tc]);
////                }
////            }          
////            
////            
////            cern.jet.random.engine.RandomEngine randomEngine = new cern.jet.random.engine.DRand();
////            double[] vals1 = new double[vec1.size()];
////            for (int v = 0; v < vals1.length; v++) {
////                vals1[v] = ((Double) vec1.get(v)).doubleValue();
////            }
////
////            double[] vals2 = new double[vec2.size()];
////            for (int v = 0; v < vals2.length; v++) {
////                vals2[v] = ((Double) vec2.get(v)).doubleValue();
////            }
////            int n1 = vals1.length;
////            int n2 = vals2.length;
////            double mean1 = mean(vals1);
////            double mean2 = mean(vals2);
////            double var1 = variance(vals1, mean1);
////            double var2 = variance(vals2, mean2);
////            double t = (mean1 - mean2) / Math.sqrt(var1 / n1 + var2 / n2);
////            double df = ((var1/n1+var2/n2) * (var1/n1+var2/n2)) / ( ((var1/n1) * (var1/n1)) / (n1-1) + ((var2/n2) * (var2/n2)) / (n2-1));
////            cern.jet.random.StudentT tDist = new cern.jet.random.StudentT(df, randomEngine);
////            double pValue = 1;
////            if (t < 0) {
////                pValue = tDist.cdf(t);
////                if (pValue < 2.0E-323) pValue = 2.0E-323;
////            } else {
////                pValue = tDist.cdf(-t);
////                if (pValue < 2.0E-323) pValue = 2.0E-323;
////            }
////            hgea.math.WilcoxonMannWhitney wmw = new hgea.math.WilcoxonMannWhitney();
////            double wmwPValue = wmw.returnWilcoxonMannWhitneyPValue(vals1, vals2);
////            
////            System.out.println(tc + "\t" + n1 + "\t" + n2 + "\t" + t + "\t" + pValue + "\t" + wmwPValue + "\t" + wmw.getAUC());
////            
////            
////        }
////        
////        System.exit(0);
////        
////    }
//	
//}
