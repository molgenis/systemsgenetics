package nl.umcg.westrah.binarymetaanalyzer.westrah.binarymetaanalyzer.posthoc;

import umcg.genetica.legacy.ExpressionDataset;

import java.util.HashMap;

public class NullSetGenerator {
	
	public void generateTPAndTNSet() {
	
//		//Directory holding the Z-Score matrices:
//		int nrPermutedZScoreMatrices = 10;
//		String fileDir = "/Users/lude/Documents/Genetica/eQTLGen/ZScoresMatricesMetaAnalysis-20170901/";
//
//		HashMap hashSignificantTransEQTLs = new HashMap();
//		try {
//			java.io.BufferedReader in = new java.io.BufferedReader(new java.io.FileReader(new File("/Users/lude/Downloads/eQTLsFDR0.05-PrunedLevel_2CohortFilter_ParalogueFilter_20171204.txt")));
//			String str = in.readLine();
//			while ((str = in.readLine()) != null) {
//				String[] data = str.split("\t");
//				hashSignificantTransEQTLs.put(data[1] + "\t" + data[4], null);
//			}
//			System.out.println("Number of significant trans-eQTLs in input file:\t" + hashSignificantTransEQTLs.size());
//		} catch (Exception e) {
//			System.out.println("Error:\t" + e.getMessage());
//			e.printStackTrace();
//		}
//
//		System.out.println("Defining TP and TN set of trans-eQTLs:");
//		nrPermutedZScoreMatrices = 1;
//		int[] snpNrSign = new int[20000];
//		int[] geneNrSign = new int[20000];
//		for (int fileNr=0; fileNr<nrPermutedZScoreMatrices + 1; fileNr++) {
//			String fileName = fileDir;
//			String fileNameNrSamples = fileDir;
//			if (fileNr==0) {
//				fileName += "ZScoreMatrix.txt.binary";
//				fileNameNrSamples += "ZScoreMatrixNrSamples.txt.binary";
//			} else {
//				fileName += "ZscoreMatrix-Permutation" + fileNr + ".txt.binary";
//				fileNameNrSamples += "ZScoreMatrixNrSamples.txt.binary";
//			}
//			int[] snpNrSignPerm = new int[20000];
//			int[] geneNrSignPerm = new int[20000];
//			ExpressionDataset dataset = new ExpressionDataset(fileName, "\t", null, null);
//			ExpressionDataset datasetN = new ExpressionDataset(fileNameNrSamples, "\t", null, null);
//			if (fileNr==0) {
//				for (int p=0; p<dataset.nrProbes; p++) {
//					for (int s=0; s<dataset.nrSamples; s++) {
//						if (datasetN.rawData[p][s]>0) {
//							if (Math.abs(dataset.rawData[p][s])>=4.410633) {  //Least significant trans-eQTL in significant trans-eQTL file has Z-Score of 4.4106339.
//								if (hashSignificantTransEQTLs.containsKey(dataset.probeNames[p] + "\t" + dataset.sampleNames[s])) {
//									//System.out.println("TP\t" + dataset.probeNames[p] + "\t" + dataset.sampleNames[s] + "\t" + dataset.rawData[p][s]);
//									snpNrSign[p]++;
//									geneNrSign[s]++;
//								}
//							}
//						}
//					}
//				}
//			} else {
//				ExpressionDataset datasetOriginal = new ExpressionDataset(fileDir + "ZScoreMatrix.txt.binary", "\t", null, null);
//				System.out.println(JSci.maths.ArrayMath.mass(snpNrSign));
//				System.out.println(JSci.maths.ArrayMath.mass(geneNrSign));
//				IntegerDoubleObjectSorter sorter = new IntegerDoubleObjectSorter();
//				Vector vecSNPNrSign = new Vector();
//				for (int p=0; p<dataset.nrProbes; p++) {
//					if (snpNrSign[p]>0) {
//						vecSNPNrSign.add(new IntegerDoubleObject(p, snpNrSign[p]));
//					}
//				}
//				sorter.sort(vecSNPNrSign);
//				for (int a=0; a<vecSNPNrSign.size(); a++) {
//					int p = ((IntegerDoubleObject) vecSNPNrSign.get(vecSNPNrSign.size() - 1 - a)).intValue;
//					Vector vecValues = new Vector();
//					for (int s=0; s<dataset.nrSamples; s++) {
//						if (geneNrSignPerm[s]<geneNrSign[s]) {
//							//Select an equal number of trans-eQTL genes for this gene
//							if (Math.abs(datasetOriginal.rawData[p][s])<2.8) {
//								IntegerDoubleObject object = new IntegerDoubleObject(s, dataset.rawData[p][s]);
//								vecValues.add(object);
//							}
//						}
//					}
//					//System.out.println(p + "\t" + snpNrSign[p] + "\t" + vecValues.size());
//					sorter.sort(vecValues);
//					for (int q=0; q<snpNrSign[p]; q++) {
//						IntegerDoubleObject object = (IntegerDoubleObject) vecValues.get(q);
//						snpNrSignPerm[p]++;
//						geneNrSignPerm[object.intValue]++;
//						System.out.println("TN Permutation " + fileNr + "\t" + dataset.probeNames[p] + "\t" + dataset.sampleNames[object.intValue] + "\t" + datasetOriginal.rawData[p][object.intValue]);
//					}
//				}
//			}
//		}
//		System.exit(0);
	}
	
}
