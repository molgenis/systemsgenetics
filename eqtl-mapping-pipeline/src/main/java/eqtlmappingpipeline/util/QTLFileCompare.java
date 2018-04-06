/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package eqtlmappingpipeline.util;

import eqtlmappingpipeline.binarymeta.meta.graphics.ZScorePlot;
import gnu.trove.map.hash.THashMap;
import gnu.trove.set.hash.THashSet;

import java.io.IOException;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.HashSet;
import java.util.LinkedHashSet;
import java.util.Map;
import java.util.logging.Level;
import java.util.logging.Logger;
import java.util.regex.Pattern;

import umcg.genetica.console.ConsoleGUIElems;
import umcg.genetica.io.text.TextFile;
import umcg.genetica.io.trityper.util.BaseAnnot;

/**
 * @author harmjan
 */
public class QTLFileCompare {
	
	private static Pattern SPLIT_ON_TAB = Pattern.compile("\\t");
	private static Pattern SEMI_COLON_PATTERN = Pattern.compile(";");
	private int nrShared = 0;
	private int nrOpposite = 0;
	
	public int getNrShared() {
		return nrShared;
	}
	
	public int getNrOpposite() {
		return nrOpposite;
	}
	
	public QTLFileCompare() {
	}
	
	public QTLFileCompare(String[] args) {
		
		String out = null;
		String file1 = null;
		String file2 = null;
		double fdrCut = -1;
		double pCut = -1;
		double pCut2 = 1;
		boolean matchOnGeneName = false;
		boolean matchSnpOnPos = false;
		boolean splitGeneNames = false;
		String name1 = null;
		String name2 = null;
		
		for (int i = 0; i < args.length; i++) {
			String arg = args[i];
			String val = null;
			
			if (i + 1 < args.length) {
				val = args[i + 1];
			}
			
			if (arg.equals("--out")) {
				out = val;
			} else if (arg.equals("--file1")) {
				file1 = val;
			} else if (arg.equals("--file2")) {
				file2 = val;
			} else if (arg.equals("--name1")) {
				name1 = val;
			} else if (arg.equals("--name2")) {
				name2 = val;
			} else if (arg.equals("--fdrCutoff")) {
				fdrCut = Double.parseDouble(val);
			} else if (arg.equals("--pCutoff")) {
				pCut = Double.parseDouble(val);
			} else if (arg.equals("--pCutoff2")) {
				pCut2 = Double.parseDouble(val);
			} else if (arg.equals("--genebased")) {
				matchOnGeneName = true;
				System.out.println("Performing gene based analysis");
			} else if (arg.toLowerCase().equals("--matchsnponpos")) {
				matchSnpOnPos = true;
				System.out.println("Matching snp based on position");
			} else if (arg.toLowerCase().equals("--splitgenenames")) {
				splitGeneNames = true;
				System.out.println("Splitting gene names on ;");
			}
		}
		
		if (out != null && file1 != null && file2 != null) {
			try {
				if (pCut == -1) {
					compareOverlapAndZScoreDirectionTwoEQTLFiles(file1, name1, file2, name2, out, fdrCut, matchOnGeneName, matchSnpOnPos, splitGeneNames);
				} else {
					compareOverlapAndZScoreDirectionTwoEQTLFilesMj(file1, file2, out, pCut, pCut2, matchOnGeneName, matchSnpOnPos, splitGeneNames);
				}
			} catch (IOException ex) {
				Logger.getLogger(QTLFileCompare.class.getName()).log(Level.SEVERE, null, ex);
			} catch (Exception ex) {
				Logger.getLogger(QTLFileCompare.class.getName()).log(Level.SEVERE, null, ex);
			}
			
		} else {
			printUsage();
		}
	}
	
	private void printUsage() {
		System.out.print("QTL File comparison\n" + ConsoleGUIElems.LINE);
		System.out.println("Compares two eQTL files with each other.");
		
		System.out.print("Command line options:\n" + ConsoleGUIElems.LINE);
		
		System.out.println("--out\t\tstring\t\tOutput file name\n"
				+ "--file1\t\tstring\t\tLocation of file 1\n"
				+ "--file2\t\tstring\t\tLocation of file 2\n"
				+ "--fdrCutoff\t\tdouble\t\talternative FDR cutoff\n"
				+ "--genebased\t\t\tPerform comparison on the basis of gene names (optional, defaults to probe based comparison)\n"
				+ "--matchSnpOnPos\t\tUse chr and and chr pos to match SNPs and ignore identifiers\n"
				+ "--splitGeneNames\t\tSplit gene names on ; when doing --genebased. Count as 2 effects (beta)");
	}
	
	public final void compareOverlapAndZScoreDirectionTwoEQTLFiles(String file1, String file2, String outputFile, boolean matchOnGeneName) throws IOException, Exception {
		compareOverlapAndZScoreDirectionTwoEQTLFiles(file1, file2, outputFile, -1.0d, matchOnGeneName, false, false);
	}
	
	public final void compareOverlapAndZScoreDirectionTwoEQTLFiles(String file1, String file2, String outputFile, double fdrCut, boolean matchOnGeneName, boolean matchSnpOnPos, boolean splitGeneNames) throws IOException, Exception {
		this.compareOverlapAndZScoreDirectionTwoEQTLFiles(file1, null, file2, null, outputFile, fdrCut, matchOnGeneName, matchSnpOnPos, splitGeneNames);
	}
	
	public final void compareOverlapAndZScoreDirectionTwoEQTLFiles(String file1, String datasetname1, String file2, String datasetname2, String outputFile, double fdrCut, boolean matchOnGeneName, boolean matchSnpOnPos, boolean splitGeneNames) throws IOException, Exception {
		
		double filterOnFDR = fdrCut; //Do we want to use another FDR measure? When set to -1 this is not used at all.
		
		HashMap<String, String> hashConvertProbeNames = new HashMap<String, String>(); //When comparing two eQTL files, run on different platforms, we can convert the probe names from one platform to the other, accommodating this comparison, example: hashConvertProbeNames.put(probeNameInFile1, equivalentProbeNameInFile2);
		HashSet<String> hashExcludeEQTLs = new HashSet<String>();   //We can exclude some eQTLs from the analysis. If requested, put the entire eQTL string in this HashMap for each eQTL. Does not work in combination with mathcing based on chr and pos
		HashSet<String> hashConfineAnalysisToSubsetOfProbes = new HashSet<String>(); //We can confine the analysis to only a subset of probes. If requested put the probe name in this HapMap
		HashSet<String> hashTestedSNPsThatPassedQC = null; //We can confine the analysis to only those eQTLs for which the SNP has been successfully passed QC, otherwise sometimes unfair comparisons are made. If requested, put the SNP name in this HashMap
		
		//Now load the eQTLs for file 1:
		THashMap<String, String[]> hashEQTLs = new THashMap<String, String[]>();
		THashSet<String> hashUniqueProbes = new THashSet<String>();
		THashSet<String> hashUniqueGenes = new THashSet<String>();
		
		TextFile in = new TextFile(file1, TextFile.R);
		String[] header = in.readLineElems(TextFile.tab);
		int fdrcol = 0;
		for (int i = 0; i < header.length; i++) {
			if (header[i].toLowerCase().equals("fdr")) {
				fdrcol = i;
			}
		}
		System.out.println("FDR col file 1: " + fdrcol);
		String[] data = in.readLineElemsReturnReference(SPLIT_ON_TAB);
		
		if (data.length < 5) {
			throw new IllegalStateException("QTL File does not have enough columns. Detected columns: " + data.length + " in file " + in.getFileName());
		}
		
		while (data != null) {
			if (hashConvertProbeNames.size() > 0) {
				if (hashConvertProbeNames.containsKey(data[4].trim())) {
					data[4] = hashConvertProbeNames.get(data[4].trim());
				}
			}
			if (filterOnFDR == -1 || Double.parseDouble(data[fdrcol]) <= filterOnFDR) {
				if (hashConfineAnalysisToSubsetOfProbes.isEmpty() || hashConfineAnalysisToSubsetOfProbes.contains(data[4])) {
					if (matchOnGeneName) {
						if (data[16].length() > 1) {
							if (splitGeneNames) {
								for (String gene : SEMI_COLON_PATTERN.split(data[16])) {
									if (!hashExcludeEQTLs.contains(data[1] + "\t" + data[16])) {
										hashEQTLs.put((matchSnpOnPos ? data[2] + ":" + data[3] : data[1]) + "\t" + gene, data);
										hashUniqueProbes.add(data[4]);
										hashUniqueGenes.add(gene);
									}
									
								}
							} else if (!hashExcludeEQTLs.contains(data[1] + "\t" + data[16])) {
								hashEQTLs.put((matchSnpOnPos ? data[2] + ":" + data[3] : data[1]) + "\t" + data[16], data);
								hashUniqueProbes.add(data[4]);
								hashUniqueGenes.add(data[16]);
								//log.write("Added eQTL from original file " + (matchSnpOnPos ? data[2] + ":" + data[3] : data[1]) + "\t" + data[16]);
							}
							
						}
					} else if (!hashExcludeEQTLs.contains(data[1] + "\t" + data[4])) {
						hashEQTLs.put((matchSnpOnPos ? data[2] + ":" + data[3] : data[1]) + "\t" + data[4], data);
						hashUniqueProbes.add(data[4]);
						hashUniqueGenes.add(data[16]);
						//	log.write("Added eQTL from original file " + (matchSnpOnPos ? data[2] + ":" + data[3] : data[1]) + "\t" + data[4]);
					}
				}
			}
			data = in.readLineElemsReturnReference(SPLIT_ON_TAB);
		}
		in.close();
		
		int nrUniqueProbes = hashUniqueProbes.size();
		int nrUniqueGenes = hashUniqueGenes.size();
		hashUniqueProbes = null;
		hashUniqueGenes = null;
		
		//Initialize Graphics2D for the Z-Score allelic direction comparison:
		int width = 1000;
		int height = 1000;
		int margin = 100;
		int x0 = margin;
		int x1 = width - margin;
		int y0 = margin;
		int y1 = height - margin;
		
		ZScorePlot zs = new ZScorePlot();
		
		
		String zsOutFileName = outputFile + "-ZScoreComparison.png";
		
		if (datasetname1 == null) {
			datasetname1 = "Dataset1";
		}
		if (datasetname2 == null) {
			datasetname1 = "Dataset2";
		}
		zs.init(2, new String[]{datasetname1, datasetname2}, false, zsOutFileName);
		
		//Variables holding variousStatistics:
		int nreQTLsIdenticalDirection = 0;
		int nreQTLsOppositeDirection = 0;
		HashMap<String, Integer> hashEQTLNrTimesAssessed = new HashMap<String, Integer>();
		ArrayList<String> vecEQTLNrTimesAssessed = new ArrayList<String>();
		
		HashMap<String, String[]> hashEQTLs2 = new HashMap<String, String[]>();
		HashSet<String> hashUniqueProbes2 = new HashSet<String>();
		HashSet<String> hashUniqueGenes2 = new HashSet<String>();
		HashSet<String> hashUniqueProbesOverlap = new HashSet<String>();
		HashSet<String> hashUniqueGenesOverlap = new HashSet<String>();
		
		int counterFile2 = 0;
		int overlap = 0;
		ArrayList<Double> vecX = new ArrayList<Double>();
		ArrayList<Double> vecY = new ArrayList<Double>();
		
		//Vector holding all opposite allelic effects:
		LinkedHashSet<String> vecOppositeEQTLs = new LinkedHashSet<String>();
		//Vector holding identifiers observed.
		THashSet<String> identifiersUsed = new THashSet<String>();
		
		//Now process file 2:
		in = new TextFile(file2, TextFile.R);
		header = in.readLineElems(TextFile.tab);
		fdrcol = 0;
		for (int i = 0; i < header.length; i++) {
			if (header[i].toLowerCase().equals("fdr")) {
				fdrcol = i;
			}
		}
		System.out.println("FDR col file2: " + fdrcol);
		data = null;
		TextFile identicalOut = new TextFile(outputFile + "-eQTLsWithIdenticalDirecton.txt.gz", TextFile.W);
		TextFile log = new TextFile(outputFile + "-eQTLComparisonLog.txt", TextFile.W);
		while ((data = in.readLineElemsReturnReference(SPLIT_ON_TAB)) != null) {
			
			if (filterOnFDR == -1 || Double.parseDouble(data[fdrcol]) <= filterOnFDR) {
				
				if (hashConvertProbeNames.size() > 0) {
					if (hashConvertProbeNames.containsKey(data[4].trim())) {
						data[4] = hashConvertProbeNames.get(data[4].trim());
					}
				}
				if (hashConfineAnalysisToSubsetOfProbes.isEmpty() || hashConfineAnalysisToSubsetOfProbes.contains(data[4])) {
					if (matchOnGeneName) {
						if (!hashExcludeEQTLs.contains(data[1] + "\t" + data[16])) {
							if (data[16].length() > 1) {
								
								if (splitGeneNames) {
									for (String gene : SEMI_COLON_PATTERN.split(data[16])) {
										
										hashUniqueProbes2.add(data[4]);
										hashUniqueGenes2.add(gene);
										if (!hashEQTLs2.containsKey((matchSnpOnPos ? data[2] + ":" + data[3] : data[1]) + "\t" + gene)) {
											hashEQTLs2.put((matchSnpOnPos ? data[2] + ":" + data[3] : data[1]) + "\t" + gene, data);
											counterFile2++;
										}
										
									}
								} else {
									
									hashUniqueProbes2.add(data[4]);
									hashUniqueGenes2.add(data[16]);
									if (!hashEQTLs2.containsKey((matchSnpOnPos ? data[2] + ":" + data[3] : data[1]) + "\t" + data[16])) {
										hashEQTLs2.put((matchSnpOnPos ? data[2] + ":" + data[3] : data[1]) + "\t" + data[16], data);
										counterFile2++;
									}
								}
							}
						}
					} else if (!hashExcludeEQTLs.contains(data[1] + "\t" + data[4])) {
						//hashEQTLs2.put(data[1] + "\t" + data[4], str);
						hashUniqueProbes2.add(data[4]);
						hashUniqueGenes2.add(data[16]);
						counterFile2++;
					}
					String[] eQTL = null;
					String identifier = null;
					if (matchOnGeneName) {
						
						if (data.length > 16 && data[16].length() > 1) {
							if (splitGeneNames) {
								//NB Plotting and processing of all QTLs here is not okay!
								for (String gene : SEMI_COLON_PATTERN.split(data[16])) {
									if (!hashExcludeEQTLs.contains(data[1] + "\t" + gene)) {
										identifier = (matchSnpOnPos ? data[2] + ":" + data[3] : data[1]) + "\t" + gene;
										if (hashEQTLs.containsKey(identifier)) {
											eQTL = hashEQTLs.get(identifier);
										}
									}
								}
							} else if (!hashExcludeEQTLs.contains(data[1] + "\t" + data[16])) {
								identifier = (matchSnpOnPos ? data[2] + ":" + data[3] : data[1]) + "\t" + data[16];
								if (hashEQTLs.containsKey(identifier)) {
									eQTL = hashEQTLs.get(identifier);
								}
							}
						}
					} else if (!hashExcludeEQTLs.contains(data[1] + "\t" + data[4])) {
						identifier = (matchSnpOnPos ? data[2] + ":" + data[3] : data[1]) + "\t" + data[4];
						if (hashEQTLs.containsKey(identifier)) {
							eQTL = hashEQTLs.get(identifier);
						}
					}
					
					if (eQTL == null) {
						
						//The eQTL, present in file 2 is not present in file 1:
						//if (Double.parseDouble(data[0] < 1E-4) {
						if (hashTestedSNPsThatPassedQC == null || hashTestedSNPsThatPassedQC.contains(data[1])) {
							log.write("eQTL Present In New file But Not In Original File:\t" + identifier + "\t" + data[0] + "\t" + data[2] + "\t" + data[3] + "\t" + data[16] + "\n");
						}
						//}
						double zScore2 = Double.parseDouble(data[10]);
						int posX = 500 + (int) 0;
						int posY = 500 - (int) Math.round(zScore2 * 10);
						zs.draw(null, zScore2, 0, 1);
						
					} else {
						identifiersUsed.add(identifier);
						String[] eQtlData = eQTL;
						boolean identicalProbe = true;
						String probe = data[4];
						String probeFound = eQtlData[4];
						if (!probe.equals(probeFound)) {
							identicalProbe = false;
						}
						
						hashUniqueProbesOverlap.add(data[4]);
						hashUniqueGenesOverlap.add(data[16]);
						if (!hashEQTLNrTimesAssessed.containsKey(identifier)) {
							hashEQTLNrTimesAssessed.put(identifier, 1);
							vecEQTLNrTimesAssessed.add(identifier);
						} else {
							hashEQTLNrTimesAssessed.put(identifier, 1 + hashEQTLNrTimesAssessed.get(identifier));
						}
						String alleles = eQtlData[8];
						String alleleAssessed = eQtlData[9];
						
						String correlations[] = (eQtlData[17]).split(";");
						double correlation = 0;
						int numCorr1 = 0;
						for (int c = 0; c < correlations.length; c++) {
							try {
								if (!correlations[c].equals("-")) {
									correlation += Double.parseDouble(correlations[c]);
									numCorr1++;
								}
							} catch (Exception e) {
							}
						}
						
						correlation /= (double) numCorr1;
//                       if(numCorr1 == 0){
//                           System.out.println("Warning: no correlations defined for eqtl file 1");
//                       }
						double zScore = Double.parseDouble(eQtlData[10]);
//                        double pValue = Double.parseDouble(eQtlData[0]);
						String alleles2 = data[8];
						String alleleAssessed2 = data[9];
						double zScore2 = Double.parseDouble(data[10]);

//                        double pValue2 = Double.parseDouble(data[0]);
						String correlations2[] = data[17].split(";");
						double correlation2 = 0;
						
						boolean alleleflipped = false;
						if (!alleleAssessed.equals(data[9])) {
							if (data[9].equals(eQtlData[8].split("/")[0])) {
								alleleflipped = true;
							} else {
//                               System.out.println("WTF BBQ!");
							}
						}
						
						int numCorr2 = 0;
						for (int c = 0; c < correlations2.length; c++) {
							try {
								if (!correlations2[c].equals("-")) {
									
									correlation2 += (Double.parseDouble(correlations2[c]));
									
									numCorr2++;
								}
							} catch (NumberFormatException e) {
							}
						}
//                       if(numCorr2 == 0){
//                           System.out.println("Warning: no correlations defined for eqtl file 2");
//                       }
						correlation2 /= (double) numCorr2;
						if (alleleflipped) {
							correlation2 = -correlation2;
						}
						boolean sameDirection = false;
						int nrIdenticalAlleles = 0;
						if (alleles.length() > 2 && alleles2.length() > 2) {
							for (int a = 0; a < 3; a++) {
								for (int b = 0; b < 3; b++) {
									if (a != 1 && b != 1) {
										if (alleles.getBytes()[a] == alleles2.getBytes()[b]) {
											nrIdenticalAlleles++;
										}
									}
								}
							}
						}
						
						if (nrIdenticalAlleles == 0) {
							alleles2 = (char) BaseAnnot.getComplement((byte) alleles2.charAt(0)) + "/" + (char) BaseAnnot.getComplement((byte) alleles2.charAt(2));
							alleleAssessed2 = BaseAnnot.getComplement(alleleAssessed2);
							if (alleles.length() > 2 && alleles2.length() > 2) {
								for (int a = 0; a < 3; a++) {
									for (int b = 0; b < 3; b++) {
										if (a != 1 && b != 1) {
											if (alleles.getBytes()[a] == alleles2.getBytes()[b]) {
												nrIdenticalAlleles++;
											}
										}
									}
								}
							}
						}
						
						if (nrIdenticalAlleles != 2) {
							log.write("Error! SNPs have incompatible alleles!!:\t" + alleles + "\t" + alleles2 + "\t" + identifier + "\n");
						} else {
							overlap++;
							if (!alleleAssessed.equals(alleleAssessed2)) {
								zScore2 = -zScore2;
								//                           correlation2 = -correlation2;
								alleleAssessed2 = alleleAssessed;
							}
							
							//Recode alleles:
							// if contains T, but no A, take complement
							//                        if (alleles.contains("T") && !alleles.contains("A")) {
							//                            alleles = BaseAnnot.getComplement(alleles);
							//                            alleleAssessed = BaseAnnot.getComplement(alleleAssessed);
							//                            alleleAssessed2 = BaseAnnot.getComplement(alleleAssessed2);
							//                        }
							if (zScore2 * zScore > 0) {
								sameDirection = true;
							}
							
							//                       if(correlation != correlation2 && (numCorr1 > 0 && numCorr2 > 0)){
							//                           if(Math.abs(correlation - correlation2) > 0.00001){
							//                               System.out.println("Correlations are different: "+lineno+"\t"+correlation +"\t"+correlation2+"\t"+str);
							//                           }
							//
							//                       }
							zs.draw(zScore, zScore2, 0, 1);
							if (!sameDirection) {
								nreQTLsOppositeDirection++;
								
								String oppositeEQTL;
								
								if (matchOnGeneName) {
									oppositeEQTL = data[1] + "\t" + data[16];
									
								} else {
									oppositeEQTL = data[1] + "\t" + data[4];
								}
								
								oppositeEQTL += '\t' + alleles + '\t' + alleleAssessed + '\t' + zScore + '\t' + alleles2 + '\t' + alleleAssessed2 + '\t' + zScore2;
								
								if (!vecOppositeEQTLs.contains(oppositeEQTL)) {
									vecOppositeEQTLs.add(oppositeEQTL);
								}
								
								//                            int posX = 500 + (int) Math.round(zScore * 10);
								//                            int posY = 500 - (int) Math.round(zScore2 * 10);
								vecX.add(zScore);
								vecY.add(zScore2);
								
							} else {
								// write to output
								
								String identicalEQTL;
								
								if (matchOnGeneName) {
									identicalEQTL = data[1] + "\t" + data[16];
									
								} else {
									identicalEQTL = data[1] + "\t" + data[4];
								}
								
								identicalOut.writeln(identicalEQTL + '\t' + alleles + '\t' + alleleAssessed + '\t' + zScore + '\t' + alleles2 + '\t' + alleleAssessed2 + '\t' + zScore2);
								nreQTLsIdenticalDirection++;
								if (alleles.length() > 2 && !alleles.equals("A/T") && !alleles.equals("T/A") && !alleles.equals("C/G") && !alleles.equals("G/C")) {
									//                                int posX = 500 + (int) Math.round(zScore * 10);
									//                                int posY = 500 - (int) Math.round(zScore2 * 10);
									vecX.add(zScore);
									vecY.add(zScore2);
								}
							}
						}
					}
				}
			}
		}
		
		identicalOut.close();
		in.close();
		
		log.write("\n/// Writing missing QTLs observed in original file but not in the new file ////\n\n");
		for (Map.Entry<String, String[]> QTL : hashEQTLs.entrySet()) {
			if (!identifiersUsed.contains(QTL.getKey())) {
				//The eQTL, present in file 1 is not present in file 2:
				//if (Double.parseDouble(data[0] < 1E-4) {
				if (hashTestedSNPsThatPassedQC == null || hashTestedSNPsThatPassedQC.contains(data[1])) {
					log.write("eQTL Present In Original file But Not In New File:\t" + QTL.getKey() + "\t" + QTL.getValue()[0] + "\t" + QTL.getValue()[2] + "\t" + QTL.getValue()[3] + "\t" + QTL.getValue()[16] + "\n");
				}
//                }
				double zScore = Double.parseDouble(QTL.getValue()[10]);
//                int posX = 500 + (int) 0;
//                int posY = 500 - (int) Math.round(zScore * 10);
				zs.draw(zScore, null, 0, 1);
			}
		}
		
		log.close();
		
		double[] valsX = new double[vecX.size()];
		double[] valsY = new double[vecX.size()];
		for (int v = 0; v < valsX.length; v++) {
			valsX[v] = ((Double) vecX.get(v)).doubleValue();
			valsY[v] = ((Double) vecY.get(v)).doubleValue();
		}
		if (valsX.length > 2) {
			double correlation = JSci.maths.ArrayMath.correlation(valsX, valsY);
			double r2 = correlation * correlation;

            /*
             * randomEngine = new cern.jet.random.engine.DRand();
             tDistColt = new cern.jet.random.StudentT(olsY.length - 4, randomEngine);
             */
			cern.jet.random.tdouble.engine.DoubleRandomEngine randomEngine = new cern.jet.random.tdouble.engine.DRand();
			cern.jet.random.tdouble.StudentT tDistColt = new cern.jet.random.tdouble.StudentT(valsX.length - 2, randomEngine);
			double pValuePearson = 1;
			double tValue = correlation / (Math.sqrt((1 - r2) / (double) (valsX.length - 2)));
			if (tValue < 0) {
				pValuePearson = tDistColt.cdf(tValue);
			} else {
				pValuePearson = tDistColt.cdf(-tValue);
			}
			pValuePearson *= 2;
			System.out.println("\nCorrelation between the Z-Scores of the overlapping set of eQTLs:\t" + correlation + "\tP-Value:\t" + pValuePearson);
		}
		
		TextFile out = new TextFile(outputFile + "-OppositeEQTLs.txt", TextFile.W);
		for (String oppositeEQTL : vecOppositeEQTLs) {
			out.write(oppositeEQTL);
			out.append('\n');
		}
		out.close();
		
		in.close();
		
		zs.write(zsOutFileName);
		
		TextFile outSummary = new TextFile(outputFile + "-Summary.txt", TextFile.W);
		
		System.out.println("");
		System.out.println("Nr of eQTLs:\t" + hashEQTLs.size() + "\tin file:\t" + file1 + "\tNrUniqueProbes:\t" + nrUniqueProbes + "\tNrUniqueGenes:\t" + nrUniqueGenes);
		outSummary.writeln("Nr of eQTLs:\t" + hashEQTLs.size() + "\tin file:\t" + file1 + "\tNrUniqueProbes:\t" + nrUniqueProbes + "\tNrUniqueGenes:\t" + nrUniqueGenes);
		
		System.out.println("Nr of eQTLs:\t" + counterFile2 + "\tin file:\t" + file2 + "\tNrUniqueProbes:\t" + hashUniqueProbes2.size() + "\tNrUniqueGenes:\t" + hashUniqueGenes2.size());
		outSummary.writeln("Nr of eQTLs:\t" + counterFile2 + "\tin file:\t" + file2 + "\tNrUniqueProbes:\t" + hashUniqueProbes2.size() + "\tNrUniqueGenes:\t" + hashUniqueGenes2.size());
		
		System.out.println("Overlap:\t" + overlap + "\tNrUniqueProbesOverlap:\t" + hashUniqueProbesOverlap.size() + "\tNrUniqueGenesOverlap:\t" + hashUniqueGenesOverlap.size());
		outSummary.writeln("Overlap:\t" + overlap + "\tNrUniqueProbesOverlap:\t" + hashUniqueProbesOverlap.size() + "\tNrUniqueGenesOverlap:\t" + hashUniqueGenesOverlap.size());
		
		System.out.println("");
		outSummary.writeln();
		
		System.out.println("Nr eQTLs with identical direction:\t" + nreQTLsIdenticalDirection);
		outSummary.writeln("Nr eQTLs with identical direction:\t" + nreQTLsIdenticalDirection);
		
		double proportionOppositeDirection = 100d * (double) nreQTLsOppositeDirection / (double) (nreQTLsOppositeDirection + nreQTLsIdenticalDirection);
		String proportionOppositeDirectionString = (new java.text.DecimalFormat("0.00;-0.00", new java.text.DecimalFormatSymbols(java.util.Locale.US))).format(proportionOppositeDirection);
		
		System.out.println("Nr eQTLs with opposite direction:\t" + nreQTLsOppositeDirection + "\t(" + proportionOppositeDirectionString + "%)");
		outSummary.writeln("Nr eQTLs with opposite direction:\t" + nreQTLsOppositeDirection + "\t(" + proportionOppositeDirectionString + "%)");
		
		outSummary.close();
		
		nrShared = hashUniqueProbesOverlap.size();
		nrOpposite = nreQTLsOppositeDirection;
		
	}
	
	public final void compareOverlapAndZScoreDirectionTwoEQTLFilesMj(String file1, String file2, String outputFile, double pCut, double pCut2, boolean matchOnGeneName, boolean matchSnpOnPos, boolean splitGeneNames) throws IOException, Exception {
		
		HashMap<String, String> hashConvertProbeNames = new HashMap<String, String>(); //When comparing two eQTL files, run on different platforms, we can convert the probe names from one platform to the other, accommodating this comparison, example: hashConvertProbeNames.put(probeNameInFile1, equivalentProbeNameInFile2);
		HashSet<String> hashExcludeEQTLs = new HashSet<String>();   //We can exclude some eQTLs from the analysis. If requested, put the entire eQTL string in this HashMap for each eQTL. Does not work in combination with mathcing based on chr and pos
		HashSet<String> hashConfineAnalysisToSubsetOfProbes = new HashSet<String>(); //We can confine the analysis to only a subset of probes. If requested put the probe name in this HapMap
		HashSet<String> hashTestedSNPsThatPassedQC = null; //We can confine the analysis to only those eQTLs for which the SNP has been successfully passed QC, otherwise sometimes unfair comparisons are made. If requested, put the SNP name in this HashMap
		
		//Now load the eQTLs for file 1:
		THashMap<String, String[]> hashEQTLs = new THashMap<String, String[]>();
		THashSet<String> hashUniqueProbes = new THashSet<String>();
		THashSet<String> hashUniqueGenes = new THashSet<String>();
		
		TextFile in = new TextFile(file1, TextFile.R);
		in.readLine();
		String[] data = in.readLineElemsReturnReference(SPLIT_ON_TAB);
		
		if (data.length < 5) {
			throw new IllegalStateException("QTL File does not have enough columns. Detected columns: " + data.length + " in file " + in.getFileName());
		}
		
		while (data != null) {
			if (hashConvertProbeNames.size() > 0) {
				if (hashConvertProbeNames.containsKey(data[4].trim())) {
					data[4] = hashConvertProbeNames.get(data[4].trim());
				}
			}
			if (Double.parseDouble(data[0]) <= pCut) {
				if (hashConfineAnalysisToSubsetOfProbes.isEmpty() || hashConfineAnalysisToSubsetOfProbes.contains(data[4])) {
					if (matchOnGeneName) {
						if (data[16].length() > 1) {
							
							if (splitGeneNames) {
								for (String gene : SEMI_COLON_PATTERN.split(data[16])) {
									if (!hashExcludeEQTLs.contains(data[1] + "\t" + data[16])) {
										hashEQTLs.put((matchSnpOnPos ? data[2] + ":" + data[3] : data[1]) + "\t" + gene, data);
										hashUniqueProbes.add(data[4]);
										hashUniqueGenes.add(gene);
									}
									
								}
							} else if (!hashExcludeEQTLs.contains(data[1] + "\t" + data[16])) {
								hashEQTLs.put((matchSnpOnPos ? data[2] + ":" + data[3] : data[1]) + "\t" + data[16], data);
								hashUniqueProbes.add(data[4]);
								hashUniqueGenes.add(data[16]);
								//log.write("Added eQTL from original file " + (matchSnpOnPos ? data[2] + ":" + data[3] : data[1]) + "\t" + data[16]);
							}
							
						}
					} else if (!hashExcludeEQTLs.contains(data[1] + "\t" + data[4])) {
						hashEQTLs.put((matchSnpOnPos ? data[2] + ":" + data[3] : data[1]) + "\t" + data[4], data);
						hashUniqueProbes.add(data[4]);
						hashUniqueGenes.add(data[16]);
						//	log.write("Added eQTL from original file " + (matchSnpOnPos ? data[2] + ":" + data[3] : data[1]) + "\t" + data[4]);
					}
				}
			}
			data = in.readLineElemsReturnReference(SPLIT_ON_TAB);
		}
		in.close();
		
		int nrUniqueProbes = hashUniqueProbes.size();
		int nrUniqueGenes = hashUniqueGenes.size();
		hashUniqueProbes = null;
		hashUniqueGenes = null;
		
		//Initialize Graphics2D for the Z-Score allelic direction comparison:
		int width = 1000;
		int height = 1000;
		int margin = 100;
		int x0 = margin;
		int x1 = width - margin;
		int y0 = margin;
		int y1 = height - margin;
		
		ZScorePlot zs = new ZScorePlot();
		String zsOutFileName = outputFile + "-ZScoreComparison.pdf";
		zs.init(2, new String[]{"Dataset1", "Dataset2"}, true, zsOutFileName);
		
		//Variables holding variousStatistics:
		int nreQTLsIdenticalDirection = 0;
		int nreQTLsOppositeDirection = 0;
		HashMap<String, Integer> hashEQTLNrTimesAssessed = new HashMap<String, Integer>();
		ArrayList<String> vecEQTLNrTimesAssessed = new ArrayList<String>();
		
		HashMap<String, String[]> hashEQTLs2 = new HashMap<String, String[]>();
		HashSet<String> hashUniqueProbes2 = new HashSet<String>();
		HashSet<String> hashUniqueGenes2 = new HashSet<String>();
		HashSet<String> hashUniqueProbesOverlap = new HashSet<String>();
		HashSet<String> hashUniqueGenesOverlap = new HashSet<String>();
		
		int counterFile2 = 0;
		int overlap = 0;
		ArrayList<Double> vecX = new ArrayList<Double>();
		ArrayList<Double> vecY = new ArrayList<Double>();
		
		//Vector holding all opposite allelic effects:
		LinkedHashSet<String> vecOppositeEQTLs = new LinkedHashSet<String>();
		//Vector holding identifiers observed.
		THashSet<String> identifiersUsed = new THashSet<String>();
		
		//Now process file 2:
		in = new TextFile(file2, TextFile.R);
		in.readLine();
		
		data = null;
		TextFile comparisonOut = new TextFile(outputFile + "-compared.txt", TextFile.W);
		TextFile log = new TextFile(outputFile + "-eQTLComparisonLog.txt", TextFile.W);
		while ((data = in.readLineElemsReturnReference(SPLIT_ON_TAB)) != null) {
			if (Double.parseDouble(data[0]) <= pCut2) {
				if (hashConvertProbeNames.size() > 0) {
					if (hashConvertProbeNames.containsKey(data[4].trim())) {
						data[4] = hashConvertProbeNames.get(data[4].trim());
					}
				}
				if (hashConfineAnalysisToSubsetOfProbes.isEmpty() || hashConfineAnalysisToSubsetOfProbes.contains(data[4])) {
					if (matchOnGeneName) {
						if (!hashExcludeEQTLs.contains(data[1] + "\t" + data[16])) {
							if (data[16].length() > 1) {
								
								if (splitGeneNames) {
									for (String gene : SEMI_COLON_PATTERN.split(data[16])) {
										
										hashUniqueProbes2.add(data[4]);
										hashUniqueGenes2.add(gene);
										if (!hashEQTLs2.containsKey((matchSnpOnPos ? data[2] + ":" + data[3] : data[1]) + "\t" + gene)) {
											hashEQTLs2.put((matchSnpOnPos ? data[2] + ":" + data[3] : data[1]) + "\t" + gene, data);
											counterFile2++;
										}
										
									}
								} else {
									
									hashUniqueProbes2.add(data[4]);
									hashUniqueGenes2.add(data[16]);
									if (!hashEQTLs2.containsKey((matchSnpOnPos ? data[2] + ":" + data[3] : data[1]) + "\t" + data[16])) {
										hashEQTLs2.put((matchSnpOnPos ? data[2] + ":" + data[3] : data[1]) + "\t" + data[16], data);
										counterFile2++;
									}
								}
							}
						}
					} else if (!hashExcludeEQTLs.contains(data[1] + "\t" + data[4])) {
						//hashEQTLs2.put(data[1] + "\t" + data[4], str);
						hashUniqueProbes2.add(data[4]);
						hashUniqueGenes2.add(data[16]);
						counterFile2++;
					}
					String[] eQTL = null;
					String identifier = null;
					if (matchOnGeneName) {
						
						if (data.length > 16 && data[16].length() > 1) {
							if (splitGeneNames) {
								//NB Plotting and processing of all QTLs here is not okay!
								for (String gene : SEMI_COLON_PATTERN.split(data[16])) {
									if (!hashExcludeEQTLs.contains(data[1] + "\t" + gene)) {
										identifier = (matchSnpOnPos ? data[2] + ":" + data[3] : data[1]) + "\t" + gene;
										if (hashEQTLs.containsKey(identifier)) {
											eQTL = hashEQTLs.get(identifier);
										}
									}
								}
							} else if (!hashExcludeEQTLs.contains(data[1] + "\t" + data[16])) {
								identifier = (matchSnpOnPos ? data[2] + ":" + data[3] : data[1]) + "\t" + data[16];
								if (hashEQTLs.containsKey(identifier)) {
									eQTL = hashEQTLs.get(identifier);
								}
							}
						}
					} else if (!hashExcludeEQTLs.contains(data[1] + "\t" + data[4])) {
						identifier = (matchSnpOnPos ? data[2] + ":" + data[3] : data[1]) + "\t" + data[4];
						if (hashEQTLs.containsKey(identifier)) {
							eQTL = hashEQTLs.get(identifier);
						}
					}
					
					if (eQTL == null) {
						
						//The eQTL, present in file 2 is not present in file 1:
						//if (Double.parseDouble(data[0] < 1E-4) {
						if (hashTestedSNPsThatPassedQC == null || hashTestedSNPsThatPassedQC.contains(data[1])) {
							log.write("eQTL Present In New file But Not In Original File:\t" + identifier + "\t" + data[0] + "\t" + data[2] + "\t" + data[3] + "\t" + data[16] + "\n");
						}
						//}
						double zScore2 = Double.parseDouble(data[10]);
						int posX = 500 + (int) 0;
						int posY = 500 - (int) Math.round(zScore2 * 10);
						zs.draw(null, zScore2, 0, 1);
						
					} else {
						identifiersUsed.add(identifier);
						String[] eQtlData = eQTL;
						boolean identicalProbe = true;
						String probe = data[4];
						String probeFound = eQtlData[4];
						if (!probe.equals(probeFound)) {
							identicalProbe = false;
						}
						
						hashUniqueProbesOverlap.add(data[4]);
						hashUniqueGenesOverlap.add(data[16]);
						if (!hashEQTLNrTimesAssessed.containsKey(identifier)) {
							hashEQTLNrTimesAssessed.put(identifier, 1);
							vecEQTLNrTimesAssessed.add(identifier);
						} else {
							hashEQTLNrTimesAssessed.put(identifier, 1 + hashEQTLNrTimesAssessed.get(identifier));
						}
						String alleles = eQtlData[8];
						String alleleAssessed = eQtlData[9];
						
						String correlations[] = (eQtlData[17]).split(";");
						double correlation = 0;
						int numCorr1 = 0;
						for (int c = 0; c < correlations.length; c++) {
							try {
								if (!correlations[c].equals("-")) {
									correlation += Double.parseDouble(correlations[c]);
									numCorr1++;
								}
							} catch (Exception e) {
							}
						}
						
						correlation /= (double) numCorr1;
//                       if(numCorr1 == 0){
//                           System.out.println("Warning: no correlations defined for eqtl file 1");
//                       }
						double zScore = Double.parseDouble(eQtlData[10]);
//                        double pValue = Double.parseDouble(eQtlData[0]);
						String alleles2 = data[8];
						String alleleAssessed2 = data[9];
						double zScore2 = Double.parseDouble(data[10]);

//                        double pValue2 = Double.parseDouble(data[0]);
						String correlations2[] = data[17].split(";");
						double correlation2 = 0;
						
						boolean alleleflipped = false;
						if (!alleleAssessed.equals(data[9])) {
							if (data[9].equals(eQtlData[8].split("/")[0])) {
								alleleflipped = true;
							} else {
//                               System.out.println("WTF BBQ!");
							}
						}
						
						int numCorr2 = 0;
						for (int c = 0; c < correlations2.length; c++) {
							try {
								if (!correlations2[c].equals("-")) {
									
									correlation2 += (Double.parseDouble(correlations2[c]));
									
									numCorr2++;
								}
							} catch (NumberFormatException e) {
							}
						}
//                       if(numCorr2 == 0){
//                           System.out.println("Warning: no correlations defined for eqtl file 2");
//                       }
						correlation2 /= (double) numCorr2;
						if (alleleflipped) {
							correlation2 = -correlation2;
						}
						boolean sameDirection = false;
						int nrIdenticalAlleles = 0;
						if (alleles.length() > 2 && alleles2.length() > 2) {
							for (int a = 0; a < 3; a++) {
								for (int b = 0; b < 3; b++) {
									if (a != 1 && b != 1) {
										if (alleles.getBytes()[a] == alleles2.getBytes()[b]) {
											nrIdenticalAlleles++;
										}
									}
								}
							}
						}
						
						if (nrIdenticalAlleles == 0) {
							alleles2 = (char) BaseAnnot.getComplement((byte) alleles2.charAt(0)) + "/" + (char) BaseAnnot.getComplement((byte) alleles2.charAt(2));
							alleleAssessed2 = BaseAnnot.getComplement(alleleAssessed2);
							if (alleles.length() > 2 && alleles2.length() > 2) {
								for (int a = 0; a < 3; a++) {
									for (int b = 0; b < 3; b++) {
										if (a != 1 && b != 1) {
											if (alleles.getBytes()[a] == alleles2.getBytes()[b]) {
												nrIdenticalAlleles++;
											}
										}
									}
								}
							}
						}
						
						if (nrIdenticalAlleles != 2) {
							log.write("Error! SNPs have incompatible alleles!!:\t" + alleles + "\t" + alleles2 + "\t" + identifier + "\n");
						} else {
							overlap++;
							if (!alleleAssessed.equals(alleleAssessed2)) {
								zScore2 = -zScore2;
								//                           correlation2 = -correlation2;
								alleleAssessed2 = alleleAssessed;
							}
							
							//Recode alleles:
							// if contains T, but no A, take complement
							//                        if (alleles.contains("T") && !alleles.contains("A")) {
							//                            alleles = BaseAnnot.getComplement(alleles);
							//                            alleleAssessed = BaseAnnot.getComplement(alleleAssessed);
							//                            alleleAssessed2 = BaseAnnot.getComplement(alleleAssessed2);
							//                        }
							if (zScore2 * zScore > 0) {
								sameDirection = true;
							}
							
							//                       if(correlation != correlation2 && (numCorr1 > 0 && numCorr2 > 0)){
							//                           if(Math.abs(correlation - correlation2) > 0.00001){
							//                               System.out.println("Correlations are different: "+lineno+"\t"+correlation +"\t"+correlation2+"\t"+str);
							//                           }
							//
							//                       }
							zs.draw(zScore, zScore2, 0, 1);
							if (!sameDirection) {
								nreQTLsOppositeDirection++;
								
								String oppositeEQTL;
								
								if (matchOnGeneName) {
									oppositeEQTL = data[1] + "\t" + data[16];
									
								} else {
									oppositeEQTL = data[1] + "\t" + data[4];
								}
								
								oppositeEQTL += '\t' + alleles + '\t' + alleleAssessed + '\t' + zScore + '\t' + alleles2 + '\t' + alleleAssessed2 + '\t' + zScore2;
								
								if (!vecOppositeEQTLs.contains(oppositeEQTL)) {
									vecOppositeEQTLs.add(oppositeEQTL);
								}
								
								//                            int posX = 500 + (int) Math.round(zScore * 10);
								//                            int posY = 500 - (int) Math.round(zScore2 * 10);
								vecX.add(zScore);
								vecY.add(zScore2);
								
							} else {
								// write to output
								
								String identicalEQTL;
								
								if (matchOnGeneName) {
									identicalEQTL = data[1] + "\t" + data[16];
									
								} else {
									identicalEQTL = data[1] + "\t" + data[4];
								}
								
								comparisonOut.writeln(identicalEQTL + '\t' + alleles + '\t' + alleleAssessed + '\t' + zScore + '\t' + alleles2 + '\t' + alleleAssessed2 + '\t' + zScore2);
								nreQTLsIdenticalDirection++;
								if (alleles.length() > 2 && !alleles.equals("A/T") && !alleles.equals("T/A") && !alleles.equals("C/G") && !alleles.equals("G/C")) {
									//                                int posX = 500 + (int) Math.round(zScore * 10);
									//                                int posY = 500 - (int) Math.round(zScore2 * 10);
									vecX.add(zScore);
									vecY.add(zScore2);
								}
							}
						}
					}
				}
			}
		}
		
		in.close();
		
		log.write("\n/// Writing missing QTLs observed in original file but not in the new file ////\n\n");
		for (Map.Entry<String, String[]> QTL : hashEQTLs.entrySet()) {
			if (!identifiersUsed.contains(QTL.getKey())) {
				//The eQTL, present in file 1 is not present in file 2:
				//if (Double.parseDouble(data[0] < 1E-4) {
				if (hashTestedSNPsThatPassedQC == null || hashTestedSNPsThatPassedQC.contains(data[1])) {
					log.write("eQTL Present In Original file But Not In New File:\t" + QTL.getKey() + "\t" + QTL.getValue()[0] + "\t" + QTL.getValue()[2] + "\t" + QTL.getValue()[3] + "\t" + QTL.getValue()[16] + "\n");
				}
//                }
				double zScore = Double.parseDouble(QTL.getValue()[10]);
//                int posX = 500 + (int) 0;
//                int posY = 500 - (int) Math.round(zScore * 10);
				zs.draw(zScore, null, 0, 1);
			}
		}
		
		log.close();
		
		double[] valsX = new double[vecX.size()];
		double[] valsY = new double[vecX.size()];
		for (int v = 0; v < valsX.length; v++) {
			valsX[v] = ((Double) vecX.get(v)).doubleValue();
			valsY[v] = ((Double) vecY.get(v)).doubleValue();
		}
		if (valsX.length > 2) {
			double correlation = JSci.maths.ArrayMath.correlation(valsX, valsY);
			double r2 = correlation * correlation;

            /*
             * randomEngine = new cern.jet.random.engine.DRand();
             tDistColt = new cern.jet.random.StudentT(olsY.length - 4, randomEngine);
             */
			cern.jet.random.tdouble.engine.DoubleRandomEngine randomEngine = new cern.jet.random.tdouble.engine.DRand();
			cern.jet.random.tdouble.StudentT tDistColt = new cern.jet.random.tdouble.StudentT(valsX.length - 2, randomEngine);
			double pValuePearson = 1;
			double tValue = correlation / (Math.sqrt((1 - r2) / (double) (valsX.length - 2)));
			if (tValue < 0) {
				pValuePearson = tDistColt.cdf(tValue);
			} else {
				pValuePearson = tDistColt.cdf(-tValue);
			}
			pValuePearson *= 2;
			System.out.println("\nCorrelation between the Z-Scores of the overlapping set of eQTLs:\t" + correlation + "\tP-Value:\t" + pValuePearson);
		}
		
		for (String oppositeEQTL : vecOppositeEQTLs) {
			comparisonOut.write(oppositeEQTL);
			comparisonOut.append('\n');
		}
		comparisonOut.close();
		
		in.close();
		
		zs.write(zsOutFileName);
		
		TextFile outSummary = new TextFile(outputFile + "-Summary.txt", TextFile.W);
		
		System.out.println("");
		System.out.println("Nr of eQTLs:\t" + hashEQTLs.size() + "\tin file:\t" + file1 + "\tNrUniqueProbes:\t" + nrUniqueProbes + "\tNrUniqueGenes:\t" + nrUniqueGenes);
		outSummary.writeln("Nr of eQTLs:\t" + hashEQTLs.size() + "\tin file:\t" + file1 + "\tNrUniqueProbes:\t" + nrUniqueProbes + "\tNrUniqueGenes:\t" + nrUniqueGenes);
		
		System.out.println("Nr of eQTLs:\t" + counterFile2 + "\tin file:\t" + file2 + "\tNrUniqueProbes:\t" + hashUniqueProbes2.size() + "\tNrUniqueGenes:\t" + hashUniqueGenes2.size());
		outSummary.writeln("Nr of eQTLs:\t" + counterFile2 + "\tin file:\t" + file2 + "\tNrUniqueProbes:\t" + hashUniqueProbes2.size() + "\tNrUniqueGenes:\t" + hashUniqueGenes2.size());
		
		System.out.println("Overlap:\t" + overlap + "\tNrUniqueProbesOverlap:\t" + hashUniqueProbesOverlap.size() + "\tNrUniqueGenesOverlap:\t" + hashUniqueGenesOverlap.size());
		outSummary.writeln("Overlap:\t" + overlap + "\tNrUniqueProbesOverlap:\t" + hashUniqueProbesOverlap.size() + "\tNrUniqueGenesOverlap:\t" + hashUniqueGenesOverlap.size());
		
		System.out.println("");
		outSummary.writeln();
		
		System.out.println("Nr eQTLs with identical direction:\t" + nreQTLsIdenticalDirection);
		outSummary.writeln("Nr eQTLs with identical direction:\t" + nreQTLsIdenticalDirection);
		
		double proportionOppositeDirection = 100d * (double) nreQTLsOppositeDirection / (double) (nreQTLsOppositeDirection + nreQTLsIdenticalDirection);
		String proportionOppositeDirectionString = (new java.text.DecimalFormat("0.00;-0.00", new java.text.DecimalFormatSymbols(java.util.Locale.US))).format(proportionOppositeDirection);
		
		System.out.println("Nr eQTLs with opposite direction:\t" + nreQTLsOppositeDirection + "\t(" + proportionOppositeDirectionString + "%)");
		outSummary.writeln("Nr eQTLs with opposite direction:\t" + nreQTLsOppositeDirection + "\t(" + proportionOppositeDirectionString + "%)");
		
		outSummary.close();
		
		nrShared = hashUniqueProbesOverlap.size();
		nrOpposite = nreQTLsOppositeDirection;
		
	}
}
