/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package umcg.genetica.io.trityper;

import gnu.trove.map.hash.THashMap;
import gnu.trove.map.hash.TObjectIntHashMap;
import gnu.trove.set.hash.THashSet;

import java.io.IOException;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.HashSet;
import java.util.List;

import org.apache.commons.collections.primitives.ArrayIntList;
import umcg.genetica.containers.Pair;
import umcg.genetica.io.Gpio;
import umcg.genetica.io.text.TextFile;
import umcg.genetica.io.trityper.util.ChrAnnotation;
import umcg.genetica.math.stats.Descriptives;
import umcg.genetica.util.RankArray;

/**
 * @author harmjan // rows: individuals cols: probes
 */
public class TriTyperExpressionData {
	
	//Potentialy we can use TObjectIntMaps
	private int[] chrStart;
	private int[] chrStop;
	private byte[] chr;
	private String[] annotation;
	private double[][] matrix;
	private String[] individuals;
	private TObjectIntHashMap<String> individualNameToId;
	private String[] probes;
	private TObjectIntHashMap<String> probeNameToId;
	private THashMap<String, ArrayIntList> annotationToProbeId;
	private THashSet<String> includeIndividuals;
	private boolean confineToProbesThatMapToAnyChromosome;
	private Integer confineToProbesMappingOnChromosome;
	private THashSet<String> probesConfine;
	private double[] probeOriginalMean;
	private double[] probeOriginalVariance;
	private double[] probeMean;
	private double[] probeVariance;
	private String m_platform;
	private Pair<List<String>, List<List<String>>> pathwayDefinitions;
	
	public TriTyperExpressionData() {
	}
	
	public TriTyperExpressionData(String loc, String probeAnnotationLoc, String platform, boolean cistrans) throws IOException {
		
		this.load(loc, probeAnnotationLoc, platform, cistrans);
	}
	
	/**
	 * @return the rowNames
	 */
	public String[] getRowNames() {
		return individuals;
	}
	
	/**
	 * @param rowNames the rowNames to set
	 */
	public void setRowNames(String[] rowNames) {
		this.individuals = rowNames;
	}
	
	/**
	 * @return the rowNameToId
	 */
	public TObjectIntHashMap<String> getRowNameToId() {
		return individualNameToId;
	}
	
	/**
	 * @param rowNameToId the rowNameToId to set
	 */
	public void setRowNameToId(TObjectIntHashMap<String> rowNameToId) {
		this.individualNameToId = rowNameToId;
	}
	
	/**
	 * @return the colNames
	 */
	public String[] getColNames() {
		return probes;
	}
	
	/**
	 * @param colNames the colNames to set
	 */
	public void setColNames(String[] colNames) {
		this.probes = colNames;
	}
	
	/**
	 * @return the colNameToId
	 */
	public TObjectIntHashMap<String> getColNameToId() {
		return probeNameToId;
	}
	
	/**
	 * @param colNameToId the colNameToId to set
	 */
	public void setColNameToId(TObjectIntHashMap<String> colNameToId) {
		this.probeNameToId = colNameToId;
	}
	
	/**
	 * @return the matrix
	 */
	public double[][] getMatrix() {
		return matrix;
	}
	
	/**
	 * @param matrix the matrix to set
	 */
	public void setMatrix(double[][] matrix) {
		this.matrix = matrix;
	}
	
	public void setIncludeIndividuals(THashSet<String> includedIndividuals) {
		this.includeIndividuals = includedIndividuals;
	}
	
	public void setConfineToProbesThatMapToAnyChromosome(boolean chr) {
		this.confineToProbesThatMapToAnyChromosome = chr;
	}
	
	public void setConfineToProbesThatMapToChromosome(Integer chr) {
		this.confineToProbesMappingOnChromosome = chr;
	}
	
	public void confineToProbes(THashSet<String> probes) {
		this.probesConfine = probes;
	}
	
	public final boolean load(String loc, String probeAnnotationLoc, String platform, boolean cistrans) throws IOException {
		if (!Gpio.exists(loc)) {
			throw new IOException("! Error loading expression data: " + loc + " does not exist");
		}
		System.out.println("Loading expression data from: " + loc);
		
		TextFile in = new TextFile(loc, TextFile.R);
		
		// detect whether TriTyper dataset
		// load probeAnnotation otherwise
		int numProbes = 0;
		boolean trityperformat = false;
		int offset = 1;
		String[] elems = null;
		numProbes = in.countLines();
		elems = in.readLineElemsReturnReference(TextFile.tab);
		
		// header line...
		// MultipleHits    SequenceIdentity        Chr     ChrStart        ChrEnd
		if (elems[0].trim().toLowerCase().equals("probe") && elems[3].trim().toLowerCase().equals("chr") && elems[4].trim().toLowerCase().equals("chrstart")) {
			System.out.println("! Expression data is still in the old, deprecated TriTyper format.");
			trityperformat = true;
			offset = 9;
		} else if (probeAnnotationLoc == null && !cistrans) {
			throw new IOException("ERROR: Probe annotation is not specified. Please specify probe annotation or provide expression data in TriTyper format!");
		}
		
		// load the probe annotation, if any present
		HashMap<String, Byte> hashProbeChr = null;
		HashMap<String, Integer> hashProbeChrStart = null;
		HashMap<String, Integer> hashProbeChrStop = null;
		HashMap<String, String> hashAnnot = null;
		
		if (probeAnnotationLoc != null) {
			if (trityperformat) {
				System.out.println("! WARNING: overriding probe annotation from TriTyper file with probe annotation!");
				trityperformat = false;
				offset = 9;
			}
			
			System.out.println("Loading probe annotation from: " + probeAnnotationLoc);
			
			hashProbeChr = new HashMap<String, Byte>();
			hashProbeChrStart = new HashMap<String, Integer>();
			hashProbeChrStop = new HashMap<String, Integer>();
			hashAnnot = new HashMap<String, String>();
			TextFile an = new TextFile(probeAnnotationLoc, TextFile.R);
			an.readLineElemsReturnReference(TextFile.tab);
			String[] anelems = an.readLineElemsReturnReference(TextFile.tab);
			boolean filterplatform = true;
			if (m_platform == null) {
				filterplatform = false;
			}
			
			int numAnnotLoaded = 0;
			
			while (anelems != null) {
				if (anelems.length >= 6 && (!filterplatform || (filterplatform && anelems[0].equals(m_platform)))) {
					String probe = new String(anelems[1]).intern();
					Byte bchr = ChrAnnotation.parseChr(anelems[3]);
					
					Integer ichrStart = -1;
					Integer ichrStop = -1;
					try {
						ichrStart = Integer.parseInt(anelems[4]);
						ichrStop = Integer.parseInt(anelems[5]);
					} catch (NumberFormatException e) {
					}
					
					String annot = new String(anelems[2].getBytes("UTF-8"));
					
					hashProbeChr.put(probe, bchr);
					hashProbeChrStart.put(probe, ichrStart);
					hashProbeChrStop.put(probe, ichrStop);
					hashAnnot.put(probe, annot);
					numAnnotLoaded++;
				}
				anelems = an.readLineElemsReturnReference(TextFile.tab);
			}
			
			an.close();
			
			if (numAnnotLoaded == 0) {
				throw new IOException("Error: no probe annotation available for expression platform: " + platform);
			} else {
				System.out.println("Loaded annotation for " + numAnnotLoaded + " probes.");
			}
		}
		
		int numIndsIncluded = 0;
		TObjectIntHashMap<String> indToId = new TObjectIntHashMap<String>(elems.length, 1f, -9);
		boolean[] includeCol = new boolean[elems.length];
		for (int pos = offset; pos < elems.length; pos++) {
			if (includeIndividuals == null || includeIndividuals.contains(elems[pos])) {
				numIndsIncluded++;
				includeCol[pos] = true;
			} else {
				includeCol[pos] = false;
			}
		}
		
		individuals = new String[numIndsIncluded];
		int indNo = 0;
		for (int pos = offset; pos < elems.length; pos++) {
			if (includeCol[pos]) {
				individuals[indNo] = new String(elems[pos].getBytes("UTF-8"));
				indToId.put(elems[pos], indNo);
				indNo++;
			}
		}
		
		individualNameToId = indToId;
		
		System.out.println("Found gene expression data for " + numIndsIncluded + " individuals");
		
		ArrayList<Integer> tmpChrStart = new ArrayList<Integer>();
		ArrayList<Integer> tmpChrEnd = new ArrayList<Integer>();
		ArrayList<Byte> tmpChr = new ArrayList<Byte>();
		ArrayList<String> tmpProbe = new ArrayList<String>();
		ArrayList<float[]> tmpRaw = new ArrayList<float[]>();
		ArrayList<String> tmpAnnotation = new ArrayList<String>();
		
		annotationToProbeId = new THashMap<String, ArrayIntList>();
		int probeNr = 0;
		elems = in.readLineElemsReturnReference(TextFile.tab);
		
		if (!trityperformat) {
			if (hashProbeChr == null && !cistrans) {
				System.out.println("WARNING: no probe annotation available.");
			}
		}
		
		int probesIncluded = 0;
		int probesExcluded = 0;
		HashSet<Integer> probesWithMissingValues = new HashSet<Integer>();
		while (elems != null) {
			boolean printreason = true;
			if (elems.length > 1) {
				String probe = new String(elems[0].getBytes("UTF-8"));
				String annotstr = null;
				Byte probeChr = null;
				Integer probeChrStart = null;
				Integer probeChrEnd = null;
				boolean includeprobe = false;
				String reason = "";
				
				if (trityperformat && elems.length > 9) {
					// Probe	MultipleHits	SequenceIdentity	Chr	ChrStart	ChrEnd	Ensembl	HUGO	DIP
					Byte bchr = ChrAnnotation.parseChr(elems[3]);
					if (probesConfine == null || probesConfine.contains(probe)) {
						if (!confineToProbesThatMapToAnyChromosome || bchr >= 1) {
							if (confineToProbesMappingOnChromosome == null || confineToProbesMappingOnChromosome.equals((int) bchr)) {
								probeChr = bchr;
								probeChrStart = Integer.parseInt(elems[4]);
								probeChrEnd = Integer.parseInt(elems[5]);
								annotstr = new String(elems[7].getBytes("UTF-8"));
								includeprobe = true;
							} else {
								probesExcluded++;
								reason += "\tprobe does not map to chromosome:" + confineToProbesMappingOnChromosome + "\t" + bchr;
							}
							
						} else {
							probesExcluded++;
							reason += "\tprobe does not map to any chromosome:\t" + bchr;
						}
					} else {
						probesExcluded++;
						// reason += "\tprobe is not selected during confinement.";
						printreason = false;
					}
					
				} else if (!trityperformat) {
					if (hashProbeChr != null && !cistrans) {
						probeChr = hashProbeChr.get(probe);
						probeChrStart = hashProbeChrStart.get(probe);
						probeChrEnd = hashProbeChrStop.get(probe);
						annotstr = hashAnnot.get(probe);
						if (probeChr == null) {
							reason += "Probe annotation not loaded for probe:\t" + probe;
						} else if ((probesConfine == null || probesConfine.contains(probe))
								&& (!confineToProbesThatMapToAnyChromosome || probeChr >= 1)
								&& (confineToProbesMappingOnChromosome == null || confineToProbesMappingOnChromosome.equals((int) probeChr))) {
							includeprobe = true;
						} else {
							reason += "Probe not selected during confinement.";
						}
					} else if (!cistrans) {
						throw new IOException("ERROR: probe annotation not loaded?");
					} else if (hashProbeChr != null && cistrans) {
						probeChr = hashProbeChr.get(probe);
						probeChrStart = hashProbeChrStart.get(probe);
						probeChrEnd = hashProbeChrStop.get(probe);
						annotstr = hashAnnot.get(probe);
						includeprobe = true;
					} else {
						if ((probesConfine == null || probesConfine.contains(probe))) {
							includeprobe = true;
							probeChr = -1;
							probeChrStart = -1;
							probeChrEnd = -1;
							annotstr = "-";
						} else {
							printreason = false;
						}
					}
				} else {
					System.err.println("Error: data is in TriTyper format, but is probably wrongly formated (because the number of columns <10!)");
					System.exit(-1);
				}
				
				if (!cistrans && includeprobe && probeChr == null && probeChrStart == null && probeChrEnd == null && !(probesConfine != null && !probesConfine.contains(probe))) {
					includeprobe = false;
					reason += "WARNING: probe\t" + probe + "\thas no annotation at all! Will exclude this probe from further use.";
				}
				
				if (includeprobe) {
					probesIncluded++;
					tmpProbe.add(probe.intern());
					if (probeChr == null) {
						probeChr = -1;
					}
					tmpChr.add(probeChr);
					if (probeChrStart == null) {
						probeChrStart = -1;
					}
					tmpChrStart.add(probeChrStart);
					if (probeChrEnd == null) {
						probeChrEnd = -1;
					}
					tmpChrEnd.add(probeChrEnd);
					
					if (annotstr != null && annotstr.length() > 0) {
						ArrayIntList annotationToProbe = annotationToProbeId.get(new String(elems[7].getBytes("UTF-8")));
						if (annotationToProbe == null) {
							annotationToProbe = new ArrayIntList();
						}
						
						annotationToProbe.add(probeNr);
						annotationToProbeId.put(probe, annotationToProbe);
					}
					if (annotstr == null) {
						annotstr = "-";
					}
					tmpAnnotation.add(annotstr.intern());
					
					int samplePos = 0;
					float[] tmpDt = new float[numIndsIncluded];
					for (int pos = offset; pos < elems.length; pos++) {
						if (includeCol[pos]) {
							try {
								tmpDt[samplePos] = Float.parseFloat(elems[pos]);
							} catch (NumberFormatException e) {
								System.err.println("WARNING: missing value for column:\t" + pos + "\tprobe:\t" + probe);
								tmpDt[samplePos] = Float.NaN;
								probesWithMissingValues.add(probeNr);
							}
							samplePos++;
						}
					}
					
					tmpRaw.add(tmpDt);
					probeNr++;
					if (probeNr % 10000 == 0) {
						System.out.println(probeNr + " probes parsed.");
					}
				} else if (printreason) {
					probesExcluded++;
					//System.out.println("Probe\t" + probe + "\texcluded. Reason:\t" + reason);
				}
			}
			
			elems = in.readLineElemsReturnReference(TextFile.tab);
		}
		in.close();
		
		if (probesWithMissingValues.size() > 0) {
			System.err.println("WARNING: " + probesWithMissingValues.size() + " probes with missing values detected. Will replace missing values with average probe value.");
			ArrayList<float[]> newRaw = new ArrayList<float[]>();
			for (int i = 0; i < tmpRaw.size(); i++) {
				if (probesWithMissingValues.contains(i)) {
					float[] data = tmpRaw.get(i);
					double sum = 0;
					int nrNotMissing = 0;
					for (int j = 0; j < data.length; j++) {
						if (!Float.isNaN(data[j])) {
							sum += data[j];
							nrNotMissing++;
						}
					}
					
					sum /= nrNotMissing;
					
					for (int j = 0; j < data.length; j++) {
						if (Float.isNaN(data[j])) {
							data[j] = (float) sum;
						}
					}
				} else {
					newRaw.add(tmpRaw.get(i));
				}
			}
		}
		
		System.out.println("Probes selected:\t" + probesIncluded);
		System.out.println("Probes not selected:\t" + probesExcluded);
		
		if (probeNr == 0) {
			System.err.println("ERROR: No gene expression data loaded for (no probes matching criteria): " + loc);
			return false;
		}
		
		this.chr = new byte[probeNr];
		this.annotation = new String[probeNr];
		this.chrStart = new int[probeNr];
		this.chrStop = new int[probeNr];
		this.probes = new String[probeNr];
		this.matrix = new double[probeNr][numIndsIncluded];
		
		probeNameToId = new TObjectIntHashMap<String>(probeNr, 1f, -9);
		
		for (int i = 0; i < probeNr; i++) {
			probes[i] = tmpProbe.get(i);
			probeNameToId.put(probes[i], i);
			annotation[i] = tmpAnnotation.get(i);
			if (annotation[i] == null) {
				annotation[i] = "-";
			}
			if (tmpChr.get(i) == null) {
				System.err.println("No chromosome annotation loaded for probe " + tmpProbe.get(i));
				chr[i] = -1;
				chrStart[i] = -1;
				chrStop[i] = -1;
			} else {
				chr[i] = tmpChr.get(i);
				chrStart[i] = tmpChrStart.get(i);
				chrStop[i] = tmpChrEnd.get(i);
			}
			
		}
		
		for (int i = 0; i < numIndsIncluded; i++) {
			for (int p = 0; p < probeNr; p++) {
				matrix[p][i] = tmpRaw.get(p)[i];
			}
		}
		
		tmpChrStart = null;
		tmpChrEnd = null;
		tmpChr = null;
		tmpProbe = null;
		tmpRaw = null;
		tmpAnnotation = null;
		
		if (pathwayDefinitions != null) {
			// do stuff.
			System.out.println("Now performing Pathway merger magic. *stars and unicorns*");
			List<String> pathwayNames = pathwayDefinitions.getLeft();
			List<List<String>> genesInPathways = pathwayDefinitions.getRight();
			double[][] pathwayData = new double[pathwayNames.size()][this.individuals.length];
			boolean[] usePathway = new boolean[pathwayNames.size()];
			int nrPathwaysToUse = 0;
			for (int p = 0, len = pathwayNames.size(); p < len; p++) {
				
				List<String> ensGenesForPathway = genesInPathways.get(p);
				HashSet<String> hashGenesForPathway = new HashSet<String>(ensGenesForPathway);
				int nrProbesFound = 0;
				for (int probe = 0; probe < matrix.length; probe++) {
					String annotationForProbe = annotation[probe];
					String[] annStrings = annotationForProbe.split(";");
					boolean hasRightAnnotation = false;
					for (String s : annStrings) {
						if (hashGenesForPathway.contains(s)) {
							hasRightAnnotation = true;
						}
					}
					if (hasRightAnnotation) {
						// do stuff :-) :P
						for (int i = 0; i < individuals.length; i++) {
							pathwayData[p][i] += matrix[probe][i];
						}
						nrProbesFound++;
					}
				}
				
				if (nrProbesFound > 1) {
					for (int i = 0; i < individuals.length; i++) {
						pathwayData[p][i] /= nrProbesFound;
					}
				}
				
				if (nrProbesFound >= 10) {
					usePathway[p] = true;
					nrPathwaysToUse++;
				}
			}
			
			System.out.println("Final number of pathways used in this dataset: " + nrPathwaysToUse);
			
			matrix = new double[nrPathwaysToUse][0];
			probes = new String[nrPathwaysToUse];
			probeNameToId = new TObjectIntHashMap<String>(nrPathwaysToUse, 1f, -9);
			probeMean = new double[nrPathwaysToUse];
			probeVariance = new double[nrPathwaysToUse];
			probeOriginalMean = new double[nrPathwaysToUse];
			probeOriginalVariance = new double[nrPathwaysToUse];
			chr = new byte[nrPathwaysToUse];
			chrStart = new int[nrPathwaysToUse];
			chrStop = new int[nrPathwaysToUse];
			int nrPathwaysUsed = 0;
			for (int p = 0, len = pathwayNames.size(); p < len; p++) {
				if (usePathway[p]) {
					matrix[nrPathwaysUsed] = pathwayData[p];
					probes[nrPathwaysUsed] = pathwayNames.get(p);
					probeNameToId.put(probes[nrPathwaysUsed], nrPathwaysUsed);
					nrPathwaysUsed++;
				}
			}
			
		}
		
		//This is required for the correlation matrix.
		probeMean = new double[probes.length];
		probeOriginalMean = new double[probes.length];
		probeVariance = new double[probes.length];
		probeOriginalVariance = new double[probes.length];

//        calcMean();
//        
//        for (int p = 0; p < matrix.length; p++) {
//            for (int s = 0; s < matrix[p].length; s++) {
//                matrix[p][s] -= probeMean[p];
//            }
//        }
//        calcMeanAndVariance();
		setVarianceAndMean();
		System.out.println("Loaded " + matrix.length + " probes for " + individuals.length + " individuals");
		
		return true;
	}
	
	/**
	 * @return the chrStart
	 */
	public int[] getChrStart() {
		return chrStart;
	}
	
	/**
	 * @param chrStart the chrStart to set
	 */
	public void setChrStart(int[] chrStart) {
		this.chrStart = chrStart;
	}
	
	/**
	 * @return the chrStop
	 */
	public int[] getChrStop() {
		return chrStop;
	}
	
	/**
	 * @param chrStop the chrStop to set
	 */
	public void setChrStop(int[] chrStop) {
		this.chrStop = chrStop;
	}
	
	/**
	 * @return the annotation
	 */
	public String[] getAnnotation() {
		return annotation;
	}
	
	/**
	 * @param annotation the annotation to set
	 */
	public void setAnnotation(String[] annotation) {
		this.annotation = annotation;
	}
	
	/**
	 * @return the chr
	 */
	public byte[] getChr() {
		return chr;
	}
	
	/**
	 * @param chr the chr to set
	 */
	public void setChr(byte[] chr) {
		this.chr = chr;
	}
	
	public String[] getIndividuals() {
		return individuals;
	}
	
	public String[] getProbes() {
		return probes;
	}
	
	public TObjectIntHashMap<String> getIndividualToId() {
		return individualNameToId;
	}
	
	public TObjectIntHashMap<String> getProbeToId() {
		return probeNameToId;
	}
	
	public Integer getIndividualId(String key) {
		return individualNameToId.get(key);
	}
	
	public void calcMeanAndVariance() {
		//Precalculate means and variances. If data is ranked afterwards this does not satisfy the mean and variance set here.
		//This will improve calculations substantially:
		
		//Shouldn't we actually store the old results of mean and variance in the Original equivalents?
		for (int f = 0; f < probes.length; ++f) {
			probeOriginalMean[f] = Descriptives.mean(getProbeData(f));
			probeMean[f] = probeOriginalMean[f];
			probeOriginalVariance[f] = Descriptives.variance(getProbeData(f), probeMean[f]);
			probeVariance[f] = probeOriginalVariance[f];
		}
	}
	
	public void calcAndSubtractMean() {
		for (int p = 0; p < probes.length; ++p) {
			probeMean[p] = Descriptives.mean(getProbeData(p));
			for (int s = 0; s < matrix[p].length; s++) {
				matrix[p][s] -= probeMean[p];
			}
		}
	}
	
	private void setVarianceAndMean() {
		for (int p = 0; p < probes.length; ++p) {
			probeMean[p] = Descriptives.mean(getProbeData(p));
			probeVariance[p] = Descriptives.variance(getProbeData(p), probeMean[p]);
		}
	}
	
	private double[] getProbeData(int f) {
		double[] probeData = new double[individuals.length];
		System.arraycopy(matrix[f], 0, probeData, 0, individuals.length);
		return probeData;
	}
	
	public double[] getProbeVariance() {
		return probeVariance;
	}
	
	public double[] getProbeMean() {
		return probeMean;
	}
	
	public double[] getOriginalProbeVariance() {
		return probeOriginalVariance;
	}
	
	public double[] getOriginalProbeMean() {
		return probeOriginalMean;
	}
	
	/**
	 * Ranks all the expressiondata available in rawData
	 */
	public void rankAllExpressionData(boolean rankWithTies) {
		
		RankArray r = new RankArray();
//        setVarianceAndMean();
		
		for (int p = 0; p < probes.length; ++p) {
			double[] probeData = getProbeData(p);
			
			if (probeVariance[p] == 0) {
				System.out.println("Excluding probe that has no variance in expression:\t" + probes[p] + "\t" + annotation[p]);
			} else {
				if (Double.isNaN(probeMean[p]) || Double.isNaN(probeVariance[p])) {
					System.out.println("Error ranking expression data: probe mean or variance is NaN!:\t" + p + "\t" + probes[p] + "\t" + probeMean[p] + "\t" + probeVariance[p]);
					System.exit(-1);
				} else {
					probeData = r.rank(probeData, rankWithTies);
					setProbeData(p, probeData);
					probeMean[p] = Descriptives.mean(probeData);
					probeVariance[p] = Descriptives.variance(probeData, probeMean[p]);
				}
			}
		}
	}
	
	private void setProbeData(int f, double[] probeData) {
		for (int i = 0; i < individuals.length; i++) {
			matrix[f][i] = (float) probeData[i];
		}
	}
	
	public void pruneAndReorderSamples(List<String> colObjects) {
		double[][] newMatrix = new double[matrix.length][colObjects.size()];
		TObjectIntHashMap<String> indToInd = new TObjectIntHashMap<String>(colObjects.size(), 1f, -9);
		for (int i = 0; i < colObjects.size(); i++) {
			String ind = colObjects.get(i);
			indToInd.put(ind, i);
		}
		
		for (int i = 0; i < individuals.length; i++) {
			Integer newId = indToInd.get(individuals[i]);
			if (newId != -9) {
				for (int row = 0; row < matrix.length; row++) {
					newMatrix[row][newId] = matrix[row][i];
				}
			}
		}
		
		matrix = newMatrix;
		individualNameToId = indToInd;
		individuals = colObjects.toArray(new String[0]);
	}
	
	public void setPathwayDefinitions(Pair<List<String>, List<List<String>>> pathwayDefinitions) {
		this.pathwayDefinitions = pathwayDefinitions;
	}
}
