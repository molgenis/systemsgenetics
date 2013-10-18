/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package nl.systemsgenetics.eqtlpermutationtranscriptionfactoranalysis;

import java.io.IOException;
import java.util.HashMap;
import java.util.HashSet;
import java.util.Iterator;
import java.util.Map;
import java.util.TreeMap;
import java.util.regex.Pattern;
import org.molgenis.genotype.RandomAccessGenotypeData;
import org.molgenis.genotype.util.Ld;
import org.molgenis.genotype.util.LdCalculatorException;
import org.molgenis.genotype.variant.GeneticVariant;
import umcg.genetica.genomicboundaries.GenomicBoundaries;
import umcg.genetica.genomicboundaries.GenomicBoundary;
import umcg.genetica.graphics.ViolinBoxPlot;
import umcg.genetica.io.text.TextFile;
import umcg.genetica.io.trityper.EQTL;
import umcg.genetica.io.ucsc.UCSCDataObject;
import umcg.genetica.io.ucsc.WigFile;
import umcg.genetica.math.stats.WilcoxonMannWhitney;

/**
 *
 * @author Matthieu
 */
public class HistoneSiteEnrichment {
	private static final Pattern TAB_PATTERN = Pattern.compile("\\t");
	
	public static void main(String[] args) throws IOException{
		HistoneSiteEnrichment hse = new HistoneSiteEnrichment();
		hse.readHistoneDataFromText("C:\\Users\\Matthieu\\Documents\\Afstudeerstage\\Pilot\\3.histones\\Gm12864Ctcf\\Gm12864Ctcfchr21.txt");
	}
	
	/**
	 * Reads histone data from a plain text file into GenomicBoundaries.
	 * @param inputFileLocation
	 * @return
	 * @throws IOException 
	 */
	public GenomicBoundaries readHistoneDataFromText(String inputFileLocation)throws IOException{
		GenomicBoundaries<Object> histoneBoundaries = new GenomicBoundaries();
		
		String fileLine;
		String[] fileLineData;
		TextFile tf = new TextFile(inputFileLocation, false);
		while((fileLine=tf.readLine())!=null){
			fileLineData = TAB_PATTERN.split(fileLine);
			
			String chr = new String(fileLineData[0]);
			int startPos = Integer.parseInt(fileLineData[1]);
			int stopPos = Integer.parseInt(fileLineData[2]);
			double bindingValue = Double.valueOf(fileLineData[3]);
			
			histoneBoundaries.addBoundary(chr, startPos, stopPos, bindingValue);
		}
		tf.close();
		return histoneBoundaries;
	}
	
	
	/**
	 * Reads histone data from a wig file into GenomicBoundaries.
	 * @param wigFileLocation
	 * @return
	 * @throws IOException 
	 */
	public GenomicBoundaries readHistoneSiteDataFromWig(String wigFileLocation) throws IOException{
		GenomicBoundaries<Object> histoneSiteBoundaries = new GenomicBoundaries();
		
		WigFile wf = new WigFile(wigFileLocation, false);
		long totalN = wf.size();
		long n = 0;
		while(n < totalN){
			UCSCDataObject ucscdo = wf.parseLn();
			String chr = Byte.toString(ucscdo.getChr());
			histoneSiteBoundaries.addBoundary(chr, ucscdo.getPositionStart(), ucscdo.getPositionEnd(), ucscdo.getValue());
		}
		wf.close();
		return histoneSiteBoundaries;
	}
	
	
	public void aap(double[] A, double[] B){
		WilcoxonMannWhitney wcmw = new WilcoxonMannWhitney();
		double testPvalue = wcmw.returnWilcoxonMannWhitneyPValue(A, B);
		double wilcoxAuc = wcmw.getAUC();
	}
	
	
	//Draw some boxplots of the data.
	public void drawBoxPlot(double[] eQtlHitScores, double[] permutationHitScores, String fileType, String outputFile)throws IOException{
		//Create some things.
		double[][][] dataToPlot = new double[1][2][0];
		dataToPlot[0][0] = eQtlHitScores;
		dataToPlot[0][1] = permutationHitScores;
		
		String[][] plotLabels = new String[][]{{"Real", "Chance"}};
		String[] datasets = new String[]{"eQTLs", "Permutations"};
		
		//Draw the plot.
		ViolinBoxPlot vbp = new ViolinBoxPlot();
		if( fileType.equalsIgnoreCase("png") ){
			vbp.draw(dataToPlot, datasets, plotLabels, "Binding", ViolinBoxPlot.Output.PNG, outputFile);
		}
		else{
			vbp.draw(dataToPlot, datasets, plotLabels, "Binding", ViolinBoxPlot.Output.PDF, outputFile);
		}
	}
	
	
	/**
	 * Performs the Histone site Enrichment.
	 * @param topEffectData
	 * @param nonTopEffectData
	 * @param countsMap
	 * @param genotypeData
	 * @param boundaries
	 * @param r2CutOff
	 * @throws IOException 
	 */
	public void performAnalysis(HashMap<String, EQTL> topEffectData, HashMap<String, HashSet<Integer>> nonTopEffectData, HashMap<String, Integer> countsMap,
			RandomAccessGenotypeData genotypeData, GenomicBoundaries<Object> boundaries, double r2CutOff)throws IOException{
		//First read the data.
		readHistoneDataFromText("");
		
		Ld ld = null;
		Iterator<Map.Entry<String, EQTL>> topEffectIterator = topEffectData.entrySet().iterator();
		while(topEffectIterator.hasNext()){
			Map.Entry<String, EQTL> topEffectDataEntry = (Map.Entry) topEffectIterator.next();
			EQTL eqtl = topEffectDataEntry.getValue();
			
			//FETCH THE EQTL FROM GENOTYPE DATA.
			String rsChr = eqtl.getRsChr().toString();
			int rsChrPos = eqtl.getRsChrPos();
			String rsProbe = eqtl.getProbe();
			GeneticVariant eQtlVariant = genotypeData.getSnpVariantByPos(rsChr.toString(), rsChrPos);
			
			
			if(eQtlVariant != null){
				HashSet<Integer> nonTopEffectsPos = nonTopEffectData.get(rsProbe);
				
				if(nonTopEffectsPos != null){
					for(int eqtlPos : nonTopEffectsPos){
						GeneticVariant snpVariantByPos = genotypeData.getSnpVariantByPos(rsChr, eqtlPos);

						if(snpVariantByPos != null){
							if(eQtlVariant.isBiallelic() && snpVariantByPos.isBiallelic()){
								try {
									ld = eQtlVariant.calculateLd(snpVariantByPos);
								} catch (LdCalculatorException ex) {
									System.out.println("Error in LD calculation: " + ex.getMessage());
									System.exit(1);
								}

								if(ld.getR2() >= r2CutOff){
									//SEARCH THROUGH REGULOMEDB
									if(boundaries.isInBoundary(rsChr, eqtlPos, 0)){
										GenomicBoundary boundary = boundaries.getBoundary(rsChr, eqtlPos, 0);
										String annotation = (String) boundary.getAnnotation();
										
										if(countsMap.containsKey(annotation)){
											countsMap.put(annotation, countsMap.get(annotation)+1);
										}
										else{
											countsMap.put(annotation, 1);
										}
									}

								}
							}
						}
					}
				}
				
				
				//CODE TO SEARCH FOR THE TOP EFFECT.
				if(boundaries.isInBoundary(rsChr, rsChrPos, 0)){
					GenomicBoundary boundary = boundaries.getBoundary(rsChr, rsChrPos, 0);
					String annotation = (String) boundary.getAnnotation();

					if(countsMap.containsKey(annotation)){
						countsMap.put(annotation, countsMap.get(annotation)+1);
					}
					else{
						countsMap.put(annotation, 1);
					}
				}
			}
		}
	}
	
}
