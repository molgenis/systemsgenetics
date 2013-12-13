/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package nl.systemsgenetics.eqtlpermutationtranscriptionfactoranalysis;

import java.io.File;
import java.io.IOException;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.Iterator;
import java.util.List;
import java.util.Map;
import java.util.Set;
import java.util.TreeMap;
import java.util.regex.Pattern;
import org.molgenis.genotype.RandomAccessGenotypeData;
import org.molgenis.genotype.sampleFilter.SampleIncludedFilter;
import org.molgenis.genotype.trityper.TriTyperGenotypeData;
import org.molgenis.genotype.variantFilter.VariantIdIncludeFilter;
import umcg.genetica.genomicboundaries.GenomicBoundaries;
import umcg.genetica.io.gtf.GffElement;
import umcg.genetica.io.gtf.GtfReader;
import umcg.genetica.io.regulomedb.RegulomeDbEntry;
import umcg.genetica.io.regulomedb.RegulomeDbFile;
import umcg.genetica.io.regulomedb.RegulomeDbFiles;
import umcg.genetica.io.regulomedb.RegulomeDbSupportingData;
import umcg.genetica.io.text.TextFile;
import umcg.genetica.io.ucsc.UCSCDataObject;
import umcg.genetica.io.ucsc.WigFile;

/**
 *
 * @author Matthieu
 */
public class EnrichmentDataReader {
	private static final Pattern TAB_PATTERN = Pattern.compile("\t");
	
	public EnrichmentDataReader(){
		
	}
	
	public RandomAccessGenotypeData readEQtlGenotypeData(String genotypeData, Set<String> variantIdFilter) throws IOException{
		//Provide a Set<String> containing rsID of all significant eQTLs.
		RandomAccessGenotypeData gonlImputedBloodGenotypeData = new TriTyperGenotypeData( new File(genotypeData), 630000, new VariantIdIncludeFilter(variantIdFilter), new SampleIncludedFilter());
		return gonlImputedBloodGenotypeData;
	}
	
	
	public RandomAccessGenotypeData readEQtlGenotypeDataV2(String genotypeData, Set<String> variantIdFilter) throws IOException{
		//Provide a Set<String> containing rsID of all significant eQTLs.
		RandomAccessGenotypeData gonlImputedBloodGenotypeData = new TriTyperGenotypeData( new File(genotypeData), 630000, null, null);
		return gonlImputedBloodGenotypeData;
	}
	
	
	
	public GenomicBoundaries readRegulomeDbDataV2(ArrayList<RegulomeDbFile> regulomeDbFiles){
		GenomicBoundaries<Object> regulomeDbData = new GenomicBoundaries();
		
		RegulomeDbFiles regulomeDbFilesData = new RegulomeDbFiles(regulomeDbFiles);
		Iterator<RegulomeDbEntry> regulomeDbDataIterator = regulomeDbFilesData.iterator();
		int n = 0;
		
		while(regulomeDbDataIterator.hasNext()){
			RegulomeDbEntry rdbe = regulomeDbDataIterator.next();
			String rdbeChr = rdbe.getChr();
			int rdbePos = rdbe.getChrPos();
			String[] transcriptionFactors = getTranscriptionFactors(rdbe);
			
			if(transcriptionFactors.length > 0){
				regulomeDbData.addBoundary(rdbeChr, rdbePos, rdbePos, transcriptionFactors);
			}
		}
		return regulomeDbData;
	}
	
	private String[] getTranscriptionFactors(RegulomeDbEntry rdbe){
		Map<String, List<RegulomeDbSupportingData>> supportData = rdbe.getSupportData();
		ArrayList<String> tfs = new ArrayList<String>();
		
		Iterator it = supportData.entrySet().iterator();
		while(it.hasNext()){
			Map.Entry pairs = (Map.Entry) it.next();

			List<RegulomeDbSupportingData> aap = (List<RegulomeDbSupportingData>) pairs.getValue();
			for(RegulomeDbSupportingData rdbsd : aap){

				//Check if the annotation is protein_binding.
				if(rdbsd.getSupportClass().equalsIgnoreCase("Protein_Binding")){
					tfs.add(rdbsd.getSupportValue());
				}
			}
		}
		return tfs.toArray( new String[tfs.size()] );
	}
	
	
	
	
	
	public GenomicBoundaries<Object> readEncodeAnnotationData(String encodeFile)throws IOException{
		GenomicBoundaries<Object> encodeBoundaries = new GenomicBoundaries();
		
		GtfReader gtfr = new GtfReader(new File(encodeFile));
		Iterator<GffElement> encodeIterator = gtfr.iterator();
		while(encodeIterator.hasNext()){
			GffElement encodeEntry = encodeIterator.next();
			
			encodeBoundaries.addBoundary(encodeEntry.getSeqname(), encodeEntry.getStart(), encodeEntry.getEnd(),
					encodeEntry.getAttributeValue("transcript_type"));
		}
		gtfr.close();
		return encodeBoundaries;
	}
	
	
	
	public GenomicBoundaries<Object> readRepeatData(String repeatFile)throws IOException{
		GenomicBoundaries<Object> repeatRegions = new GenomicBoundaries();
		
		String fileLine;
		String[] fileLineData;
		
		TextFile repeatReader = new TextFile(repeatFile, false);
		while( (fileLine = repeatReader.readLine())!=null ){
			fileLineData = TAB_PATTERN.split(fileLine);
			
			String chr = new String(fileLineData[0]);
			int start = Integer.parseInt(fileLineData[1]);
			int stop = Integer.parseInt(fileLineData[2]);
			String repeatClass = new String(fileLineData[3]);
			//String repeatFamily = new String(fileLineData[4]);
			//String repeatSubFamily = new String(fileLineData[5]);
			//String supportData = repeatClass+"|"+repeatFamily+"|"+repeatSubFamily;
			
			repeatRegions.addBoundary(chr, start, stop, repeatClass);
		}
		return repeatRegions;
	}
	
	
	
	
	
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
	
	//This method might not be neede as the hypersensitivity wig data has the exact same format as the histone wig data.
	public GenomicBoundaries readHypersensitivityDataFromText(String dnaseFileLocation)throws IOException{
		GenomicBoundaries<Object> dnaseSiteBoundaries = new GenomicBoundaries();
		
		String fileLine;
		String[] fileLineData;
		TextFile tf = new TextFile(dnaseFileLocation, false);
		while((fileLine=tf.readLine())!=null){
			fileLineData = TAB_PATTERN.split(fileLine);
			String chr = new String(fileLineData[0]);
			int start = Integer.parseInt(fileLineData[1]);
			int stop = Integer.parseInt(fileLineData[2]);
		}
		return dnaseSiteBoundaries;
	}
	
	
	public GenomicBoundaries readEncodeMethylationData(String methylationFileLocation)throws IOException{
		GenomicBoundaries<Object> encodeMethylationBoundaries = new GenomicBoundaries();
		
		return encodeMethylationBoundaries;
	}
	
	
	public GenomicBoundaries readHistoneNarrowPeakFileData(String peakDataFile)throws IOException{
		GenomicBoundaries<Object> peakData = new GenomicBoundaries();
		
		String fileLine;
		String[] fileLineData;
		TextFile tf = new TextFile(peakDataFile, false);
		while((fileLine=tf.readLine())!=null){
			fileLineData = TAB_PATTERN.split(fileLine);
			String chromosome = new String(fileLineData[0]);
			int startPos = Integer.parseInt(fileLineData[1]);
			int stopPos = Integer.parseInt(fileLineData[2]);
			double peakScore = Double.valueOf(fileLineData[6]);
			
			peakData.addBoundary(chromosome, startPos, stopPos, peakScore);
		}
		tf.close();
		return peakData;
	}
	
	
	public HashMap<String, TreeMap<String, String>> readGwasCatalogData(String gwasCatalogLocation)throws IOException{
		HashMap<String, TreeMap<String, String>> gwasData = new HashMap<String, TreeMap<String, String>>();
		TreeMap<String, String> tmp;
		
		String fileLine;
		String[] fileLineData;
		TextFile tf = new TextFile(gwasCatalogLocation, false);
		while((fileLine = tf.readLine())!=null){
			fileLineData = TAB_PATTERN.split(fileLine);
			//fileLineData[11];
			//fileLineData[23];
			//fileLineData[7];
			
			if(gwasData.containsKey(fileLineData[11])){
				tmp = gwasData.get(fileLineData[11]);
				tmp.put(fileLineData[23], fileLineData[7]);
			}
			
			else{
				tmp = new TreeMap<String, String>();
				tmp.put(fileLineData[23], fileLineData[7]);
				gwasData.put(fileLineData[11], tmp);
			}
		}
		tf.close();
		return gwasData;
	}
	
	
	
	public HashMap<String, TreeMap<Integer, String[]>> readRegulomeDbData(ArrayList<RegulomeDbFile> regulomeDbFiles){
		HashMap<String, TreeMap<Integer, String[]>> regulomeDbData = new HashMap<String, TreeMap<Integer, String[]>>();
		TreeMap<Integer, String[]> tmp;
		
		RegulomeDbFiles regulomeDbFilesData = new RegulomeDbFiles(regulomeDbFiles);
		Iterator<RegulomeDbEntry> regulomeDbDataIterator = regulomeDbFilesData.iterator();
		int n = 0;
		
		while(regulomeDbDataIterator.hasNext()){
			RegulomeDbEntry rdbe = regulomeDbDataIterator.next();
			String rdbeChr = rdbe.getChr();
			int rdbePos = rdbe.getChrPos();
			String[] transcriptionFactors = getTranscriptionFactors(rdbe);
			
			if(transcriptionFactors.length > 0){
				if(regulomeDbData.containsKey(rdbeChr)){
					tmp = regulomeDbData.get(rdbeChr);
					tmp.put(rdbePos, transcriptionFactors);
				}
				else{
					tmp = new TreeMap<Integer, String[]>();
					tmp.put(rdbePos, transcriptionFactors);
					regulomeDbData.put(rdbeChr, tmp);
				}
				n++;
			}
		}
		return regulomeDbData;
	}
}
