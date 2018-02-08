/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package nl.umcg.westrah.binarymetaanalyzer;

import java.io.IOException;
import java.util.ArrayList;
import java.util.logging.Level;
import java.util.logging.Logger;

import org.apache.commons.configuration.ConfigurationException;
import org.apache.commons.configuration.XMLConfiguration;
import umcg.genetica.io.Gpio;

/**
 * @author harm-jan
 */
public class BinaryMetaAnalysisSettings {
	
	private int nrPermutations = 10;
	private int startPermutations = 00;
	private boolean useAbsoluteZscore = false;
	private int finalEQTLBufferMaxLength = 1000000;
	private int nrOfBins = 100;
	private boolean includeSNPsWithoutProperMapping = true;
	private boolean includeProbesWithoutProperMapping = true;
	
	public Analysis analysisType;
	
	private int cisdistance = 250000;
	private int transdistance = 5000000;
	private boolean makezscoreplot = true;
	private String probetranslationfile;
	private ArrayList<String> datasetnames;
	private ArrayList<String> datasetPrefix;
	private ArrayList<String> datasetlocations;
	private ArrayList<String> datasetannotations;
	private ArrayList<Integer> selectedProbes;
	private String output;
	private boolean makezscoretable = false;
	private boolean confineSNPs = false;
	
	private int probeDatasetPresenceThreshold = 0;
	private int snpDatasetPresenceThreshold = 0;
	private int probeAndSNPPresenceFilterSampleThreshold = 0;
	private int runonlypermutation;
	private int nrThresds;
	private String probeselection;
	private String snpselection;
	private XMLConfiguration config;
	private String snpprobeselection;
	
	private String snpAnnotationFile;
	
	public Analysis getAnalysisType() {
		return analysisType;
	}
	
	public void copyToOutputDir() throws ConfigurationException {
		config.save(output + "settings.xml");
	}
	
	
	public enum Analysis {
		CIS, TRANS, CISTRANS
	}
	
	public void parse(String settings, String texttoreplace, String replacetextwith)  {
		try {
			config = new XMLConfiguration(settings);
			
			nrPermutations = config.getInt("defaults.permutations", 0);
			startPermutations = config.getInt("defaults.startpermutation", 0);
			snpAnnotationFile = config.getString("defaults.snpannotation");
			useAbsoluteZscore = config.getBoolean("defaults.absolutezscore", false);
			finalEQTLBufferMaxLength = config.getInt("defaults.finalnreqtls", 100000);
			cisdistance = config.getInt("defaults.cisprobedistance", 250000);
			transdistance = config.getInt("defaults.transprobedistance", 5000000);
			includeProbesWithoutProperMapping = config.getBoolean("defaults.includeprobeswithoutmapping", true);
			includeSNPsWithoutProperMapping = config.getBoolean("defaults.includesnpswithoutmapping", true);
			makezscoreplot = config.getBoolean("defaults.makezscoreplot", true);
			makezscoretable = config.getBoolean("defaults.makezscoretable", false);
			confineSNPs = config.getBoolean("defaults.confineSNPsToSNPsPresentInAllDatasets", false);
			probetranslationfile = config.getString("defaults.probetranslationfile");
			output = config.getString("defaults.output");
			if (!Gpio.exists(output)) {
				try {
					Gpio.createDir(output);
					copyToOutputDir();
				} catch (IOException e) {
					e.printStackTrace();
				}
			}
			
			probeDatasetPresenceThreshold = config.getInt("defaults.minimalnumberofdatasetsthatcontainprobe", 0);
			snpDatasetPresenceThreshold = config.getInt("defaults.minimalnumberofdatasetsthatcontainsnp", 0);
			probeAndSNPPresenceFilterSampleThreshold = config.getInt("defaults.snpprobeselectsamplesizethreshold", -1);
			
			runonlypermutation = config.getInt("defaults.runonlypermutation", -1);
			nrThresds = config.getInt("defaults.threads", 1);
			boolean cis = config.getBoolean("defaults.cis", false);
			boolean trans = config.getBoolean("defaults.trans", false);
			if (cis && !trans) {
				analysisType = Analysis.CIS;
			} else if (!cis && trans) {
				analysisType = Analysis.TRANS;
			} else {
				analysisType = Analysis.CISTRANS;
			}
			
			probeselection = config.getString("defaults.probeselection");
			
			if (probeselection != null && probeselection.trim().length() == 0) {
				probeselection = null;
			}
			
			snpselection = config.getString("defaults.snpselection");
			
			if (snpselection != null && snpselection.trim().length() == 0) {
				snpselection = null;
			}
			
			snpprobeselection = config.getString("defaults.snpprobeselection");
			
			if (snpprobeselection != null && snpprobeselection.trim().length() == 0) {
				snpprobeselection = null;
			} else {
				System.out.println("SNP PROBE SELECTION: " + snpprobeselection);
			}
			
			int i = 0;
			
			String dataset = "";
			datasetnames = new ArrayList<String>();
			datasetlocations = new ArrayList<String>();
			datasetannotations = new ArrayList<String>();
			datasetPrefix = new ArrayList<String>();
			
			while (dataset != null) {
				dataset = config.getString("datasets.dataset(" + i + ").name");  // see if a dataset is defined
				if (dataset != null) {
					
					datasetnames.add(dataset);
					String prefix = config.getString("datasets.dataset(" + i + ").prefix");  // see if a dataset is defined
					
					if (prefix == null) {
						prefix = "Dataset";
					}
					datasetPrefix.add(prefix);
					String datasetlocation = config.getString("datasets.dataset(" + i + ").location");  // see if a dataset is defined
					if (texttoreplace != null && replacetextwith != null && datasetlocation.contains(texttoreplace)) {
						datasetlocation = datasetlocation.replace(texttoreplace, replacetextwith);
					}
					String datasetannotation = config.getString("datasets.dataset(" + i + ").expressionplatform");  // see if a dataset is defined
					
					datasetlocations.add(datasetlocation);
					datasetannotations.add(datasetannotation);
				}
				i++;
			}
			
			// parse datasets
		} catch (ConfigurationException e) {
		
		}
	}
	
	public int getStartPermutations() {
		return startPermutations;
	}
	
	public void setStartPermutations(int startPermutations) {
		this.startPermutations = startPermutations;
	}
	
	public boolean isConfineSNPs() {
		return confineSNPs;
	}
	
	public void setConfineSNPs(boolean confineSNPs) {
		this.confineSNPs = confineSNPs;
	}
	
	/**
	 * @return the nrPermutations
	 */
	public int getNrPermutations() {
		return nrPermutations;
	}
	
	/**
	 * @param nrPermutations the nrPermutations to set
	 */
	public void setNrPermutations(int nrPermutations) {
		this.nrPermutations = nrPermutations;
	}
	
	/**
	 * @return the useAbsoluteZscore
	 */
	public boolean isUseAbsoluteZscore() {
		return useAbsoluteZscore;
	}
	
	/**
	 * @param useAbsoluteZscore the useAbsoluteZscore to set
	 */
	public void setUseAbsoluteZscore(boolean useAbsoluteZscore) {
		this.useAbsoluteZscore = useAbsoluteZscore;
	}
	
	/**
	 * @return the finalEQTLBufferMaxLength
	 */
	public int getFinalEQTLBufferMaxLength() {
		return finalEQTLBufferMaxLength;
	}
	
	/**
	 * @param finalEQTLBufferMaxLength the finalEQTLBufferMaxLength to set
	 */
	public void setFinalEQTLBufferMaxLength(int finalEQTLBufferMaxLength) {
		this.finalEQTLBufferMaxLength = finalEQTLBufferMaxLength;
	}
	
	/**
	 * @return the nrOfBins
	 */
	public int getNrOfBins() {
		return nrOfBins;
	}
	
	/**
	 * @param nrOfBins the nrOfBins to set
	 */
	public void setNrOfBins(int nrOfBins) {
		this.nrOfBins = nrOfBins;
	}
	
	/**
	 * @return the includeSNPsWithoutProperMapping
	 */
	public boolean isIncludeSNPsWithoutProperMapping() {
		return includeSNPsWithoutProperMapping;
	}
	
	/**
	 * @param includeSNPsWithoutProperMapping the
	 *                                        includeSNPsWithoutProperMapping to set
	 */
	public void setIncludeSNPsWithoutProperMapping(boolean includeSNPsWithoutProperMapping) {
		this.includeSNPsWithoutProperMapping = includeSNPsWithoutProperMapping;
	}
	
	/**
	 * @return the includeProbesWithoutProperMapping
	 */
	public boolean isIncludeProbesWithoutProperMapping() {
		return includeProbesWithoutProperMapping;
	}
	
	/**
	 * @param includeProbesWithoutProperMapping the
	 *                                          includeProbesWithoutProperMapping to set
	 */
	public void setIncludeProbesWithoutProperMapping(boolean includeProbesWithoutProperMapping) {
		this.includeProbesWithoutProperMapping = includeProbesWithoutProperMapping;
	}
	
	/**
	 * @return the cisdistance
	 */
	public int getCisdistance() {
		return cisdistance;
	}
	
	/**
	 * @param cisdistance the cisdistance to set
	 */
	public void setCisdistance(int cisdistance) {
		this.cisdistance = cisdistance;
	}
	
	/**
	 * @return the transdistance
	 */
	public int getTransdistance() {
		return transdistance;
	}
	
	/**
	 * @param transdistance the transdistance to set
	 */
	public void setTransdistance(int transdistance) {
		this.transdistance = transdistance;
	}
	
	/**
	 * @return the makezscoreplot
	 */
	public boolean isMakezscoreplot() {
		return makezscoreplot;
	}
	
	/**
	 * @param makezscoreplot the makezscoreplot to set
	 */
	public void setMakezscoreplot(boolean makezscoreplot) {
		this.makezscoreplot = makezscoreplot;
	}
	
	/**
	 * @return the probetranslationfile
	 */
	public String getProbetranslationfile() {
		return probetranslationfile;
	}
	
	/**
	 * @param probetranslationfile the probetranslationfile to set
	 */
	public void setProbetranslationfile(String probetranslationfile) {
		this.probetranslationfile = probetranslationfile;
	}
	
	/**
	 * @return the datasetnames
	 */
	public ArrayList<String> getDatasetnames() {
		return datasetnames;
	}
	
	/**
	 * @param datasetnames the datasetnames to set
	 */
	public void setDatasetnames(ArrayList<String> datasetnames) {
		this.datasetnames = datasetnames;
	}
	
	/**
	 * @return the datasetlocations
	 */
	public ArrayList<String> getDatasetlocations() {
		return datasetlocations;
	}
	
	/**
	 * @param datasetlocations the datasetlocations to set
	 */
	public void setDatasetlocations(ArrayList<String> datasetlocations) {
		this.datasetlocations = datasetlocations;
	}
	
	/**
	 * @return the datasetannotations
	 */
	public ArrayList<String> getDatasetannotations() {
		return datasetannotations;
	}
	
	/**
	 * @param datasetannotations the datasetannotations to set
	 */
	public void setDatasetannotations(ArrayList<String> datasetannotations) {
		this.datasetannotations = datasetannotations;
	}
	
	/**
	 * @return the output
	 */
	public String getOutput() {
		return output;
	}
	
	/**
	 * @param output the output to set
	 */
	public void setOutput(String output) {
		this.output = output;
	}
	
	/**
	 * @return the makezscoretable
	 */
	public boolean isMakezscoretable() {
		return makezscoretable;
	}
	
	/**
	 * @param makezscoretable the makezscoretable to set
	 */
	public void setMakezscoretable(boolean makezscoretable) {
		this.makezscoretable = makezscoretable;
	}
	
	/**
	 * @return the probeDatasetPresenceThreshold
	 */
	public int getProbeDatasetPresenceThreshold() {
		return probeDatasetPresenceThreshold;
	}
	
	/**
	 * @param probeDatasetPresenceThreshold the probeDatasetPresenceThreshold to
	 *                                      set
	 */
	public void setProbeDatasetPresenceThreshold(int probeDatasetPresenceThreshold) {
		this.probeDatasetPresenceThreshold = probeDatasetPresenceThreshold;
	}
	
	/**
	 * @return the snpDatasetPresenceThreshold
	 */
	public int getSnpDatasetPresenceThreshold() {
		return snpDatasetPresenceThreshold;
	}
	
	/**
	 * @param snpDatasetPresenceThreshold the snpDatasetPresenceThreshold to set
	 */
	public void setSnpDatasetPresenceThreshold(int snpDatasetPresenceThreshold) {
		this.snpDatasetPresenceThreshold = snpDatasetPresenceThreshold;
	}
	
	/**
	 * @return the probeAndSNPPresenceFilterSampleThreshold
	 */
	public int getProbeAndSNPPresenceFilterSampleThreshold() {
		return probeAndSNPPresenceFilterSampleThreshold;
	}
	
	/**
	 * @param probeAndSNPPresenceFilterSampleThreshold the
	 *                                                 probeAndSNPPresenceFilterSampleThreshold to set
	 */
	public void setProbeAndSNPPresenceFilterSampleThreshold(int probeAndSNPPresenceFilterSampleThreshold) {
		this.probeAndSNPPresenceFilterSampleThreshold = probeAndSNPPresenceFilterSampleThreshold;
	}
	
	/**
	 * @return the runonlypermutation
	 */
	public int getRunonlypermutation() {
		return runonlypermutation;
	}
	
	/**
	 * @param runonlypermutation the runonlypermutation to set
	 */
	public void setRunonlypermutation(int runonlypermutation) {
		this.runonlypermutation = runonlypermutation;
	}
	
	/**
	 * @return the nrThresds
	 */
	public int getNrThreads() {
		return nrThresds;
	}
	
	/**
	 * @param nrThresds the nrThresds to set
	 */
	public void setNrThresds(int nrThresds) {
		this.nrThresds = nrThresds;
	}
	
	public ArrayList<String> getDatasetPrefix() {
		return datasetPrefix;
	}
	
	/**
	 * @return the probeselection
	 */
	public String getProbeselection() {
		return probeselection;
	}
	
	/**
	 * @param probeselection the probeselection to set
	 */
	public void setProbeselection(String probeselection) {
		this.probeselection = probeselection;
	}
	
	public String getSNPSelection() {
		return snpselection;
	}
	
	public String getSNPProbeSelection() {
		return snpprobeselection;
	}
	
	void save() {
		try {
			config.save(output + "metasettings.xml");
		} catch (ConfigurationException ex) {
			Logger.getLogger(BinaryMetaAnalysisSettings.class.getName()).log(Level.SEVERE, null, ex);
		}
		
	}
	
	public String getSNPAnnotationFile() {
		return snpAnnotationFile;
	}
}
