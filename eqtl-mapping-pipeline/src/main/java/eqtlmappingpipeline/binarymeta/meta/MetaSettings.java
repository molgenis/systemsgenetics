/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package eqtlmappingpipeline.binarymeta.meta;

import org.apache.commons.configuration.ConfigurationException;
import org.apache.commons.configuration.XMLConfiguration;

import java.util.ArrayList;
import java.util.logging.Level;
import java.util.logging.Logger;

/**
 * @author harm-jan
 */
public class MetaSettings {

	private int nrPermutations = 10;
	private boolean useAbsoluteZscore = false;
	private int finalEQTLBufferMaxLength = 1000000;
	private int nrOfBins = 100;
	private double fdrthreshold = 0.05;
	private boolean includeSNPsWithoutProperMapping = true;
	private boolean includeProbesWithoutProperMapping = true;
	private boolean cis = true;
	private boolean trans = true;
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
	private int probeDatasetPresenceThreshold = 0;
	private int snpDatasetPresenceThreshold = 0;
	private int probeAndSNPPresenceFilterSampleThreshold = 0;
	private int runonlypermutation;
	private int nrThresds;
	private String probeselection;
	private String snpselection;
	private XMLConfiguration config;
	private String snpprobeselection;

	public void parse(String settings, String texttoreplace, String replacetextwith) {
		try {
			config = new XMLConfiguration(settings);

			nrPermutations = config.getInt("defaults.permutations", 0);

			useAbsoluteZscore = config.getBoolean("defaults.absolutezscore", false);
			finalEQTLBufferMaxLength = config.getInt("defaults.finalnreqtls", 100000);
			fdrthreshold = config.getDouble("defaults.fdrthreshold", 0.05);
			cisdistance = config.getInt("defaults.cisprobedistance", 250000);
			transdistance = config.getInt("defaults.transprobedistance", 5000000);
			includeProbesWithoutProperMapping = config.getBoolean("defaults.includeprobeswithoutmapping", true);
			includeSNPsWithoutProperMapping = config.getBoolean("defaults.includesnpswithoutmapping", true);
			makezscoreplot = config.getBoolean("defaults.makezscoreplot", true);
			makezscoretable = config.getBoolean("defaults.makezscoretable", false);
			probetranslationfile = config.getString("defaults.probetranslationfile");
			String outputStr = config.getString("defaults.output");

			System.out.println("outputstr: " + outputStr);

			if (texttoreplace != null && replacetextwith != null && outputStr.contains(texttoreplace)) {
				outputStr = outputStr.replaceAll(texttoreplace, replacetextwith);
				System.out.println("outputstr: " + outputStr);
			}
			output = outputStr;
			System.out.println("outputstr: " + outputStr);
//			System.exit(-1);


			probeDatasetPresenceThreshold = config.getInt("defaults.minimalnumberofdatasetsthatcontainprobe", 0);
			snpDatasetPresenceThreshold = config.getInt("defaults.minimalnumberofdatasetsthatcontainsnp", 0);
			probeAndSNPPresenceFilterSampleThreshold = config.getInt("defaults.snpprobeselectsamplesizethreshold", -1);

			runonlypermutation = config.getInt("defaults.runonlypermutation", -1);
			nrThresds = config.getInt("defaults.threads", 0);
			cis = config.getBoolean("defaults.cis", false);
			trans = config.getBoolean("defaults.trans", false);

			probeselection = config.getString("defaults.probeselection");

			if (probeselection != null && probeselection.trim().length() == 0) {
				probeselection = null;
			}
			snpselection = config.getString("defaults.snpselection");

			if (snpselection != null && snpselection.trim().length() == 0) {
				snpselection = null;
			}

			if (texttoreplace != null && replacetextwith != null && snpselection.contains(texttoreplace)) {
				snpselection = snpselection.replaceAll(texttoreplace, replacetextwith);
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
			e.printStackTrace();
		}
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
	 * @return the fdrthreshold
	 */
	public double getFdrthreshold() {
		return fdrthreshold;
	}

	/**
	 * @param fdrthreshold the fdrthreshold to set
	 */
	public void setFdrthreshold(double fdrthreshold) {
		this.fdrthreshold = fdrthreshold;
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
	 * @return the cis
	 */
	public boolean isCis() {
		return cis;
	}

	/**
	 * @param cis the cis to set
	 */
	public void setCis(boolean cis) {
		this.cis = cis;
	}

	/**
	 * @return the trans
	 */
	public boolean isTrans() {
		return trans;
	}

	/**
	 * @param trans the trans to set
	 */
	public void setTrans(boolean trans) {
		this.trans = trans;
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
	public int getNrThresds() {
		return nrThresds;
	}

	/**
	 * @param nrThresds the nrThresds to set
	 */
	public void setNrThresds(int nrThresds) {
		this.nrThresds = nrThresds;
	}

	ArrayList<String> getDatasetPrefix() {
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
			Logger.getLogger(MetaSettings.class.getName()).log(Level.SEVERE, null, ex);
		}

	}
}

/*
 * String[] locations = { //
 * "/Volumes/Data2/MarjoleinHomeAccount/marjolein/Results/2011-09-30-EGCUT-Trans-40PCs-GeneticVectorsNotRemoved-SNPs4PCVectorsNotRemoved/", //
 * "/Volumes/Data2/MarjoleinHomeAccount/marjolein/Results/2011-10-20.RotterdamStudy.TRANS.40PCs.YES-GVR.YES-GWAS-PCs/",
 * "/Volumes/Data2/MarjoleinHomeAccount/marjolein/Results/2011-10-25-Groningen-BloodH8v2-TRANS-40PCs-4GWASPCs-GeneticPCsNotRemoved/", //
 * "/Volumes/Data2/MarjoleinHomeAccount/marjolein/Results/2011-10-25-Groningen-BloodHT12-TRANS-40PCs-4GWASPCs-GeneticPCsNotRemoved/", //
 * "/Volumes/Data2/MarjoleinHomeAccount/marjolein/Results/2011-10-20_SHIP-TREND_CIS-TRANS_40PCs_YES-GVR_YES-GWASPCs/", // //
 * "/Volumes/Data2/MarjoleinHomeAccount/marjolein/BloodHT12-40PC-GeneticVectorsNotRemoved/", //
 * "/Volumes/Data2/MarjoleinHomeAccount/marjolein/BloodH8-40PC-GeneticVectorsNotRemoved/" };
 *
 * String[] datasetnames = { "Dataset", //Rotterdam "Dataset", //SHIP "Dataset", //Rotterdam "Dataset", //SHIP "Dataset" //SHIP // "Dataset", //EST //
 * "BloodH8v2-Imputed", //HT12 // "BloodHT12-Imputed" //H8 }; // String[] annotationused = { // "HumanHT-12_V3_0_R2_11283641_A.txt", //
 * "HumanHT-12_V4_0_R1_15002873_B.txt", "H8v2ConvToHT12", // "HumanHT-12_V3_0_R2_11283641_A.txt", // "HumanHT-12_V3_0_R2_11283641_A.txt" ////
 * "HumanHT-12_V3_0_R2_11283641_A.txt", // "H8v2ConvToHT12", // "HumanHT-12_V3_0_R2_11283641_A.txt" };
 *
 * // // HumanHT-12_V3_0_R2_11283641_A.txt	HumanHT-12_V3_0_R3_11283641_A.txt	HumanHT-12_V4_0_R1_15002873_B.txt	HumanHT-12_V4_0_R2_15002873_B.txt
 * HumanHT-12_V4_0_R2_15002873_B_WGDASL.txt // // HumanRef-8_V2_0_R4_11223162_A.txt	HUMANREF-8_V3_0_R1_11282963_A_WGDASL.txt // //
 * HumanRef-8_V3_0_R2_11282963_A.txt	HumanRef-8_V3_0_R3_11282963_A.txt	HumanWG-6_V2_0_R4_11223189_A.txt	HumanWG-6_V3_0_R2_11282955_A.txt
 * HumanWG-6_V3_0_R3_11282955_A.txt	H8v2ConvToHT12 // String output =
 * "/Volumes/Data2/MarjoleinHomeAccount/marjolein/Results/2011-10-25-METATEST.GroningenH8v2-40PCs.YES-GVR.YES-GWAS-PCS/"; // //
 */
