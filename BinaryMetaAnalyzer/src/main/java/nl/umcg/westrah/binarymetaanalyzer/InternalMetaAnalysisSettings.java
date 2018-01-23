/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package nl.umcg.westrah.binarymetaanalyzer;

import java.util.logging.Level;
import java.util.logging.Logger;
import org.apache.commons.configuration.ConfigurationException;
import org.apache.commons.configuration.XMLConfiguration;

/**
 *
 * @author harm-jan
 */
public class InternalMetaAnalysisSettings {

    private int nrPermutations = 10;
    private int startPermutations = 0;
    private String datasetname;
    private String datasetPrefix;
    private String datasetlocation;
    private String output;
    private String probeToFeature;
    private String zScoreMergeOption;
    
    private XMLConfiguration config;
    
    public void parse(String settings, String texttoreplace, String replacetextwith) {
        try {
            config = new XMLConfiguration(settings);

            nrPermutations = config.getInt("defaults.permutations", 0);
            startPermutations = config.getInt("defaults.startpermutation", 0);
            output = config.getString("defaults.output");
            probeToFeature = config.getString("defaults.probeToFeature");
            datasetname = config.getString("defaults.dataset.name");  // see if a dataset is defined
            datasetPrefix = config.getString("defaults.dataset.prefix");  // see if a dataset is defined
            zScoreMergeOption = config.getString("defaults.zScoreMerge").toLowerCase();
            if (datasetPrefix == null) {
                datasetPrefix = "Dataset";
            }
            String datasetloc = config.getString("defaults.dataset.location");  // see if a dataset is defined
            if (texttoreplace != null && replacetextwith != null && datasetloc.contains(texttoreplace)) {
                datasetloc = datasetloc.replace(texttoreplace, replacetextwith);
            }
            datasetlocation = datasetloc;
            
            // parse datasets
        } catch (ConfigurationException e) {
        
        }
    }

    public String getProbeToFeature() {
        return probeToFeature;
    }

    public void setProbeToFeature(String probeToFeature) {
        this.probeToFeature = probeToFeature;
    }

    public int getStartPermutations() {
        return startPermutations;
    }

    public void setStartPermutations(int startPermutations) {
        this.startPermutations = startPermutations;
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

    public String getzScoreMergeOption() {
        return zScoreMergeOption;
    }

    public void setzScoreMergeOption(String zScoreMergeOption) {
        this.zScoreMergeOption = zScoreMergeOption;
    }

    public XMLConfiguration getConfig() {
        return config;
    }

    public void setConfig(XMLConfiguration config) {
        this.config = config;
    }

    

    /**
     * @return the datasetname
     */
    public String getDatasetname() {
        return datasetname;
    }

    /**
     * @param datasetname the datasetnames to set
     */
    public void setDatasetname(String datasetname) {
        this.datasetname = datasetname;
    }

    /**
     * @return the datasetlocations
     */
    public String getDatasetlocation() {
        return datasetlocation;
    }

    /**
     * @param datasetlocation the datasetlocations to set
     */
    public void setDatasetlocations(String datasetlocation) {
        this.datasetlocation = datasetlocation;
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


    public String getDatasetPrefix() {
        return datasetPrefix;
    }

    
    void save() {
        try {
            config.save(output + "metasettings.xml");
        } catch (ConfigurationException ex) {
            Logger.getLogger(InternalMetaAnalysisSettings.class.getName()).log(Level.SEVERE, null, ex);
        }

    }

}
