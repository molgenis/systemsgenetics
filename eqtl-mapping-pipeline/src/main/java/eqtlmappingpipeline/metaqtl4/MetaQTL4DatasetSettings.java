/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package eqtlmappingpipeline.metaqtl4;

import org.molgenis.genotype.RandomAccessGenotypeDataReaderFormats;

/**
 *
 * @author harmjan
 */
public class MetaQTL4DatasetSettings {

    private final String name;
    private final String genotypeLocation;
    private final String traitDataLocation;
    private final String traitPlatform;
    private final String genotypeToTraitCoupling;

    public MetaQTL4DatasetSettings(String name, String genotypeDataLocation, String traitDataLocation, String genotypeToTraitCoupling, String traitplatform) {
        this.name = name;
        this.genotypeLocation = genotypeDataLocation;
        this.traitDataLocation = traitDataLocation;
        this.genotypeToTraitCoupling = genotypeToTraitCoupling;
        this.traitPlatform = traitplatform;
    }
  
    public String getName() {
        return name;
    }

    public String getGenotypeLocation() {
        return genotypeLocation;
    }

    public String getTraitDataLocation() {
        return traitDataLocation;
    }

    public String getAnnotationFile() {
        return traitPlatform;
    }

    public String getGenotypeToTraitCoupling() {
        return genotypeToTraitCoupling;
    }
    
    public RandomAccessGenotypeDataReaderFormats getGenotypeFormat(){
        return RandomAccessGenotypeDataReaderFormats.TRITYPER;
    }

    public String getTraitPlatform() {
        return traitPlatform;
    }
}
