/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package org.molgenis.genotype.variant;

import java.util.ArrayList;
import java.util.Collections;
import java.util.LinkedHashSet;
import java.util.List;
import java.util.Map;
import org.molgenis.genotype.Allele;
import org.molgenis.genotype.Alleles;
import org.molgenis.genotype.util.FixedSizeIterable;
import org.molgenis.genotype.util.MafCalculator;
import org.molgenis.genotype.util.ProbabilitiesConvertor;
import org.molgenis.genotype.variant.id.GeneticVariantId;
import org.molgenis.genotype.variant.sampleProvider.SampleVariantsProvider;

/**
 *
 * @author Patrick Deelen
 */
public class ReadOnlyGeneticVariantTriTyper extends AbstractGeneticVariant {

    private final GeneticVariantId variantId;
    private final int startPos;
    private final String sequenceName;
    private final SampleVariantsProvider sampleVariantsProvider;
    private final int indexOfVariantInTriTyperData;
    private Alleles alleles;
	private final GeneticVariantMeta variantMeta;

    public ReadOnlyGeneticVariantTriTyper(String variantId, int startPos, String sequenceName, SampleVariantsProvider sampleVariantsProvider, int indexOfVariantInTriTyperData, GeneticVariantMeta variantMeta) {
        this.variantId = GeneticVariantId.createVariantId(variantId);
        this.startPos = startPos;
        this.sequenceName = sequenceName;
        this.sampleVariantsProvider = sampleVariantsProvider;
        this.indexOfVariantInTriTyperData = indexOfVariantInTriTyperData;
		this.variantMeta = variantMeta;
    }

    @Override
    public String getPrimaryVariantId() {
        return variantId.getPrimairyId();
    }

    @Override
    public List<String> getAlternativeVariantIds() {
        return variantId.getAlternativeIds();
    }

    @Override
    public List<String> getAllIds() {
        return variantId.getVariantIds();
    }

    @Override
    public GeneticVariantId getVariantId() {
        return variantId;
    }

    @Override
    public int getStartPos() {
        return startPos;
    }

    @Override
    public String getSequenceName() {
        return sequenceName;
    }

    @Override
    public final Alleles getVariantAlleles() {
        if (alleles == null) {
            getSampleVariants();
        }
        return alleles;
    }

    @Override
    public int getAlleleCount() {
        return this.getVariantAlleles().getAlleleCount();
    }

    @Override
    public Allele getRefAllele() {
        return this.getVariantAlleles().get(0);
    }

    @Override
    public List<Alleles> getSampleVariants() {
        List<Alleles> sampleVariantAlleles = Collections.unmodifiableList(sampleVariantsProvider.getSampleVariants(this));
        
        if (this.alleles == null) {
            //set alleles here

            LinkedHashSet<Allele> variantAlleles = new LinkedHashSet<Allele>(2);

            for (Alleles alleles2 : sampleVariantAlleles) {
                for (Allele allele : alleles2) {
                    if (allele != allele.ZERO) {
                        variantAlleles.add(allele);
                    }
                }
            }

            this.alleles = Alleles.createAlleles(new ArrayList<Allele>(variantAlleles));


        }
        return sampleVariantAlleles;
    }

    @Override
    public Map<String, ?> getAnnotationValues() {
        return Collections.emptyMap();
    }

    @Override
    public double getMinorAlleleFrequency() {
        return MafCalculator.calculateMaf(this.getVariantAlleles(), this.getRefAllele(), this.getSampleVariants()).getFreq();
    }

    @Override
    public Allele getMinorAllele() {
        return MafCalculator.calculateMaf(this.getVariantAlleles(), this.getRefAllele(), this.getSampleVariants()).getMinorAllele();
    }

    @Override
    public float[] getSampleDosages() {
        return sampleVariantsProvider.getSampleDosage(this);
    }

    @Override
    public SampleVariantsProvider getSampleVariantsProvider() {
        return sampleVariantsProvider;
    }

    @Override
    public byte[] getSampleCalledDosages() {
        return sampleVariantsProvider.getSampleCalledDosage(this);
    }

    @Override
    public List<Boolean> getSamplePhasing() {
        return sampleVariantsProvider.getSamplePhasing(this);
    }

    @Override
    public float[][] getSampleGenotypeProbilities() {
        return sampleVariantsProvider.getSampleProbilities(this);
    }

    @Override
    public double[][] getSampleGenotypeProbabilitiesBgen() {
        throw new UnsupportedOperationException("Not supported yet.");
    }

    @Override
    public int hashCode() {
        // TriTyper genotypes always have a primary ID, we should hash that.
        return this.getPrimaryVariantId().hashCode();
    }

    @Override
    public boolean equals(Object obj) {
        if (this == obj) {
            return true;
        }
        if (!(obj instanceof ReadOnlyGeneticVariantTriTyper)) {
            return false;
        }
        ReadOnlyGeneticVariantTriTyper other = (ReadOnlyGeneticVariantTriTyper) obj;
        
        // in trityper data, all SNPs have an unique identifier.
        // this is not always true for other datatypes, however...
       
        if (this.getPrimaryVariantId() != null && other.getPrimaryVariantId() != null) {
            if (this.getPrimaryVariantId().equals(other.getPrimaryVariantId())) {
                return true;
            } else {
                return false;
            }
        }
        
        if (getSequenceName() == null) {
            if (other.getSequenceName() != null) {
                return false;
            }
        } else if (!getSequenceName().equals(other.getSequenceName())) {
            return false;
        }
        if (getStartPos() != other.getStartPos()) {
            return false;
        }

        // If we get here pos and sequence are identical
        if (getSampleVariantsProvider() == null) {
            if (other.getSampleVariantsProvider() != null) {
                return false;
            }
        } else if (!getSampleVariantsProvider().equals(other.getSampleVariantsProvider())) {
            //For trityper id must be identical
            return this.variantId.getPrimairyId().equals(other.getPrimaryVariantId());
        }
        return true;
    }

    public int getIndexOfVariantInTriTyperData() {
        return indexOfVariantInTriTyperData;
    }
	
	@Override
	public GeneticVariantMeta getVariantMeta()
	{
		return variantMeta;
	}	

	@Override
	public FixedSizeIterable<GenotypeRecord> getSampleGenotypeRecords() {
		return sampleVariantsProvider.getSampleGenotypeRecords(this);
	}
	
	@Override
	public Alleles getAlternativeAlleles() {
		ArrayList<Allele> altAlleles = new ArrayList<>(this.getVariantAlleles().getAlleles());
		altAlleles.remove(this.getRefAllele());
		return Alleles.createAlleles(altAlleles);
	}
	
}
