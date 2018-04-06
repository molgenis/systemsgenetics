package org.molgenis.genotype.modifiable;

import java.io.IOException;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.HashSet;
import java.util.Iterator;
import java.util.List;
import java.util.Map;
import org.molgenis.genotype.AbstractRandomAccessGenotypeData;
import org.molgenis.genotype.Allele;
import org.molgenis.genotype.Alleles;
import org.molgenis.genotype.GenotypeDataException;
import org.molgenis.genotype.RandomAccessGenotypeData;
import org.molgenis.genotype.Sample;
import org.molgenis.genotype.Sequence;
import org.molgenis.genotype.annotation.Annotation;
import org.molgenis.genotype.annotation.SampleAnnotation;
import org.molgenis.genotype.variant.GeneticVariant;
import org.molgenis.genotype.variant.id.GeneticVariantId;
import org.molgenis.genotype.variant.sampleProvider.SampleVariantsProvider;
import org.molgenis.genotype.variant.sampleProvider.SwappingSampleVariantsProvider;

public class ModifiableGenotypeDataInMemory extends AbstractRandomAccessGenotypeData implements ModifiableGenotypeData {

    private final RandomAccessGenotypeData sourceGenotypeData;
    private final HashMap<GeneticVariant, GeneticVariantId> idUpdates;
    private final HashMap<GeneticVariant, Allele> refAlleleUpdate;
    private final HashMap<GeneticVariant, Alleles> allelesUpdate;
    private final HashMap<GeneticVariant, SampleVariantsProvider> variantProviderUpdates;
    private final HashSet<ModifiableGeneticVariant> filteredOutVariants;
    private final HashSet<GeneticVariant> swappedVariants;

    private final HashMap<SampleVariantsProvider, SampleVariantsProvider> swappingSampleVariantProviders;

    public ModifiableGenotypeDataInMemory(RandomAccessGenotypeData sourceGenotypeData) {
        super();
        this.sourceGenotypeData = sourceGenotypeData;
        this.idUpdates = new HashMap<GeneticVariant, GeneticVariantId>();
        this.refAlleleUpdate = new HashMap<GeneticVariant, Allele>();
        this.allelesUpdate = new HashMap<GeneticVariant, Alleles>();
        this.variantProviderUpdates = new HashMap<GeneticVariant, SampleVariantsProvider>();
        this.swappingSampleVariantProviders = new HashMap<SampleVariantsProvider, SampleVariantsProvider>();
        this.filteredOutVariants = new HashSet<ModifiableGeneticVariant>();
        this.swappedVariants = new HashSet<>();
    }

    @Override
    public List<String> getSeqNames() {
        return sourceGenotypeData.getSeqNames();
    }

    @Override
    public Iterable<Sequence> getSequences() {
        return sourceGenotypeData.getSequences();
    }

    @Override
    public Sequence getSequenceByName(String name) {
        return sourceGenotypeData.getSequenceByName(name);
    }

    @Override
    public Iterable<GeneticVariant> getVariantsByPos(String seqName, int startPos) {
        return ModifiableGeneticVariantIterator.createGeneticVariantIterableBackByModifiable(sourceGenotypeData
                .getVariantsByPos(seqName, startPos).iterator(), this, filteredOutVariants);
    }

    @Override
    public GeneticVariant getSnpVariantByPos(String seqName, int startPos) {
        return getModifiableSnpVariantByPos(seqName, startPos);
    }

    @Override
    public Iterable<GeneticVariant> getSequenceGeneticVariants(String seqName) {
        return ModifiableGeneticVariantIterator.createGeneticVariantIterableBackByModifiable(sourceGenotypeData
                .getSequenceGeneticVariants(seqName).iterator(), this, filteredOutVariants);
    }

    @Override
    public List<Annotation> getVariantAnnotations() {
        return sourceGenotypeData.getVariantAnnotations();
    }

    @Override
    public Annotation getVariantAnnotation(String annotationId) {
        return sourceGenotypeData.getVariantAnnotation(annotationId);
    }

    @Override
    public List<Sample> getSamples() {
        return sourceGenotypeData.getSamples();
    }

    @Override
    public Iterator<GeneticVariant> iterator() {
        return ModifiableGeneticVariantIterator.createGeneticVariantIterableBackByModifiable(
                sourceGenotypeData.iterator(), this, filteredOutVariants).iterator();
    }

    @Override
    public synchronized GeneticVariantId getUpdatedId(ModifiableGeneticVariant geneticVariant) {
        return idUpdates.get(geneticVariant.getOriginalVariant());
    }

    @Override
    public synchronized Allele getUpdatedRef(ModifiableGeneticVariant geneticVariant) {
        return refAlleleUpdate.get(geneticVariant.getOriginalVariant());
    }

    @Override
    public synchronized SampleVariantsProvider getUpdatedSampleVariantProvider(ModifiableGeneticVariant geneticVariant) {
        return variantProviderUpdates.get(geneticVariant.getOriginalVariant());
    }

    @Override
    public synchronized void updateVariantId(ModifiableGeneticVariant geneticVariant,
            GeneticVariantId newGeneticVariantId) {

        GeneticVariant originalGeneticVariant = geneticVariant.getOriginalVariant();

        if (originalGeneticVariant.getVariantId().equals(newGeneticVariantId)) {
            idUpdates.remove(originalGeneticVariant);
            return;
        }
        idUpdates.put(originalGeneticVariant, newGeneticVariantId);
    }

    @Override
    public synchronized void updateVariantPrimaryId(ModifiableGeneticVariant geneticVariant, String newPrimaryId) {

        GeneticVariant originalGeneticVariant = geneticVariant.getOriginalVariant();

        GeneticVariantId currentId = idUpdates.get(originalGeneticVariant);

        if (currentId != null) {
            if (newPrimaryId != null && currentId.getPrimairyId().equals(newPrimaryId)) {
                return;
            }
        } else {
            currentId = originalGeneticVariant.getVariantId();
            if (currentId.getPrimairyId() != null && newPrimaryId != null && currentId.getPrimairyId().equals(newPrimaryId)) {
                return;
            }
        }

        String oldPrimairyId = currentId.getPrimairyId();

        // Create alternative alleles based on old alternatives
        ArrayList<String> newAlternativeAlleles = new ArrayList<String>(currentId.getAlternativeIds());
        // Remove new primary if it is one of the old alternative
        newAlternativeAlleles.remove(newPrimaryId);
        // add the old primary ID to the alternative ID list
        newAlternativeAlleles.add(oldPrimairyId);

        updateVariantId(geneticVariant, GeneticVariantId.createVariantId(newPrimaryId, newAlternativeAlleles));

    }

    @Override
    public synchronized void swapGeneticVariant(ModifiableGeneticVariant geneticVariant) {

        GeneticVariant originalGeneticVariant = geneticVariant.getOriginalVariant();

        if (swappedVariants.contains(originalGeneticVariant)) {
            throw new GenotypeDataException("Cannot swap same variant twice");
        }

        swappedVariants.add(originalGeneticVariant);

        Alleles variantAlleles = getUpdatedAlleles(geneticVariant);
        if (variantAlleles == null) {
            variantAlleles = originalGeneticVariant.getVariantAlleles();
        }

        Allele refAllele = getUpdatedRef(geneticVariant);
        if (refAllele == null) {
            refAllele = originalGeneticVariant.getRefAllele();
        }

        SampleVariantsProvider sampleVariantProvider = getUpdatedSampleVariantProvider(geneticVariant);
        if (sampleVariantProvider == null) {
            sampleVariantProvider = originalGeneticVariant.getSampleVariantsProvider();
        }

        SampleVariantsProvider swappingSampleVariantsProvider = swappingSampleVariantProviders
                .get(sampleVariantProvider);
        if (swappingSampleVariantsProvider == null) {
            swappingSampleVariantsProvider = new SwappingSampleVariantsProvider(sampleVariantProvider);
            //Enabling this will cause bug.
//			if (sampleVariantProvider.cacheSize() > 0)
//			{
//				swappingSampleVariantsProvider = new CachedSampleVariantProvider(swappingSampleVariantsProvider,
//						sampleVariantProvider.cacheSize());
//			}
            swappingSampleVariantProviders.put(sampleVariantProvider, swappingSampleVariantsProvider);
        }

        allelesUpdate.put(originalGeneticVariant, variantAlleles.getComplement());
        variantProviderUpdates.put(originalGeneticVariant, swappingSampleVariantsProvider);

        if (refAllele != null) {
            refAlleleUpdate.put(originalGeneticVariant, refAllele.getComplement());
        }

    }

    @Override
    public synchronized void updateRefAllele(ModifiableGeneticVariant geneticVariant, Allele newRefAllele) {

        GeneticVariant originalGeneticVariant = geneticVariant.getOriginalVariant();

        // If no update do nothing expect if there was already a previous
        // update. Reverting back is complicated because that would require
        // recoding the alleles. Might undo intentional changes in ordering of
        // alternative alleles
        if (originalGeneticVariant.getRefAllele() == newRefAllele
                && !refAlleleUpdate.containsKey(originalGeneticVariant)) {
            return;
        }

        Alleles variantAlleles = getUpdatedAlleles(geneticVariant);
        if (variantAlleles == null) {
            variantAlleles = originalGeneticVariant.getVariantAlleles();
        }

        if (!variantAlleles.contains(newRefAllele)) {
            throw new GenotypeDataException("Can not update to reference allele (" + newRefAllele
                    + ") is not a found in supplied alleles " + variantAlleles.getAllelesAsString()
                    + " for variant with ID: " + originalGeneticVariant.getPrimaryVariantId() + " at: "
                    + originalGeneticVariant.getSequenceName() + ":" + originalGeneticVariant.getStartPos());
        }

        // reference allele is changed so can never be the first allele so lets
        // do the reorder dance :)
        ArrayList<Allele> allelesWithoutRef = new ArrayList<Allele>(variantAlleles.getAlleles());
        allelesWithoutRef.remove(newRefAllele);
        allelesWithoutRef.add(0, newRefAllele);

        allelesUpdate.put(originalGeneticVariant, Alleles.createAlleles(allelesWithoutRef));
        refAlleleUpdate.put(originalGeneticVariant, newRefAllele);

    }

    @Override
    public synchronized Alleles getUpdatedAlleles(ModifiableGeneticVariant geneticVariant) {
        return allelesUpdate.get(geneticVariant.getOriginalVariant());
    }

    @Override
    public Iterable<ModifiableGeneticVariant> getModifiableSequenceGeneticVariants(String seqName) {
        Iterator<GeneticVariant> originalIterator = sourceGenotypeData.getSequenceGeneticVariants(seqName).iterator();
        return ModifiableGeneticVariantIterator.createModifiableGeneticVariantIterable(originalIterator, this,
                filteredOutVariants);
    }

    @Override
    public Iterable<ModifiableGeneticVariant> getModifiableVariantsByPos(String seqName, int startPos) {
        Iterator<GeneticVariant> originalIterator = sourceGenotypeData.getVariantsByPos(seqName, startPos).iterator();
        return ModifiableGeneticVariantIterator.createModifiableGeneticVariantIterable(originalIterator, this,
                filteredOutVariants);
    }

    @Override
    public ModifiableGeneticVariant getModifiableSnpVariantByPos(String seqName, int startPos) {
        GeneticVariant originalVariant = sourceGenotypeData.getSnpVariantByPos(seqName, startPos);
        if (originalVariant == null) {
            return null;
        }

        ModifiableGeneticVariant modifiableVariant = new ModifiableGeneticVariant(originalVariant, this);
        if (filteredOutVariants.contains(modifiableVariant)) {
            return null;
        } else {
            return modifiableVariant;
        }

    }

    @Override
    public Iterable<ModifiableGeneticVariant> getModifiableGeneticVariants() {
        return ModifiableGeneticVariantIterator.createModifiableGeneticVariantIterable(sourceGenotypeData.iterator(),
                this, filteredOutVariants);
    }

    @Override
    public void excludeVariant(ModifiableGeneticVariant geneticVariant) {
        filteredOutVariants.add(geneticVariant);
    }

    @Override
    public int getExcludedVariantCount() {
        return filteredOutVariants.size();
    }

    @Override
    public List<SampleAnnotation> getSampleAnnotations() {
        return sourceGenotypeData.getSampleAnnotations();
    }

    @Override
    public Annotation getSampleAnnotation(String annotationId) {
        return sourceGenotypeData.getSampleAnnotation(annotationId);
    }

    @Override
    public Iterable<GeneticVariant> getVariantsByRange(String seqName, int rangeStart, int rangeEnd) {
        return ModifiableGeneticVariantIterator.createGeneticVariantIterableBackByModifiable(sourceGenotypeData
                .getVariantsByRange(seqName, rangeStart, rangeEnd).iterator(), this, filteredOutVariants);
    }

    @Override
    public void close() throws IOException {
    }

    @Override
    public Map<String, Annotation> getVariantAnnotationsMap() {
        return sourceGenotypeData.getVariantAnnotationsMap();
    }

    @Override
    public Map<String, SampleAnnotation> getSampleAnnotationsMap() {
        return sourceGenotypeData.getSampleAnnotationsMap();
    }

    @Override
    public boolean isOnlyContaingSaveProbabilityGenotypes() {
        return sourceGenotypeData.isOnlyContaingSaveProbabilityGenotypes();
    }

    @Override
    public boolean isSwapped(GeneticVariant geneticVariant) {
        return swappedVariants.contains(geneticVariant);
    }

}
