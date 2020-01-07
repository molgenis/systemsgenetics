package eqtlmappingpipeline;

import java.io.File;
import java.util.List;
import org.apache.commons.math3.stat.descriptive.SummaryStatistics;
import org.molgenis.genotype.Allele;
import org.molgenis.genotype.Alleles;
import org.molgenis.genotype.GenotypedDataWriterFormats;
import org.molgenis.genotype.RandomAccessGenotypeData;
import org.molgenis.genotype.editable.EditableGenotypeDosageSampleVariantsProvider;
import org.molgenis.genotype.editable.GenotypeDataCustomVariantProvider;
import org.molgenis.genotype.multipart.MultiPartGenotypeData;
import org.molgenis.genotype.util.CalledDosageConvertor;
import org.molgenis.genotype.variant.GeneticVariant;
import org.molgenis.genotype.vcf.VcfGenotypeData;

/**
 *
 * @author Patrick Deelen
 */
public class CreateGeuvadisTestTrityper {

    public void createGeuvadisTestTrityper() throws Exception {
		
		//GeuvadisTestData/Geuvadis_chr1.vcf.gz
		File chr1File = new File((this.getClass().getResource("/GeuvadisTestData/Geuvadis_chr1.vcf.gz")).toURI());
		File chr2File = new File((this.getClass().getResource("/GeuvadisTestData/Geuvadis_chr2.vcf.gz")).toURI());
		
		RandomAccessGenotypeData genotypeData = new MultiPartGenotypeData(new VcfGenotypeData(chr1File, 100, 0.8), new VcfGenotypeData(chr2File, 100, 0.8));
		
		int sampleCount = genotypeData.getSamples().size();
		
		GenotypeDataCustomVariantProvider<EditableGenotypeDosageSampleVariantsProvider> modifiedGenotypeData = new GenotypeDataCustomVariantProvider<EditableGenotypeDosageSampleVariantsProvider>(genotypeData, genotypeData.getSamples(), new EditableGenotypeDosageSampleVariantsProvider(sampleCount));
		
		for(GeneticVariant variant : genotypeData){
			
			List<Alleles> alleles = variant.getSampleVariants();
			float[] dosages = variant.getSampleDosages();
			
			float[] convertedDosages = CalledDosageConvertor.convertCalledAllelesToDosage(alleles, variant.getVariantAlleles(), variant.getVariantAlleles().get(0));
			
			SummaryStatistics sumStat = new SummaryStatistics();
			for(int s = 0 ; s < sampleCount ; ++s){
				if(dosages[s] >= 0){
					sumStat.addValue(dosages[s]);
				}
			}
			float meanDosage = (float) sumStat.getMean();
			
			for(int s = 0 ; s < sampleCount ; ++s){
				if(dosages[s] < 0){
					if(alleles.get(s).contains(Allele.ZERO)){
						dosages[s] = meanDosage;
					} else {
						dosages[s] = convertedDosages[s];
					}
				} else {
					//this is a safty check
					if(dosages[s] < 1f && convertedDosages[s] > 1f || dosages[s] > 1f && convertedDosages[s] < 1f){
						throw new Exception();
					}
				}
			}
			
			
			modifiedGenotypeData.getSampleVariantProvider().setDosageAndGenotypes(variant, alleles, dosages);
			
		}
		
		
		File folder = new File((this.getClass().getResource("/GeuvadisTestData/")).toURI());
		File outputFile = new File(folder, "trityper");
		
		System.out.println(outputFile.getAbsolutePath());
		
		GenotypedDataWriterFormats.TRITYPER.createGenotypeWriter(modifiedGenotypeData, sampleCount).write(outputFile.getAbsolutePath());
		
		
    }
	
}
