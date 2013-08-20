package nl.umcg.deelenp.genotypealigner;

import java.util.ArrayList;
import java.util.Collection;
import java.util.Collections;
import java.util.TreeMap;

import org.apache.log4j.Logger;
import org.molgenis.genotype.RandomAccessGenotypeData;
import org.molgenis.genotype.modifiable.ModifiableGeneticVariant;
import org.molgenis.genotype.modifiable.ModifiableGenotypeData;
import org.molgenis.genotype.modifiable.ModifiableGenotypeDataInMemory;
import org.molgenis.genotype.util.Ld;
import org.molgenis.genotype.util.LdCalculator;
import org.molgenis.genotype.util.LdCalculatorException;
import org.molgenis.genotype.variant.GeneticVariant;

public class Aligner {
	
	private static Logger LOGGER = Logger.getLogger(GenotypeAligner.class);
	
	public ModifiableGenotypeData alignToRef(RandomAccessGenotypeData study, RandomAccessGenotypeData ref, double minLdToIncludeAlign, double minSnpsToAlignOn, int flankSnpsToConsider, boolean ldCheck, boolean updateId) throws LdCalculatorException{
		
		ModifiableGenotypeData aligendStudyData = new ModifiableGenotypeDataInMemory(study); 

		//The included study variants after the first loop
		ArrayList<ModifiableGeneticVariant> studyVariantList = new ArrayList<ModifiableGeneticVariant>();
		
		//The ref variants for the non-excluded study variants after the first loop in the exact same order as the study variants. 
		ArrayList<GeneticVariant> refVariantList = new ArrayList<GeneticVariant>(); 
		
		int iterationCounter = 0;
		int nonGcNonAtSnpsEncountered = 0;
		int nonGcNonAtSnpsSwapped = 0;
		
		//In this loop we filter the variants present in the reference and swap the AG, AC, TC, TG SNPs.
		studyVariants:
		for(ModifiableGeneticVariant studyVariant : aligendStudyData.getModifiableGeneticVariants()){
			
			++iterationCounter;
			
			if(iterationCounter % 1000 == 0){
				LOGGER.info("Iteration 1 - " + iterationCounter + " variants processed");
				System.out.println("Iteration 1 - " + iterationCounter + " variants processed");
			}
			
			if (! studyVariant.isMapped()){
				LOGGER.warn("Excluding variant: " + studyVariant.getPrimaryVariantId() + " Has no mapping");
				studyVariant.exclude();
				continue studyVariants;
			}
			
			if( ! studyVariant.isSnp()){
				LOGGER.warn("Excluding variant: " + studyVariant.getPrimaryVariantId() + " currently only SNPs are suppored. Feel free to contact the autors.");
				studyVariant.exclude();
				continue studyVariants;
			}
			
			if( ! studyVariant.isBiallelic()){
				LOGGER.warn("Excluding variant: " + studyVariant.getPrimaryVariantId() + " only biallelic variants currently not suppored. Feel free to contact the autors.");
				studyVariant.exclude();
				continue studyVariants;
			}
			
			if (! (studyVariant.getMinorAlleleFrequency() > 0)){
				LOGGER.warn("Excluding variant: " + studyVariant.getPrimaryVariantId() + " has a MAF of 0 in the study data");
				studyVariant.exclude();
				continue studyVariants;
			}
			
			Iterable<GeneticVariant> potentialRefVariants = ref.getVariantsByPos(studyVariant.getSequenceName(), studyVariant.getStartPos());
			
			GeneticVariant refVariant;
			
			if(!potentialRefVariants.iterator().hasNext()){
				
				LOGGER.warn("Excluding variant: " + studyVariant.getPrimaryVariantId() + " no variants at this position (" + studyVariant.getSequenceName() + ":" + studyVariant.getStartPos() + ") in the reference data");
				studyVariant.exclude();
				continue studyVariants;
				
			} else {
				
				refVariant = deterimeBestMatchRefVariant(studyVariant, potentialRefVariants);
				if(refVariant == null){
					//Not writing to log. The function already wrote exact cause why there is no best match reference variant selected
					studyVariant.exclude();
					continue studyVariants;
				}
			}
			
			if (! (refVariant.getMinorAlleleFrequency() > 0)){
				LOGGER.warn("Excluding variant: " + refVariant.getPrimaryVariantId() + " has a MAF of 0 in the reference data");
				studyVariant.exclude();
				continue studyVariants;
			}
			
			//If we get here we have found a variant is our reference data on the same position with comparable alleles.
			
			
			if(updateId && ! studyVariant.getPrimaryVariantId().equals(refVariant.getPrimaryVariantId())){
				LOGGER.debug("Updating primary variant ID of " + studyVariant.getPrimaryVariantId() + " to: " + refVariant.getPrimaryVariantId()); 
				studyVariant.updatePrimaryId(refVariant.getPrimaryVariantId());
			}
			
			
			if( ! studyVariant.isAtOrGcSnp()){
				
				++nonGcNonAtSnpsEncountered;
				//Non ambiguous SNP. We can just swap asses the strand of the reference data to determine if we need to swap
				
				if( ! studyVariant.getVariantAlleles().sameAlleles(refVariant.getVariantAlleles())){
					++nonGcNonAtSnpsSwapped;
					//Alleles do not match we need to swap our study data.
					studyVariant.swap();
					//no need to check if there is a match now. We would not have gotten here if alleles where not comparable.
					if(LOGGER.isDebugEnabled()){
						LOGGER.debug("Swapped strand of non AT and non GC SNP: " + studyVariant.getPrimaryVariantId() + " based on non ambiguous alleles. After swap study maf: " + studyVariant.getMinorAlleleFrequency() + " (" + studyVariant.getMinorAllele() + ") ref maf: " + refVariant.getMinorAlleleFrequency() + " (" + refVariant.getMinorAllele() + ")");
					}
				}
				
			}
			
			//Create list with the study variants
			studyVariantList.add(studyVariant);
			
			//Save the ref variants for the non excluded study variant 
			//This will help us in the second iteration over the study variants
			refVariantList.add(refVariant);
			
			
		}
		
		LOGGER.info("Iteration 1 - Completed, non AT and non GC SNPs are aligned " + nonGcNonAtSnpsEncountered + " found and " + nonGcNonAtSnpsSwapped + " swapped");
		System.out.println("Iteration 1 - Completed, non AT and non GC SNPs are aligned " + nonGcNonAtSnpsEncountered + " found and " + nonGcNonAtSnpsSwapped + " swapped");
		
		int removedSnpsBasedOnLdCheck = 0;
		
		//The order of the variants is not guaranteed by the genotype reader nor by all file types.
		//So order the variant here to make it easy to select 50 variants upstream or downstream
		Collections.sort(studyVariantList);
		Collections.sort(refVariantList);
		
		LOGGER.debug("Sorting of variant lists completed");
		
		//If we want to check if LD patterns matches we do this now for the AG, AC, TC, TG SNPs
		//This we give use the most optimal results for the alignment of GC and AT SNPs.
		if(ldCheck){
			
			iterationCounter = 0;
			
			//Optional second loop
			for(int variantIndex = 0 ; variantIndex < studyVariantList.size() ; ++variantIndex){
			
				++iterationCounter;
				
				if(iterationCounter % 1000 == 0){
					LOGGER.info("Iteration 2 - " + iterationCounter + " variants processed");
					System.out.println("Iteration 2 - " + iterationCounter + " variants processed");
				}
				
				ModifiableGeneticVariant studyVariant = studyVariantList.get(variantIndex);
				GeneticVariant refVariant = refVariantList.get(variantIndex);
				
				//Here we only do LD check for AG, AC, TC, TG SNPs
				if( ! studyVariant.isAtOrGcSnp()){
					
					//Correlate the haps with both these snps between study and ref
					correlationResults hapCor = correlateHaplotypes(minLdToIncludeAlign,
							flankSnpsToConsider, studyVariantList, refVariantList,
							variantIndex, studyVariant, refVariant);
					
					//Use at least min number of snps before we can draw conclusion
					if(hapCor.getTotalCor() < minSnpsToAlignOn){
						LOGGER.warn("Excluding variant: " + studyVariant.getPrimaryVariantId() + " Not enough non AT / GC in LD to check LD pattern.");
						studyVariant.exclude();
						continue;
					}
					
					if(hapCor.getPosCor() < hapCor.getNegCor()){
						//negative correlations between haplotypes more often observed. Excluding this variants
						++removedSnpsBasedOnLdCheck;
						LOGGER.warn("Excluding variant: " + studyVariant.getPrimaryVariantId() + " non AT / GC SNP with inconsistency in LD pattern \n" +
								"\tStudy: " + studyVariant.getVariantAlleles() + " maf: " + studyVariant.getMinorAlleleFrequency() + " (" + studyVariant.getMinorAllele() + ")\n" +
								"\tRef: " + refVariant.getVariantAlleles() + " maf: " + refVariant.getMinorAlleleFrequency() + " (" + refVariant.getMinorAllele() + ")\n" +
								"\tTotal snps used: " + hapCor.getTotalCor() + " pos cor: " + hapCor.getPosCor());
						studyVariant.exclude();
						continue;
					}
				}
				
			}
			
			LOGGER.info("Iteration 2 - Completed, non AT and non GC SNPs are LD checked");
			System.out.println("Iteration 2 - Completed, non AT and non GC SNPs are LD checked ");
			LOGGER.info("Excluded " + removedSnpsBasedOnLdCheck + " non AT and non GC SNPs based on inconsistencies in LD pattern");
			
			
		} else {
			System.out.println("Iteration 2 - Skipped, non AT and non GC SNPs are not LD checked ");
			LOGGER.info("Iteration 2 - Skipped, non AT and non GC SNPs are not LD checked ");
		}
		
		
		
		iterationCounter = 0;
		int GcAtSnpsEncountered = 0;
		int swapBasedOnLdCount = 0;
		removedSnpsBasedOnLdCheck = 0;
		
		//Third loop over the included variants. Now that the other variants are fixed we can focus on the GC and AT SNPs.
		for(int variantIndex = 0 ; variantIndex < studyVariantList.size() ; ++variantIndex){
			 
			++iterationCounter;
			
			if(iterationCounter % 1000 == 0){
				LOGGER.info("Iteration 3 - " + iterationCounter + " variants processed (" + GcAtSnpsEncountered +  " GC or AT SNPs checked)");
				System.out.println("Iteration 3 - " + iterationCounter + " variants processed (" + GcAtSnpsEncountered +  " GC or AT SNPs checked)");
			}
			
			ModifiableGeneticVariant studyVariant = studyVariantList.get(variantIndex);
			GeneticVariant refVariant = refVariantList.get(variantIndex);

			//Only do LD alignment on AT and GC SNPs.
			if( studyVariant.isAtOrGcSnp()){
				
				++GcAtSnpsEncountered;
				
				//Correlate the haps with both these snps between study and ref
				correlationResults hapCor = correlateHaplotypes(minLdToIncludeAlign,
						flankSnpsToConsider, studyVariantList, refVariantList,
						variantIndex, studyVariant, refVariant);
				
				//Use at least min number of snps before we can draw conclusion
				if(hapCor.getTotalCor() < minSnpsToAlignOn){
					LOGGER.warn("Excluding variant: " + studyVariant.getPrimaryVariantId() + " Not enough non AT / GC in LD to assess strand based on LD. Pos cor " + hapCor.getPosCor() + " neg cor " + hapCor.getNegCor());
					studyVariant.exclude();
					continue;
				}
				
				
				if(hapCor.getPosCor() < hapCor.getNegCor()){
					//negative correlation more often observed. We need to swap the strand of this SNP.
					studyVariant.swap();
					++swapBasedOnLdCount;
					
					if(LOGGER.isDebugEnabled()){
						LOGGER.debug("Swapped strand of AT or GC SNP: " + studyVariant.getPrimaryVariantId() + " based on LD. After swap study maf: " + studyVariant.getMinorAlleleFrequency() + " (" + studyVariant.getMinorAllele() + ") ref maf: " + refVariant.getMinorAlleleFrequency() + " (" + refVariant.getMinorAllele() + ")");
					}
					
					if(ldCheck){
						//Ld pattern should be okay now. but we are going to do the extra check
						
						//Correlate the haps with both these snps between study and ref
						correlationResults hapCorSwapped = correlateHaplotypes(minLdToIncludeAlign,
								flankSnpsToConsider, studyVariantList, refVariantList,
								variantIndex, studyVariant, refVariant);
						
						//No need to check the count. Already done when checking unswapped LD pattern.					
						
						if(hapCorSwapped.getPosCor() < hapCorSwapped.getNegCor()){
							//LD pattern still not intact. Excluding variant
							++removedSnpsBasedOnLdCheck;
							LOGGER.warn("Excluding variant: " + studyVariant.getPrimaryVariantId() + " GC or AT SNP with inconsistency in LD pattern that is not solved by swapping");
							studyVariant.exclude();
							continue;
							}
						
					}
					
				}
				
				//No need for LD check here. If it would not have matched it would have gotten in the swapping part.
				//LD is checked again after swapping if requested.
				
				
			}
			
		}
		
		if(ldCheck){
			LOGGER.info("Iteration 3 - Completed, non AT and non GC SNPs are aligned and LD check afterwards");
			System.out.println("Iteration 3 - Completed, non AT and non GC SNPs are aligned and LD check afterwards");
		}
		else {
			LOGGER.info("Iteration 3 - Completed, non AT and non GC SNPs are aligned. Extra LD check skipped");
			System.out.println("Iteration 3 - Completed, non AT and non GC SNPs are aligned. Extra LD check skipped");
		}
			
		if(ldCheck){
			LOGGER.info("Excluded " + removedSnpsBasedOnLdCheck + " AT or GC variants based on LD patterns");
			System.out.println("Excluded " + removedSnpsBasedOnLdCheck + " AT or GC variants based on LD patterns");
		}
		LOGGER.info("Swapped " + swapBasedOnLdCount + " out " + GcAtSnpsEncountered + " AT or GC variants based on LD patterns");
		System.out.println("Swapped " + swapBasedOnLdCount + " AT or GC variants based on LD patterns");
		
		return aligendStudyData;

	}

	private correlationResults correlateHaplotypes(double minLdToIncludeAlignBase,
			int flankSnpsToConsider,
			ArrayList<ModifiableGeneticVariant> studyVariantList,
			ArrayList<GeneticVariant> refVariantList, int variantIndex,
			GeneticVariant snpStudyVariant, GeneticVariant refVariant)
			{
		
		int posCor = 0;
		int negCor = 0;
		
		//loop over potential variants variants to determine strand
		otherVariantsLoop:
		for(int otherVariantIndex = Math.max(0, variantIndex - flankSnpsToConsider) ; otherVariantIndex < variantIndex + flankSnpsToConsider && otherVariantIndex < studyVariantList.size() ; ++ otherVariantIndex){
			
			//Do not use self
			if(variantIndex == otherVariantIndex){
				continue otherVariantsLoop;
			}
			
			GeneticVariant otherSnpStudyVariant = studyVariantList.get(otherVariantIndex);
			
			//Only use variants on same chromosome
			if(!snpStudyVariant.getSequenceName().equals(otherSnpStudyVariant.getSequenceName())){
				continue otherVariantsLoop;
			}
			
			//Do not use other AT or GC for LD alignment
			if(otherSnpStudyVariant.isAtOrGcSnp()){
				continue otherVariantsLoop;
			}
			
			GeneticVariant otherRefVariant = refVariantList.get(otherVariantIndex);
			
			Ld ldStudy;
			Ld ldRef;
			try
			{
				ldStudy = LdCalculator.calculateLd(snpStudyVariant, otherSnpStudyVariant);
				ldRef = LdCalculator.calculateLd(refVariant, otherRefVariant);
			}
			catch (LdCalculatorException e)
			{
				LOGGER.warn("Error in LD calculation, skipping this comparsion when comparing haplotype structure. Following error occured: " + e.getMessage());
				continue;
			}
			
			//only use SNPs with min R2 in both study as ref
			if(ldStudy.getR2() >= minLdToIncludeAlignBase && ldRef.getR2() >= minLdToIncludeAlignBase){
				
				//Put in tree map to sort haplotypes. This can differ in the case of different reference allele
				TreeMap<String, Double> studyHapFreq = new TreeMap<String, Double>(ldStudy.getHaplotypesFreq());
				TreeMap<String, Double> refHapFreq = new TreeMap<String, Double>(ldRef.getHaplotypesFreq());
				
				double[] studyHapFreqArray = createDoubleArrayFromCollection(studyHapFreq.values());
				double[] refHapFreqArray = createDoubleArrayFromCollection(refHapFreq.values());
				
				//Correlate study haplotypes to ref haplotypes.
				double correlation = JSci.maths.ArrayMath.correlation(studyHapFreqArray, refHapFreqArray);
				
				if(correlation < 0){
					++negCor;
				} else {
					++posCor;
				}
				
			}
			
		}
		
		return new correlationResults(posCor, negCor);
	}

	private double[] createDoubleArrayFromCollection(
			Collection<Double> values) {
		

		double[] array = new double[values.size()];
		
		int i = 0;
		for(Double d : values){
			array[i] = d;
			++i;
		}
		
		return array;
	}

	private GeneticVariant deterimeBestMatchRefVariant(
			GeneticVariant studyVariant, Iterable<GeneticVariant> refVariants) {
		
		for(GeneticVariant refVariant : refVariants){
			
			if(refVariant.getVariantId().isSameId(studyVariant.getVariantId())){
				
				//test if same alleles or complement
				if(refVariant.getVariantAlleles().sameAlleles(studyVariant.getVariantAlleles()) 
						|| refVariant.getVariantAlleles().sameAlleles(studyVariant.getVariantAlleles().getComplement())){
					
					return refVariant;
					
				} else {
					
					LOGGER.warn("Excluding variant: " + studyVariant.getPrimaryVariantId() + " Found variant with same ID but alleles are not comparable.");
					return null;
				}
			}
		}
		
		GeneticVariant potentialRef = null;
		for(GeneticVariant refVariant : refVariants){
			
			//test if same alleles or complement
			if(refVariant.getVariantAlleles().sameAlleles(studyVariant.getVariantAlleles()) 
					|| refVariant.getVariantAlleles().sameAlleles(studyVariant.getVariantAlleles().getComplement())){
				
				if(potentialRef == null){
					potentialRef = refVariant;
				} else {
					LOGGER.warn("Excluding variant: " + studyVariant.getPrimaryVariantId() + " because position maps to multiple variants with same alleles. Neither of these variants have same ID as this variant. No way to know which the correspoding variant should be.");
					return null;
				}
				
			}
		}
		if(potentialRef == null){
			LOGGER.warn("Excluding variant: " + studyVariant.getPrimaryVariantId() + " There is no variant in the reference at this position with same ID or same alleles");
			return null;
		}
		return potentialRef;
		
	}
	
	private static class correlationResults{
		
		private final int posCor;
		private final int negCor;
		
		public correlationResults(int posCor, int negCor) {
			super();
			this.posCor = posCor;
			this.negCor = negCor;
		}
		
		public int getPosCor() {
			return posCor;
		}
		
		public int getNegCor() {
			return negCor;
		}
		
		public int getTotalCor() {
			return getPosCor() + getNegCor();
		}
		
	}

}
