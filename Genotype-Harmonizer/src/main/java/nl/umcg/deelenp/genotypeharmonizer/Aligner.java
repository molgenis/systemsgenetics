package nl.umcg.deelenp.genotypeharmonizer;

import static JSci.maths.ArrayMath.covariance;
import com.google.common.collect.Lists;
import java.io.BufferedWriter;
import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Collection;
import java.util.Collections;
import java.util.Iterator;
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

/**
 *
 * @author Patrick Deelen
 */
public class Aligner {

	private static final Logger LOGGER = Logger.getLogger(GenotypeHarmonizer.class);

	/**
	 *
	 * @param study data to align
	 * @param ref reference for alignment
	 * @param minLdToIncludeAlign
	 * @param minSnpsToAlignOn
	 * @param flankSnpsToConsider
	 * @param ldCheck
	 * @param updateId
	 * @param keep
	 * @param snpUpdateFile file to write ID updates to.
	 * @param maxMafForMafAlignment
	 * @param snpLogFile
	 * @param matchRefAllele
	 * @return
	 * @throws LdCalculatorException
	 * @throws java.io.IOException
	 * @throws nl.umcg.deelenp.genotypeharmonizer.GenotypeAlignmentException
	 */
	public ModifiableGenotypeData alignToRef(RandomAccessGenotypeData study, RandomAccessGenotypeData ref, double minLdToIncludeAlign, double minSnpsToAlignOn, int flankSnpsToConsider, boolean ldCheck, final boolean updateId, boolean keep, File snpUpdateFile, double maxMafForMafAlignment, File snpLogFile, boolean matchRefAllele) throws LdCalculatorException, IOException, GenotypeAlignmentException {

		ModifiableGenotypeData aligendStudyData = new ModifiableGenotypeDataInMemory(study);

		//The included study variants after the first loop
		ArrayList<ModifiableGeneticVariant> studyVariantList = new ArrayList<>();

		//The ref variants for the non-excluded study variants after the first loop in the exact same order as the study variants. 
		ArrayList<GeneticVariant> refVariantList = new ArrayList<>();

		BufferedWriter snpUpdateWriter = null;
		if (updateId) {
			snpUpdateWriter = new BufferedWriter(new FileWriter(snpUpdateFile));
			snpUpdateWriter.append("chr\tpos\toriginalId\tnewId\n");
		}

		SnpLogWriter snpLogWriter = new SnpLogWriter(snpLogFile);

		int iterationCounter = 0;
		int nonGcNonAtSnpsEncountered = 0;
		int nonGcNonAtSnpsSwapped = 0;

		//In this loop we filter the variants present in the reference and swap the AG, AC, TC, TG SNPs.
		studyVariants:
		for (ModifiableGeneticVariant studyVariant : aligendStudyData.getModifiableGeneticVariants()) {

			++iterationCounter;

			if (iterationCounter % 10000 == 0) {
				//LOGGER.info("Iteration 1 - " + GenotypeHarmonizer.DEFAULT_NUMBER_FORMATTER.format(iterationCounter) + " variants processed");
				System.out.println("Iteration 1 - " + GenotypeHarmonizer.DEFAULT_NUMBER_FORMATTER.format(iterationCounter) + " variants processed");
			}

			if (!studyVariant.isMapped()) {
				snpLogWriter.addToLog(studyVariant, SnpLogWriter.Actions.EXCLUDED, "No mapping");
				studyVariant.exclude();
				continue studyVariants;
			}

			if (studyVariant.getStartPos() == 0) {
				snpLogWriter.addToLog(studyVariant, SnpLogWriter.Actions.EXCLUDED, "No mapping");
				studyVariant.exclude();
				continue studyVariants;
			}

//			if (!studyVariant.isSnp()) {
//				snpLogWriter.addToLog(studyVariant, SnpLogWriter.Actions.EXCLUDED, "Not a SNP");
//				studyVariant.exclude();
//				continue studyVariants;
//			}
//
//			if (!studyVariant.isBiallelic()) {
//				snpLogWriter.addToLog(studyVariant, SnpLogWriter.Actions.EXCLUDED, "Not biallelic");
//				studyVariant.exclude();
//				continue studyVariants;
//			}

			Iterator<GeneticVariant> potentialRefVariants = ref.getVariantsByPos(studyVariant.getSequenceName(), studyVariant.getStartPos()).iterator();

			GeneticVariant refVariant = null;

			if (!potentialRefVariants.hasNext()) {

				if (keep) {
					//LOGGER.warn("No ref variant found for: " + studyVariant.getPrimaryVariantId() + " variant will not be aligned but will be written to output because of --keep");
				} else {
					snpLogWriter.addToLog(studyVariant, SnpLogWriter.Actions.EXCLUDED, "No variants at this position in reference");
					studyVariant.exclude();
				}

				continue studyVariants;

			} else {

				ArrayList<GeneticVariant> potentialRefVariantsList = Lists.newArrayList(potentialRefVariants);

				//Find ref based on ID
				for (GeneticVariant potentialRefVariant : potentialRefVariantsList) {

					if (potentialRefVariant.getVariantId().isSameId(studyVariant.getVariantId())) {

						//TODO only take complement for SNP variants
						//TODO how to handle if study allelles are subset of ref alleles
						
						//test if same alleles or complement
						if (potentialRefVariant.getVariantAlleles().sameAlleles(studyVariant.getVariantAlleles())
								|| potentialRefVariant.getVariantAlleles().sameAlleles(studyVariant.getVariantAlleles().getComplement())) {

							refVariant = potentialRefVariant;

						} else {
							snpLogWriter.addToLog(studyVariant, SnpLogWriter.Actions.EXCLUDED, "Found variant with same ID but alleles are not comparable");
							studyVariant.exclude();
							continue studyVariants;
						}
					}
				}

				if (refVariant == null) {

					//Find ref based on Alleles
					for (GeneticVariant potentialRefVariant : potentialRefVariantsList) {

						//TODO only take complement for SNP variants
						//TODO how to handle if study allelles are subset of ref alleles. Not possible for AT or GC but is possible for other variants
						
						//test if same alleles or complement
						if (potentialRefVariant.getVariantAlleles().sameAlleles(studyVariant.getVariantAlleles())
								|| potentialRefVariant.getVariantAlleles().sameAlleles(studyVariant.getVariantAlleles().getComplement())) {

							if (refVariant == null) {
								refVariant = potentialRefVariant;
							} else {
								snpLogWriter.addToLog(studyVariant, SnpLogWriter.Actions.EXCLUDED, "Position maps to multiple variants with same alleles. Neither of these variants have same ID as this variant. No way to know what the corresponding variant is");
								studyVariant.exclude();
								continue studyVariants;
							}

						}
					}

					if (refVariant == null) {
						if (keep) {
							//LOGGER.warn("No ref variant found for: " + studyVariant.getPrimaryVariantId() + " variant will not be aligned but will be written to output because of --keep");
						} else {
							snpLogWriter.addToLog(studyVariant, SnpLogWriter.Actions.EXCLUDED, "No variant in the reference at this position with same ID or same alleles");
							studyVariant.exclude();
						}
						continue studyVariants;
					}

				}

			}

			if (
                    updateId && 
                    !(refVariant.getPrimaryVariantId() == null && studyVariant.getPrimaryVariantId() == null) && 
                    (refVariant.getPrimaryVariantId() == null || studyVariant.getPrimaryVariantId() == null || !studyVariant.getPrimaryVariantId().equals(refVariant.getPrimaryVariantId()))) {
				snpUpdateWriter.append(studyVariant.getSequenceName());
				snpUpdateWriter.append('\t');
				snpUpdateWriter.append(String.valueOf(studyVariant.getStartPos()));
				snpUpdateWriter.append('\t');
				snpUpdateWriter.append(studyVariant.getPrimaryVariantId());
				snpUpdateWriter.append('\t');
				snpUpdateWriter.append(refVariant.getPrimaryVariantId());
				snpUpdateWriter.append('\n');

				LOGGER.debug("Updating primary variant ID of " + studyVariant.getPrimaryVariantId() + " to: " + refVariant.getPrimaryVariantId());
				studyVariant.updatePrimaryId(refVariant.getPrimaryVariantId());
			}


			if (!studyVariant.isAtOrGcSnp()) {
				
				//TODO based on previous contraits this should never happen for indels but lets make sure

				++nonGcNonAtSnpsEncountered;
				//Non ambiguous SNP. We can just swap assess the strand of the reference data to determine if we need to swap

				if (!studyVariant.getVariantAlleles().sameAlleles(refVariant.getVariantAlleles())) {
					++nonGcNonAtSnpsSwapped;
					//Alleles do not match we need to swap our study data.
					studyVariant.swap();
					//no need to check if there is a match now. We would not have gotten here if alleles where not comparable.
					snpLogWriter.addToLog(studyVariant, SnpLogWriter.Actions.SWAPPED, "");
//					if (LOGGER.isDebugEnabled()) {	
//						LOGGER.debug("Swapped strand of non AT and non GC SNP: " + studyVariant.getPrimaryVariantId() + " based on non ambiguous alleles. After swap study maf: " + studyVariant.getMinorAlleleFrequency() + " (" + studyVariant.getMinorAllele() + ") ref maf: " + refVariant.getMinorAlleleFrequency() + " (" + refVariant.getMinorAllele() + ")");
//					}
				}
				if(matchRefAllele){
                    studyVariant.updateRefAllele(refVariant.getRefAllele());
                }
			}

			//Create list with the study variants
			studyVariantList.add(studyVariant);

			//Save the ref variants for the non excluded study variant 
			//This will help us in the second iteration over the study variants
			refVariantList.add(refVariant);


		}

		if (updateId) {
			snpUpdateWriter.close();
		}

		if (iterationCounter == 0) {
			throw new GenotypeAlignmentException("No variants where found in the input genotype data. Please check your variant filter options");
		}

		LOGGER.info("Iteration 1 - Completed, non A/T and non G/C variants are aligned " + GenotypeHarmonizer.DEFAULT_NUMBER_FORMATTER.format(nonGcNonAtSnpsEncountered) + " found and " + GenotypeHarmonizer.DEFAULT_NUMBER_FORMATTER.format(nonGcNonAtSnpsSwapped) + " swapped");
		System.out.println("Iteration 1 - Completed, non A/T and non G/C variants are aligned " + GenotypeHarmonizer.DEFAULT_NUMBER_FORMATTER.format(nonGcNonAtSnpsEncountered) + " found and " + GenotypeHarmonizer.DEFAULT_NUMBER_FORMATTER.format(nonGcNonAtSnpsSwapped) + " swapped");

		if (studyVariantList.isEmpty()) {
			snpLogWriter.close();
			throw new GenotypeAlignmentException("Zero of the input variants found in reference set. Are both datasets the same genome build? Perhapse you need use --forceChr.");
		}

		int removedSnpsBasedOnLdCheck = 0;

		//The order of the variants is not guaranteed by the genotype reader nor by all file types.
		//So order the variant here to make it easy to select specified number of variants upstream or downstream
		Collections.sort(studyVariantList);
		Collections.sort(refVariantList);

		LOGGER.debug("Sorting of variant lists completed");

		//If we want to check if LD patterns matches we do this now for the AG, AC, TC, TG SNPs
		//This we give use the most optimal results for the alignment of GC and AT SNPs.
		if (ldCheck) {

			iterationCounter = 0;

			//Optional second loop
			for (int variantIndex = 0; variantIndex < studyVariantList.size(); ++variantIndex) {

				++iterationCounter;

				if (iterationCounter % 10000 == 0) {
					//LOGGER.info("Iteration 2 - " + GenotypeHarmonizer.DEFAULT_NUMBER_FORMATTER.format(iterationCounter) + " variants processed");
					System.out.println("Iteration 2 - " + GenotypeHarmonizer.DEFAULT_NUMBER_FORMATTER.format(iterationCounter) + " variants processed");
				}

				ModifiableGeneticVariant studyVariant = studyVariantList.get(variantIndex);
				GeneticVariant refVariant = refVariantList.get(variantIndex);

				//Here we only do LD check for AG, AC, TC, TG SNPs
				if (!studyVariant.isAtOrGcSnp()) {

					//TODO only do this if both the study and ref variant are bialilic. This need to be mentioned in the manual
					
					//Correlate the haps with both these snps between study and ref
					CorrelationResults hapCor = correlateHaplotypes(minLdToIncludeAlign,
							flankSnpsToConsider, studyVariantList, refVariantList,
							variantIndex, studyVariant, refVariant);

					//Use at least min number of snps before we can draw conclusion
					if (hapCor.getTotalCor() < minSnpsToAlignOn) {
						snpLogWriter.addToLog(studyVariant, SnpLogWriter.Actions.EXCLUDED, "Not enough non A/T or G/C in LD to check LD pattern");
						studyVariant.exclude();
						continue;
					}

					if (hapCor.getPosCor() < hapCor.getNegCor()) {
						//negative correlations between haplotypes more often observed. Excluding this variants
						++removedSnpsBasedOnLdCheck;
						snpLogWriter.addToLog(studyVariant, SnpLogWriter.Actions.EXCLUDED, "Non A/T or G/C SNP with inconsistency in LD pattern");
//						LOGGER.warn("Excluding variant: " + studyVariant.getPrimaryVariantId() + " non AT / GC SNP with inconsistency in LD pattern \n"
//								+ "\tStudy: " + studyVariant.getVariantAlleles() + " maf: " + studyVariant.getMinorAlleleFrequency() + " (" + studyVariant.getMinorAllele() + ")\n"
//								+ "\tRef: " + refVariant.getVariantAlleles() + " maf: " + refVariant.getMinorAlleleFrequency() + " (" + refVariant.getMinorAllele() + ")\n"
//								+ "\tTotal variants used: " + hapCor.getTotalCor() + " pos cor: " + hapCor.getPosCor());
						studyVariant.exclude();
						continue;
					}
				}
				
			}

			LOGGER.info("Iteration 2 - Completed, non A/T and non G/C SNPs are LD checked");
			System.out.println("Iteration 2 - Completed, non A/T and non G/C SNPs are LD checked ");
			LOGGER.info("Excluded " + GenotypeHarmonizer.DEFAULT_NUMBER_FORMATTER.format(removedSnpsBasedOnLdCheck) + " non A/T and non G/C SNPs based on inconsistencies in LD pattern");


		} else {
			System.out.println("Iteration 2 - Skipped, non A/T and non G/C SNPs are not LD checked ");
			LOGGER.info("Iteration 2 - Skipped, non A/T and non G/C SNPs are not LD checked ");
		}



		iterationCounter = 0;
		int GcAtSnpsEncountered = 0;
		int swapBasedOnLdCount = 0;
		removedSnpsBasedOnLdCheck = 0;

		//Third loop over the included variants. Now that the other variants are fixed we can focus on the GC and AT SNPs.
		for (int variantIndex = 0; variantIndex < studyVariantList.size(); ++variantIndex) {

			++iterationCounter;

			if (iterationCounter % 10000 == 0) {
				//LOGGER.info("Iteration 3 - " + GenotypeHarmonizer.DEFAULT_NUMBER_FORMATTER.format(iterationCounter) + " variants processed (" + GenotypeHarmonizer.DEFAULT_NUMBER_FORMATTER.format(GcAtSnpsEncountered) + " GC or AT SNPs checked)");
				System.out.println("Iteration 3 - " + GenotypeHarmonizer.DEFAULT_NUMBER_FORMATTER.format(iterationCounter) + " variants processed (" + GenotypeHarmonizer.DEFAULT_NUMBER_FORMATTER.format(GcAtSnpsEncountered) + " G/C or A/T SNPs checked)");
			}

			ModifiableGeneticVariant studyVariant = studyVariantList.get(variantIndex);
			GeneticVariant refVariant = refVariantList.get(variantIndex);

			//Only do LD alignment on AT and GC SNPs.
			if (studyVariant.isAtOrGcSnp()) {

				//TODO this is why we can't work with study GC or AT snps that are subset of alleles in reference. We can't do LD calculations
				
				++GcAtSnpsEncountered;

				//Correlate the haps with both these snps between study and ref
				CorrelationResults hapCor = correlateHaplotypes(minLdToIncludeAlign,
						flankSnpsToConsider, studyVariantList, refVariantList,
						variantIndex, studyVariant, refVariant);

				//Use at least min number of snps before we can draw conclusion, maybe use MA as backup
				if ((hapCor.getTotalCor() < minSnpsToAlignOn || hapCor.getPosCor() == hapCor.getNegCor())
						&& !ldCheck
						&& studyVariant.getMinorAlleleFrequency() <= maxMafForMafAlignment
						&& refVariant.getMinorAlleleFrequency() <= maxMafForMafAlignment) {


					//LOGGER.warn("Using minor allele to determine strand of: " + studyVariant.getPrimaryVariantId() + " study MAF: " + studyVariant.getMinorAlleleFrequency() + "(" + studyVariant.getMinorAllele() + ")" + " reference MAF: " + refVariant.getMinorAlleleFrequency() + "(" + refVariant.getMinorAllele() + ")");

					if (studyVariant.getMinorAllele() != refVariant.getMinorAllele()) {
						studyVariant.swap();
						++swapBasedOnLdCount;
						snpLogWriter.addToLog(studyVariant, SnpLogWriter.Actions.SWAPPED, "Based on minor allele, study MAF: " + studyVariant.getMinorAlleleFrequency() + "(" + studyVariant.getMinorAllele() + ")" + " reference MAF: " + refVariant.getMinorAlleleFrequency() + "(" + refVariant.getMinorAllele() + ")");
						//LOGGER.debug("Swapped " + studyVariant.getPrimaryVariantId() + " using the minor allele");
					} else if (LOGGER.isDebugEnabled()) {
						snpLogWriter.addToLog(studyVariant, SnpLogWriter.Actions.MAINTAINED, "Based on minor allele, study MAF: " + studyVariant.getMinorAlleleFrequency() + "(" + studyVariant.getMinorAllele() + ")" + " reference MAF: " + refVariant.getMinorAlleleFrequency() + "(" + refVariant.getMinorAllele() + ")");
					}

				} else if (hapCor.getTotalCor() < minSnpsToAlignOn) {

					snpLogWriter.addToLog(studyVariant, SnpLogWriter.Actions.EXCLUDED, "Not enough non A/T or non G/C in LD to assess strand based on LD. Pos cor " + hapCor.getPosCor() + " neg cor " + hapCor.getNegCor() + " MAF study: " + studyVariant.getMinorAlleleFrequency() + "(" + studyVariant.getMinorAllele() + ")" + " MAF reference: " + refVariant.getMinorAlleleFrequency() + "(" + refVariant.getMinorAllele() + ")");
					studyVariant.exclude();
					continue;

				} else if (hapCor.getPosCor() == hapCor.getNegCor()) {

					snpLogWriter.addToLog(studyVariant, SnpLogWriter.Actions.EXCLUDED, "Equal number of positive and negative correlations. Pos cor " + hapCor.getPosCor() + " neg cor " + hapCor.getNegCor() + " MAF study: " + studyVariant.getMinorAlleleFrequency() + "(" + studyVariant.getMinorAllele() + ")" + " MAF reference: " + refVariant.getMinorAlleleFrequency() + "(" + refVariant.getMinorAllele() + ")");
					studyVariant.exclude();
					continue;

				} else if (hapCor.getPosCor() < hapCor.getNegCor()) {
					//negative correlation more often observed. We need to swap the strand of this SNP.
					studyVariant.swap();
					++swapBasedOnLdCount;



					if (LOGGER.isDebugEnabled()) {
						snpLogWriter.addToLog(studyVariant, SnpLogWriter.Actions.SWAPPED, "Based on LD. Pos cor " + hapCor.getPosCor() + " neg cor " + hapCor.getNegCor() + " MAF study: " + studyVariant.getMinorAlleleFrequency() + "(" + studyVariant.getMinorAllele() + ")" + " MAF reference: " + refVariant.getMinorAlleleFrequency() + "(" + refVariant.getMinorAllele() + ")");
					} else {
						snpLogWriter.addToLog(studyVariant, SnpLogWriter.Actions.SWAPPED, "Based on LD");
					}

					if (ldCheck) {
						//Ld pattern should be okay now. but we are going to do the extra check

						//Correlate the haps with both these snps between study and ref
						CorrelationResults hapCorSwapped = correlateHaplotypes(minLdToIncludeAlign,
								flankSnpsToConsider, studyVariantList, refVariantList,
								variantIndex, studyVariant, refVariant);

						//No need to check the count. Already done when checking unswapped LD pattern.

						if (hapCorSwapped.getPosCor() < hapCorSwapped.getNegCor()) {
							//LD pattern still not intact. Excluding variant
							++removedSnpsBasedOnLdCheck;
							snpLogWriter.addToLog(studyVariant, SnpLogWriter.Actions.EXCLUDED, "G/C or A/T SNP with inconsistency in LD pattern that is not solved by swapping");
							studyVariant.exclude();
							continue;
						}

					}

				} else if (LOGGER.isDebugEnabled()) {
					snpLogWriter.addToLog(studyVariant, SnpLogWriter.Actions.MAINTAINED, "Based on LD. Pos cor " + hapCor.getPosCor() + " neg cor " + hapCor.getNegCor() + " MAF study: " + studyVariant.getMinorAlleleFrequency() + "(" + studyVariant.getMinorAllele() + ")" + " MAF reference: " + refVariant.getMinorAlleleFrequency() + "(" + refVariant.getMinorAllele() + ")");
				}
				if(matchRefAllele){
                    studyVariant.updateRefAllele(refVariant.getRefAllele());
                }
				//No need for LD check here. If it would not have matched it would have gotten in the swapping part.
				//LD is checked again after swapping if requested.


			}

		}

		if (ldCheck) {
			LOGGER.info("Iteration 3 - Completed, non A/T and non G/C SNPs are aligned and LD check afterwards");
			System.out.println("Iteration 3 - Completed, non A/T and non G/C SNPs are aligned and LD check afterwards");
		} else {
			LOGGER.info("Iteration 3 - Completed, non A/T and non G/C SNPs are aligned. Extra LD check skipped");
			System.out.println("Iteration 3 - Completed, non A/T and non G/C SNPs are aligned. Extra LD check skipped");
		}

		if (ldCheck) {
			LOGGER.info("Excluded " + GenotypeHarmonizer.DEFAULT_NUMBER_FORMATTER.format(removedSnpsBasedOnLdCheck) + " A/T or G/C variants based on LD patterns");
			System.out.println("Excluded " + GenotypeHarmonizer.DEFAULT_NUMBER_FORMATTER.format(removedSnpsBasedOnLdCheck) + " A/T or G/C variants based on LD patterns");
		}
		LOGGER.info("Swapped " + GenotypeHarmonizer.DEFAULT_NUMBER_FORMATTER.format(swapBasedOnLdCount) + " out of " + GenotypeHarmonizer.DEFAULT_NUMBER_FORMATTER.format(GcAtSnpsEncountered) + " A/T or G/C variants based on LD patterns");
		System.out.println("Swapped " + GenotypeHarmonizer.DEFAULT_NUMBER_FORMATTER.format(swapBasedOnLdCount) + " A/T or G/C variants based on LD patterns");

		snpLogWriter.close();

		return aligendStudyData;

	}

	private CorrelationResults correlateHaplotypes(double minLdToIncludeAlignBase,
			int flankSnpsToConsider,
			ArrayList<ModifiableGeneticVariant> studyVariantList,
			ArrayList<GeneticVariant> refVariantList, int variantIndex,
			GeneticVariant snpStudyVariant, GeneticVariant refVariant) {

		int posCor = 0;
		int negCor = 0;

//		if(snpStudyVariant.getPrimaryVariantId().equals("rs1001945")){
//		LOGGER.debug("Alignment of: " + snpStudyVariant.getPrimaryVariantId() + 
//				"\nstudy alleles: " + snpStudyVariant.getVariantAlleles() + " ref alleles: " + refVariant.getVariantAlleles() + "\n"
//				+ "maf study: " + snpStudyVariant.getMinorAlleleFrequency() + "(" + snpStudyVariant.getMinorAllele() + ") maf ref: " + refVariant.getMinorAlleleFrequency() + "(" + refVariant.getMinorAllele() + ")");
//		}
		//loop over potential variants variants to determine strand
		otherVariantsLoop:
		for (int otherVariantIndex = Math.max(0, variantIndex - flankSnpsToConsider); otherVariantIndex < variantIndex + flankSnpsToConsider && otherVariantIndex < studyVariantList.size(); ++otherVariantIndex) {

			//Do not use self
			if (variantIndex == otherVariantIndex) {
				continue otherVariantsLoop;
			}
			
			//TODO skip over multi allelelic

			GeneticVariant otherSnpStudyVariant = studyVariantList.get(otherVariantIndex);

			//Only use variants on same chromosome
			if (!snpStudyVariant.getSequenceName().equals(otherSnpStudyVariant.getSequenceName())) {
				continue otherVariantsLoop;
			}

			//Do not use other AT or GC for LD alignment
			if (otherSnpStudyVariant.isAtOrGcSnp()) {
				continue otherVariantsLoop;
			}

			GeneticVariant otherRefVariant = refVariantList.get(otherVariantIndex);

			Ld ldStudy;
			Ld ldRef;
			try {
				ldStudy = LdCalculator.calculateLd(snpStudyVariant, otherSnpStudyVariant);
				ldRef = LdCalculator.calculateLd(refVariant, otherRefVariant);
			} catch (LdCalculatorException e) {
				LOGGER.debug("Error in LD calculation, skipping this comparison when comparing haplotype structure. Following error occurred: " + e.getMessage());
				continue;
			}

			//only use SNPs with min R2 in both study as ref
			if ( !Double.isNaN(ldStudy.getR2()) && !Double.isNaN(ldRef.getR2()) && ldStudy.getR2() >= minLdToIncludeAlignBase && ldRef.getR2() >= minLdToIncludeAlignBase) {

				//Put in tree map to sort haplotypes. This can differ in the case of different reference allele
				TreeMap<String, Double> studyHapFreq = new TreeMap<String, Double>(ldStudy.getHaplotypesFreq());
				TreeMap<String, Double> refHapFreq = new TreeMap<String, Double>(ldRef.getHaplotypesFreq());

				double[] studyHapFreqArray = createDoubleArrayFromCollection(studyHapFreq.values());
				double[] refHapFreqArray = createDoubleArrayFromCollection(refHapFreq.values());

				//Correlate study haplotypes to ref haplotypes.
				//Note: igonore case of no variance in both that would be correlation of 1 since no added value here

				double studyHapVar = JSci.maths.ArrayMath.variance(studyHapFreqArray);
				double refHapVar = JSci.maths.ArrayMath.variance(refHapFreqArray);

				double denom = Math.sqrt(studyHapVar * refHapVar);
				if (denom != 0) {
					double correlation = covariance(studyHapFreqArray, refHapFreqArray) / denom;

					if (correlation < 0) {
						++negCor;
					} else if (correlation > 0) {
						++posCor;
					}

				}


			}

		}

		return new CorrelationResults(posCor, negCor);
	}

	private double[] createDoubleArrayFromCollection(
			Collection<Double> values) {


		double[] array = new double[values.size()];

		int i = 0;
		for (Double d : values) {
			array[i] = d;
			++i;
		}

		return array;
	}

	private static class CorrelationResults {

		private final int posCor;
		private final int negCor;

		public CorrelationResults(int posCor, int negCor) {
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
