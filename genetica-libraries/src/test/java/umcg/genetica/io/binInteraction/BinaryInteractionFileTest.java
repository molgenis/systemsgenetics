/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package umcg.genetica.io.binInteraction;

import java.io.File;
import java.io.FileNotFoundException;
import java.io.IOException;
import org.apache.mahout.math.Arrays;
import org.molgenis.genotype.Allele;
import org.testng.annotations.Test;
import umcg.genetica.io.binInteraction.gene.BinaryInteractionGeneCreator;
import umcg.genetica.io.binInteraction.variant.BinaryInteractionVariantCreator;
import static org.testng.Assert.*;
import umcg.genetica.io.binInteraction.gene.BinaryInteractionGene;
import umcg.genetica.io.binInteraction.variant.BinaryInteractionVariant;

/**
 *
 * @author Patrick Deelen
 */
public class BinaryInteractionFileTest {

	public BinaryInteractionFileTest() {
	}

	@Test
	public void test1() throws BinaryInteractionFileException, FileNotFoundException, IOException {

		File file = new File("D:\\tmp\\testInt.bin");
		BinaryInteractionCohort[] cohorts = new BinaryInteractionCohort[1];
		BinaryInteractionGeneCreator[] genes = new BinaryInteractionGeneCreator[1];
		BinaryInteractionVariantCreator[] variants = new BinaryInteractionVariantCreator[1];
		String[] covariats = new String[1];

		cohorts[0] = new BinaryInteractionCohort("cohort1", 12);
		genes[0] = new BinaryInteractionGeneCreator("Gene1", "Chr1", 10, 100);
		variants[0] = new BinaryInteractionVariantCreator("Var1", "Chr2", 4, Allele.A, Allele.T);
		covariats[0] = "CellCount";

		BinaryInteractionFileCreator creator = new BinaryInteractionFileCreator(file, variants, genes, cohorts, covariats, true, true, true, true);

		creator.addTestedVariantGene("Var1", "Gene1");
		creator.setDescription("Test file 1");
		BinaryInteractionFile createdInteractions = creator.create();

		assertFalse(createdInteractions.isReadOnly());

		double[] zscores = {1d};
		int[] samples = {11};
		double metaZ = 2d;

		createdInteractions.setQtlResults("Var1", "Gene1", new BinaryInteractionQtlZscores(zscores, samples, metaZ));


		final int[] samplesInteractionCohort = {10};
		final double[] zscoreSnpCohort = {3d};
		final double[] zscoreCovariateCohort = {4d};
		final double[] zscoreInteractionCohort = {5d};
		final double[] rSquaredCohort = {0.987654321d};
		final double[] zscoreInteractionFlippedCohort = {-5d};
		final double zscoreSnpMeta = 6d;
		final double zscoreCovariateMeta = 1E300d;
		final double zscoreInteractionMeta = 8d;
		final double zscoreInteractionFlippedMeta = -8d;

		createdInteractions.setInteractionResults(
				"Var1", "Gene1", "CellCount", new BinaryInteractionZscores(samplesInteractionCohort, zscoreSnpCohort, zscoreCovariateCohort, zscoreInteractionCohort, rSquaredCohort, zscoreInteractionFlippedCohort, zscoreSnpMeta, zscoreCovariateMeta, zscoreInteractionMeta, zscoreInteractionFlippedMeta));

//		System.out.println(Arrays.toString(createdInteractions.readQtlResults("Var1", "Gene1").getZscore()));
//		System.out.println(Arrays.toString(createdInteractions.readQtlResults("Var1", "Gene1").getSampleCounts()));
//		System.out.println(createdInteractions.readQtlResults("Var1", "Gene1").getMetaZscore());

		createdInteractions.finalizeWriting();
		BinaryInteractionFile loadedInteractions = BinaryInteractionFile.load(file);

		System.out.println(Arrays.toString(loadedInteractions.readQtlResults("Var1", "Gene1").getZscore()));
		System.out.println(Arrays.toString(loadedInteractions.readQtlResults("Var1", "Gene1").getSampleCounts()));
		System.out.println(loadedInteractions.readQtlResults("Var1", "Gene1").getMetaZscore());
		
		BinaryInteractionZscores interaction = loadedInteractions.readInteractionResults("Var1", "Gene1", "CellCount");
		System.out.println("-----");
		System.out.println(Arrays.toString(interaction.getSamplesInteractionCohort()));
		System.out.println(Arrays.toString(interaction.getZscoreSnpCohort()));
		System.out.println(Arrays.toString(interaction.getZscoreCovariateCohort()));
		System.out.println(Arrays.toString(interaction.getZscoreInteractionCohort()));
		System.out.println(Arrays.toString(interaction.getZscoreInteractionFlippedCohort()));
		System.out.println(Arrays.toString(interaction.getrSquaredCohort()));
		System.out.println(interaction.getZscoreSnpMeta());
		System.out.println(interaction.getZscoreCovariateMeta());
		System.out.println(interaction.getZscoreInteractionMeta());
		System.out.println(interaction.getZscoreInteractionFlippedMeta());

		assertEquals(loadedInteractions.getFileDescription(), "Test file 1");
		assertEquals(createdInteractions.getFileDescription(), "Test file 1");

		assertEquals(createdInteractions.getCreationDataEpoch(), loadedInteractions.getCreationDataEpoch());

		assertTrue(createdInteractions.areAllCovariatsTestedForAllVariantGenes());
		assertTrue(loadedInteractions.areAllCovariatsTestedForAllVariantGenes());

		assertTrue(createdInteractions.isFlippedZscoreStored());
		assertTrue(loadedInteractions.isFlippedZscoreStored());

		assertTrue(createdInteractions.isMetaAnalysis());
		assertTrue(loadedInteractions.isMetaAnalysis());

		assertTrue(createdInteractions.isNormalQtlStored());
		assertTrue(loadedInteractions.isNormalQtlStored());

		assertTrue(createdInteractions.isReadOnly());
		assertTrue(loadedInteractions.isReadOnly());

		assertEquals(loadedInteractions.getCohorts().get(0).getName(), "cohort1");
		assertEquals(createdInteractions.getCohorts().get(0).getName(), "cohort1");

		assertEquals(loadedInteractions.getCohorts().get(0).getSampleCount(), 12);
		assertEquals(createdInteractions.getCohorts().get(0).getSampleCount(), 12);

		BinaryInteractionVariant originalVariant = createdInteractions.getVariants().get(0);
		BinaryInteractionVariant loadedVariant = loadedInteractions.getVariants().get(0);

		assertEquals(originalVariant.getName(), "Var1");
		assertEquals(originalVariant.getChr(), "Chr2");
		assertEquals(originalVariant.getPos(), 4);
		assertEquals(originalVariant.getRefAllele(), Allele.A);
		assertEquals(originalVariant.getAltAllele(), Allele.T);
		assertEquals(originalVariant.getGeneCount(), 1);
		assertEquals(originalVariant.getGenePointers()[0], 0);
		assertEquals(originalVariant.getGenePointers().length, 1);

		assertEquals(loadedVariant.getName(), "Var1");
		assertEquals(loadedVariant.getChr(), "Chr2");
		assertEquals(loadedVariant.getPos(), 4);
		assertEquals(loadedVariant.getRefAllele(), Allele.A);
		assertEquals(loadedVariant.getAltAllele(), Allele.T);
		assertEquals(loadedVariant.getGeneCount(), 1);
		assertEquals(loadedVariant.getGenePointers()[0], 0);
		assertEquals(loadedVariant.getGenePointers().length, 1);

		BinaryInteractionGene originalGene = createdInteractions.getGenes().get(0);
		BinaryInteractionGene loadedGene = loadedInteractions.getGenes().get(0);

		assertEquals(originalGene.getName(), "Gene1");
		assertEquals(originalGene.getChr(), "Chr1");
		assertEquals(originalGene.getStart(), 10);
		assertEquals(originalGene.getEnd(), 100);
		assertEquals(originalGene.getVariantCount(), 1);
		assertEquals(originalGene.getVariantPointers()[0], 0);
		assertEquals(originalGene.getVariantPointers().length, 1);

		assertEquals(loadedGene.getName(), "Gene1");
		assertEquals(loadedGene.getChr(), "Chr1");
		assertEquals(loadedGene.getStart(), 10);
		assertEquals(loadedGene.getEnd(), 100);
		assertEquals(loadedGene.getVariantCount(), 1);
		assertEquals(loadedGene.getVariantPointers()[0], 0);
		assertEquals(loadedGene.getVariantPointers().length, 1);

		assertEquals(createdInteractions.getCovariats().get(0), "CellCount");
		assertEquals(loadedInteractions.getCovariats().get(0), "CellCount");

		assertEquals(createdInteractions.getTotalNumberInteractions(), 1l);
		assertEquals(loadedInteractions.getTotalNumberInteractions(), 1l);


	}
}