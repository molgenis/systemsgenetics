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
		BinaryInteractionQtlZscores qtlZscore = new BinaryInteractionQtlZscores(zscores, samples, metaZ);

		createdInteractions.setQtlResults("Var1", "Gene1", qtlZscore);
		assertEqualsQtlZscores(createdInteractions.readQtlResults("Var1", "Gene1"), qtlZscore);

		final int[] samplesInteractionCohort = {10};
		final double[] zscoreSnpCohort = {3d};
		final double[] zscoreCovariateCohort = {4d};
		final double[] zscoreInteractionCohort = {5d};
		final double[] rSquaredCohort = {0.987654321d};
		final double[] zscoreInteractionFlippedCohort = {-5d};
		final double zscoreSnpMeta = 6d;
		final double zscoreCovariateMeta = 1E301d;
		final double zscoreInteractionMeta = 8d;
		final double zscoreInteractionFlippedMeta = -8d;
		BinaryInteractionZscores interactionZscores = new BinaryInteractionZscores(samplesInteractionCohort, zscoreSnpCohort, zscoreCovariateCohort, zscoreInteractionCohort, rSquaredCohort, zscoreInteractionFlippedCohort, zscoreSnpMeta, zscoreCovariateMeta, zscoreInteractionMeta, zscoreInteractionFlippedMeta);
		
		createdInteractions.setInteractionResults("Var1", "Gene1", "CellCount", interactionZscores);
		assertEqualsInteractionZscores(createdInteractions.readInteractionResults("Var1", "Gene1", "CellCount"), interactionZscores);
		
		createdInteractions.finalizeWriting();
		
		assertEqualsQtlZscores(createdInteractions.readQtlResults("Var1", "Gene1"), qtlZscore);
		assertEqualsInteractionZscores(createdInteractions.readInteractionResults("Var1", "Gene1", "CellCount"), interactionZscores);
		
		BinaryInteractionFile loadedInteractions = BinaryInteractionFile.load(file);

		assertEqualsQtlZscores(loadedInteractions.readQtlResults("Var1", "Gene1"), qtlZscore);
		assertEqualsInteractionZscores(loadedInteractions.readInteractionResults("Var1", "Gene1", "CellCount"), interactionZscores);	
		
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
	
	@Test
	public void test2() throws BinaryInteractionFileException, FileNotFoundException, IOException {

		File file = new File("D:\\tmp\\testInt2.bin");
		BinaryInteractionCohort[] cohorts = new BinaryInteractionCohort[2];
		BinaryInteractionGeneCreator[] genes = new BinaryInteractionGeneCreator[2];
		BinaryInteractionVariantCreator[] variants = new BinaryInteractionVariantCreator[2];
		String[] covariats = new String[2];

		cohorts[0] = new BinaryInteractionCohort("cohort1", 12);
		cohorts[1] = new BinaryInteractionCohort("cohort2", 60000);
		genes[0] = new BinaryInteractionGeneCreator("Gene1", "Chr1", 10, 100);
		genes[1] = new BinaryInteractionGeneCreator("Gene2", "ChrX", 5354467, 5351467);
		variants[0] = new BinaryInteractionVariantCreator("Var1", "Chr2", 4, Allele.A, Allele.T);
		variants[1] = new BinaryInteractionVariantCreator("Var2", "Chr2", 5, Allele.create("AT"), Allele.T);
		covariats[0] = "CellCount";
		covariats[1] = "Age";
		
		String[] covVar1Gene1 = {"CellCount", "Age"};
		String[] covVar2Gene1 = {"CellCount"};
		String[] covVar1Gene2 = {"Age"};

		BinaryInteractionFileCreator creator = new BinaryInteractionFileCreator(file, variants, genes, cohorts, covariats, false, true, true, false);

		creator.addTestedVariantGene("Var1", "Gene1");
		creator.addTestedVariantGene("Var2", "Gene1");
		creator.addTestedVariantGene("Var1", "Gene2");
		
		creator.setDescription("Test file 2");
		
		creator.addTestedInteraction("Var1", "Gene1", covVar1Gene1);
		creator.addTestedInteraction("Var2", "Gene1", covVar2Gene1);
		creator.addTestedInteraction("Var1", "Gene2", covVar1Gene2);
		
		BinaryInteractionFile createdInteractions = creator.create();

		assertFalse(createdInteractions.isReadOnly());

		double[] zscores = {1d};
		int[] samples = {11};
		double metaZ = 2d;
		BinaryInteractionQtlZscores qtlZscore = new BinaryInteractionQtlZscores(zscores, samples, metaZ);

		createdInteractions.setQtlResults("Var1", "Gene1", qtlZscore);
		assertEqualsQtlZscores(createdInteractions.readQtlResults("Var1", "Gene1"), qtlZscore);

		final int[] samplesInteractionCohort = {10};
		final double[] zscoreSnpCohort = {3d};
		final double[] zscoreCovariateCohort = {4d};
		final double[] zscoreInteractionCohort = {5d};
		final double[] rSquaredCohort = {0.987654321d};
		final double[] zscoreInteractionFlippedCohort = {-5d};
		final double zscoreSnpMeta = 6d;
		final double zscoreCovariateMeta = 1E301d;
		final double zscoreInteractionMeta = 8d;
		final double zscoreInteractionFlippedMeta = -8d;
		BinaryInteractionZscores interactionZscores = new BinaryInteractionZscores(samplesInteractionCohort, zscoreSnpCohort, zscoreCovariateCohort, zscoreInteractionCohort, rSquaredCohort, zscoreInteractionFlippedCohort, zscoreSnpMeta, zscoreCovariateMeta, zscoreInteractionMeta, zscoreInteractionFlippedMeta);
		
		createdInteractions.setInteractionResults("Var1", "Gene1", "CellCount", interactionZscores);
		assertEqualsInteractionZscores(createdInteractions.readInteractionResults("Var1", "Gene1", "CellCount"), interactionZscores);
		
		createdInteractions.finalizeWriting();
		
		assertEqualsQtlZscores(createdInteractions.readQtlResults("Var1", "Gene1"), qtlZscore);
		assertEqualsInteractionZscores(createdInteractions.readInteractionResults("Var1", "Gene1", "CellCount"), interactionZscores);
		
		BinaryInteractionFile loadedInteractions = BinaryInteractionFile.load(file);

		assertEqualsQtlZscores(loadedInteractions.readQtlResults("Var1", "Gene1"), qtlZscore);
		assertEqualsInteractionZscores(loadedInteractions.readInteractionResults("Var1", "Gene1", "CellCount"), interactionZscores);	
		
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

	public void assertEqualsDoubleArray(double[] actual, double[] expected, double delta){
		
		assertEquals(actual.length, expected.length);
		
		for(int i = 0 ; i < actual.length ; ++i){
			assertEquals(actual[i], expected[i], delta);
		}
		
	}
	
	public void assertEqualsIntArray(int[] actual, int[] expected){
		
		assertEquals(actual.length, expected.length);
		
		for(int i = 0 ; i < actual.length ; ++i){
			assertEquals(actual[i], expected[i]);
		}
		
	}
	
	public void assertEqualsQtlZscores(BinaryInteractionQtlZscores actual, BinaryInteractionQtlZscores expected){
		
		assertEquals(actual.getCohortCount(), expected.getCohortCount());
		
		assertEqualsDoubleArray(actual.getZscore(), expected.getZscore(), 1e-10);
		assertEqualsIntArray(actual.getSampleCounts(), expected.getSampleCounts());
		assertEquals(actual.getMetaZscore(), expected.getMetaZscore());
		
	}
	
	public void assertEqualsInteractionZscores(BinaryInteractionZscores actual, BinaryInteractionZscores expected){
		
		assertEquals(actual.getCohortCount(), expected.getCohortCount());
		
		assertEqualsDoubleArray(actual.getZscoreCovariateCohort(), expected.getZscoreCovariateCohort(), 1e-10);
		assertEqualsDoubleArray(actual.getZscoreInteractionCohort(), expected.getZscoreInteractionCohort(), 1e-10);
		assertEqualsDoubleArray(actual.getZscoreInteractionFlippedCohort(), expected.getZscoreInteractionFlippedCohort(), 1e-10);
		assertEqualsDoubleArray(actual.getZscoreSnpCohort(), expected.getZscoreSnpCohort(), 1e-10);
		assertEqualsDoubleArray(actual.getrSquaredCohort(), expected.getrSquaredCohort(), 1e-10);
		
		assertEqualsIntArray(actual.getSamplesInteractionCohort(), expected.getSamplesInteractionCohort());
		
		assertEquals(actual.getZscoreCovariateMeta(), expected.getZscoreCovariateMeta());
		assertEquals(actual.getZscoreInteractionFlippedMeta(), expected.getZscoreInteractionFlippedMeta());
		assertEquals(actual.getZscoreInteractionMeta(), expected.getZscoreInteractionMeta());
		assertEquals(actual.getZscoreSnpMeta(), expected.getZscoreSnpMeta());
		
	}
	
}