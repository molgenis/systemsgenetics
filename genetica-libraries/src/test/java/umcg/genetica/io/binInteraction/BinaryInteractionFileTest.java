/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package umcg.genetica.io.binInteraction;

import java.io.File;
import java.io.FileNotFoundException;
import java.io.IOException;
import java.text.DateFormat;
import java.text.SimpleDateFormat;
import java.util.Date;
import java.util.Iterator;
import org.molgenis.genotype.Allele;
import org.testng.annotations.Test;
import umcg.genetica.io.binInteraction.gene.BinaryInteractionGeneCreator;
import umcg.genetica.io.binInteraction.variant.BinaryInteractionVariantCreator;
import static org.testng.Assert.*;
import org.testng.annotations.BeforeTest;
import umcg.genetica.io.binInteraction.gene.BinaryInteractionGene;
import umcg.genetica.io.binInteraction.variant.BinaryInteractionVariant;

/**
 *
 * @author Patrick Deelen
 */
public class BinaryInteractionFileTest {

	private File tmpOutputFolder;
	
	public BinaryInteractionFileTest() {
	}
	
	@BeforeTest
	public void setUpMethod() throws Exception {
		File tmpDir = new File(System.getProperty("java.io.tmpdir"));

		DateFormat dateFormat = new SimpleDateFormat("yyyy-MM-dd_HH-mm-ss");
		Date date = new Date();

		tmpOutputFolder = new File(tmpDir, "BinaryInteractionFileTest_" + dateFormat.format(date));

		Runtime.getRuntime().addShutdownHook(new Thread() {
			@Override
			public void run() {
				for (File file : tmpOutputFolder.listFiles()) {
						file.delete();
				}
				tmpOutputFolder.delete();
			}
		});

		tmpOutputFolder.mkdir();


	}

	@Test
	public void test1() throws BinaryInteractionFileException, FileNotFoundException, IOException {
		
		File file = new File(tmpOutputFolder, "testInt1.bin");
		BinaryInteractionCohort[] cohorts = new BinaryInteractionCohort[1];
		BinaryInteractionGeneCreator[] genes = new BinaryInteractionGeneCreator[1];
		BinaryInteractionVariantCreator[] variants = new BinaryInteractionVariantCreator[1];
		String[] covariates = new String[1];

		cohorts[0] = new BinaryInteractionCohort("cohort1", 12);
		genes[0] = new BinaryInteractionGeneCreator("Gene1", "Chr1", 10, 100);
		variants[0] = new BinaryInteractionVariantCreator("Var1", "Chr2", 4, Allele.A, Allele.T);
		covariates[0] = "CellCount";

		BinaryInteractionFileCreator creator = new BinaryInteractionFileCreator(file, variants, genes, cohorts, covariates, true, true, true, true);

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

		assertTrue(createdInteractions.areAllCovariatesTestedForAllVariantGenes());
		assertTrue(loadedInteractions.areAllCovariatesTestedForAllVariantGenes());

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

		assertEquals(createdInteractions.getCovariates().get(0), "CellCount");
		assertEquals(loadedInteractions.getCovariates().get(0), "CellCount");

		assertEquals(createdInteractions.getTotalNumberInteractions(), 1l);
		assertEquals(loadedInteractions.getTotalNumberInteractions(), 1l);

		Iterator<BinaryInteractionQueryResult> iterator = loadedInteractions.readVariantGeneResults("Var1", "Gene1");
		BinaryInteractionQueryResult interaction = iterator.next();
		assertEquals(interaction.getGeneName(), "Gene1");
		assertEquals(interaction.getVariantName(), "Var1");
		assertEquals(interaction.getCovariateName(), "CellCount");
		assertEqualsInteractionZscores(interaction.getInteractionZscores(), interactionZscores);
		assertEqualsQtlZscores(interaction.getQtlZscores(), qtlZscore);
		assertFalse(iterator.hasNext());
		
		createdInteractions.close();
		loadedInteractions.close();

	}
	
	@Test
	public void test2() throws BinaryInteractionFileException, FileNotFoundException, IOException {

		File file = new File(tmpOutputFolder, "testInt2.bin");
		BinaryInteractionCohort[] cohorts = new BinaryInteractionCohort[2];
		BinaryInteractionGeneCreator[] genes = new BinaryInteractionGeneCreator[2];
		BinaryInteractionVariantCreator[] variants = new BinaryInteractionVariantCreator[2];
		String[] covariates = new String[2];

		cohorts[0] = new BinaryInteractionCohort("cohort1", 12);
		cohorts[1] = new BinaryInteractionCohort("cohort2", 60000);
		genes[0] = new BinaryInteractionGeneCreator("Gene1", "Chr1", 10, 100);
		genes[1] = new BinaryInteractionGeneCreator("Gene2", "ChrX", 5354467, 5351467);
		variants[0] = new BinaryInteractionVariantCreator("Var1", "Chr2", 4, Allele.A, Allele.T);
		variants[1] = new BinaryInteractionVariantCreator("Var2", "Chr2", 5, Allele.create("AT"), Allele.T);
		covariates[0] = "CellCount";
		covariates[1] = "Age";
		
		String[] covVar1Gene1 = {"CellCount", "Age"};
		String[] covVar2Gene1 = {"CellCount"};
		String[] covVar1Gene2 = {"Age"};

		BinaryInteractionFileCreator creator = new BinaryInteractionFileCreator(file, variants, genes, cohorts, covariates, false, true, true, false);

		creator.addTestedVariantGene("Var1", "Gene1");
		creator.addTestedVariantGene("Var2", "Gene1");
		creator.addTestedVariantGene("Var1", "Gene2");
		
		creator.setDescription("Test file 2");
		
		creator.addTestedInteraction("Var1", "Gene1", covVar1Gene1);
		creator.addTestedInteraction("Var2", "Gene1", covVar2Gene1);
		creator.addTestedInteraction("Var1", "Gene2", covVar1Gene2);
		
		BinaryInteractionFile createdInteractions = creator.create();

		assertFalse(createdInteractions.isReadOnly());

		double[] zscores = {1d,1000d};
		int[] samples = {11,50000};
		double metaZ = 1002d;
		BinaryInteractionQtlZscores qtlZscore1 = new BinaryInteractionQtlZscores(zscores, samples, metaZ);
		
		double[] zscores2 = {2d, 1001d};
		int[] samples2 = {12, 50001};
		double metaZ2 = 1003d;
		BinaryInteractionQtlZscores qtlZscore2 = new BinaryInteractionQtlZscores(zscores2, samples2, metaZ2);
		
		double[] zscores3 = {3d,1002d};
		int[] samples3 = {12,12434};
		double metaZ3 = 1003d;
		BinaryInteractionQtlZscores qtlZscore3 = new BinaryInteractionQtlZscores(zscores3, samples3, metaZ3);

		createdInteractions.setQtlResults("Var1", "Gene1", qtlZscore1);
		createdInteractions.setQtlResults("Var2", "Gene1", qtlZscore2);
		createdInteractions.setQtlResults("Var1", "Gene2", qtlZscore3);
		assertEqualsQtlZscores(createdInteractions.readQtlResults("Var1", "Gene1"), qtlZscore1);
		assertEqualsQtlZscores(createdInteractions.readQtlResults("Var2", "Gene1"), qtlZscore2);
		assertEqualsQtlZscores(createdInteractions.readQtlResults("Var1", "Gene2"), qtlZscore3);
		

		final int[] samplesInteractionCohort = {10,10000};
		final double[] zscoreSnpCohort = {3d, 5000d};
		final double[] zscoreCovariateCohort = {4d, 1000d};
		final double[] zscoreInteractionCohort = {5d, 2000d};
		final double[] rSquaredCohort = {0.987654321d, 0.999999999};
		final double zscoreSnpMeta = 106d;
		final double zscoreCovariateMeta = 107d;
		final double zscoreInteractionMeta = 108d;
		BinaryInteractionZscores interactionZscores1 = new BinaryInteractionZscores(samplesInteractionCohort, zscoreSnpCohort, zscoreCovariateCohort, zscoreInteractionCohort, rSquaredCohort, zscoreSnpMeta, zscoreCovariateMeta, zscoreInteractionMeta);
		
		createdInteractions.setInteractionResults("Var1", "Gene1", "CellCount", interactionZscores1);
		
		final int[] samplesInteractionCohort2 = {11,10001};
		final double[] zscoreSnpCohort2 = {4d, 6000d};
		final double[] zscoreCovariateCohort2 = {5d, 1001d};
		final double[] zscoreInteractionCohort2 = {6d, 2001d};
		final double[] rSquaredCohort2 = {0.98765432d, 0.999999998};
		final double zscoreSnpMeta2 = 107d;
		final double zscoreCovariateMeta2 = 108d;
		final double zscoreInteractionMeta2 = 109d;
		BinaryInteractionZscores interactionZscores2 = new BinaryInteractionZscores(samplesInteractionCohort2, zscoreSnpCohort2, zscoreCovariateCohort2, zscoreInteractionCohort2, rSquaredCohort2, zscoreSnpMeta2, zscoreCovariateMeta2, zscoreInteractionMeta2);
		
		createdInteractions.setInteractionResults("Var1", "Gene1", "Age", interactionZscores2);
		
		final int[] samplesInteractionCohort4 = {8,8};
		final double[] zscoreSnpCohort4 = {7d, 5063d};
		final double[] zscoreCovariateCohort4 = {3d, 10302d};
		final double[] zscoreInteractionCohort4 = {82d, 20302d};
		final double[] rSquaredCohort4 = {0.987651d, 0.9999};
		final double zscoreSnpMeta4 = 103434d;
		final double zscoreCovariateMeta4 = 1021349.2d;
		final double zscoreInteractionMeta4 = 20128.21431243d;
		BinaryInteractionZscores interactionZscores4 = new BinaryInteractionZscores(samplesInteractionCohort4, zscoreSnpCohort4, zscoreCovariateCohort4, zscoreInteractionCohort4, rSquaredCohort4, zscoreSnpMeta4, zscoreCovariateMeta4, zscoreInteractionMeta4);
		
		createdInteractions.setInteractionResults("Var2", "Gene1", "CellCount", interactionZscores4);
		
		final int[] samplesInteractionCohort3 = {9,9999};
		final double[] zscoreSnpCohort3 = {5d, 5003d};
		final double[] zscoreCovariateCohort3 = {7d, 1002d};
		final double[] zscoreInteractionCohort3 = {8d, 2002d};
		final double[] rSquaredCohort3 = {0.987654321d, 0.999999999};
		final double zscoreSnpMeta3 = 105d;
		final double zscoreCovariateMeta3 = 109d;
		final double zscoreInteractionMeta3 = 208d;
		BinaryInteractionZscores interactionZscores3 = new BinaryInteractionZscores(samplesInteractionCohort3, zscoreSnpCohort3, zscoreCovariateCohort3, zscoreInteractionCohort3, rSquaredCohort3, zscoreSnpMeta3, zscoreCovariateMeta3, zscoreInteractionMeta3);
		
		createdInteractions.setInteractionResults("Var1", "Gene2", "Age", interactionZscores3);
		
		
		
		assertEqualsInteractionZscores(createdInteractions.readInteractionResults("Var1", "Gene1", "CellCount"), interactionZscores1);
		assertEqualsInteractionZscores(createdInteractions.readInteractionResults("Var1", "Gene2", "Age"), interactionZscores3);
		assertEqualsInteractionZscores(createdInteractions.readInteractionResults("Var1", "Gene1", "Age"), interactionZscores2);
		assertEqualsInteractionZscores(createdInteractions.readInteractionResults("Var2", "Gene1", "CellCount"), interactionZscores4);
				
		assertEqualsInteractionZscores(createdInteractions.readInteractionResults("Var1", "Gene2", "Age"), interactionZscores3);
		assertEqualsInteractionZscores(createdInteractions.readInteractionResults("Var1", "Gene1", "Age"), interactionZscores2);
		assertEqualsInteractionZscores(createdInteractions.readInteractionResults("Var1", "Gene2", "Age"), interactionZscores3);
		assertEqualsInteractionZscores(createdInteractions.readInteractionResults("Var1", "Gene1", "CellCount"), interactionZscores1);
		
		createdInteractions.finalizeWriting();
		
		
		BinaryInteractionFile loadedInteractions = BinaryInteractionFile.load(file);	
		
		assertEquals(loadedInteractions.getFileDescription(), "Test file 2");
		assertEquals(createdInteractions.getFileDescription(), "Test file 2");

		assertEquals(createdInteractions.getCreationDataEpoch(), loadedInteractions.getCreationDataEpoch());

		assertFalse(createdInteractions.areAllCovariatesTestedForAllVariantGenes());
		assertFalse(loadedInteractions.areAllCovariatesTestedForAllVariantGenes());

		assertFalse(createdInteractions.isFlippedZscoreStored());
		assertFalse(loadedInteractions.isFlippedZscoreStored());

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
		assertEquals(originalVariant.getGeneCount(), 2);
		assertEquals(originalVariant.getGenePointers()[0], 0);
		assertEquals(originalVariant.getGenePointers()[1], 1);
		assertEquals(originalVariant.getGenePointers().length, 2);

		assertEquals(loadedVariant.getName(), "Var1");
		assertEquals(loadedVariant.getChr(), "Chr2");
		assertEquals(loadedVariant.getPos(), 4);
		assertEquals(loadedVariant.getRefAllele(), Allele.A);
		assertEquals(loadedVariant.getAltAllele(), Allele.T);
		assertEquals(loadedVariant.getGeneCount(), 2);
		assertEquals(loadedVariant.getGenePointers()[0], 0);
		assertEquals(loadedVariant.getGenePointers()[1], 1);
		assertEquals(loadedVariant.getGenePointers().length, 2);

		//variants[1] = new BinaryInteractionVariantCreator("Var2", "Chr2", 5, Allele.create("AT"), Allele.T);
		originalVariant = createdInteractions.getVariants().get(1);
		loadedVariant = loadedInteractions.getVariants().get(1);

		assertEquals(originalVariant.getName(), "Var2");
		assertEquals(originalVariant.getChr(), "Chr2");
		assertEquals(originalVariant.getPos(), 5);
		assertEquals(originalVariant.getRefAllele(), Allele.create("AT"));
		assertEquals(originalVariant.getAltAllele(), Allele.T);
		assertEquals(originalVariant.getGeneCount(), 1);
		assertEquals(originalVariant.getGenePointers()[0], 0);
		assertEquals(originalVariant.getGenePointers().length, 1);

		assertEquals(loadedVariant.getName(), "Var2");
		assertEquals(loadedVariant.getChr(), "Chr2");
		assertEquals(loadedVariant.getPos(), 5);
		assertEquals(loadedVariant.getRefAllele(), Allele.create("AT"));
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
		assertEquals(originalGene.getVariantCount(), 2);
		assertEquals(originalGene.getVariantPointers()[0], 0);
		assertEquals(originalGene.getVariantPointers()[1], 1);
		assertEquals(originalGene.getVariantPointers().length, 2);

		assertEquals(loadedGene.getName(), "Gene1");
		assertEquals(loadedGene.getChr(), "Chr1");
		assertEquals(loadedGene.getStart(), 10);
		assertEquals(loadedGene.getEnd(), 100);
		assertEquals(loadedGene.getVariantCount(), 2);
		assertEquals(loadedGene.getVariantPointers()[0], 0);
		assertEquals(loadedGene.getVariantPointers()[1], 1);
		assertEquals(loadedGene.getVariantPointers().length, 2);
		
		originalGene = createdInteractions.getGenes().get(1);
		loadedGene = loadedInteractions.getGenes().get(1);
		
		assertEquals(originalGene.getName(), "Gene2");
		assertEquals(originalGene.getChr(), "ChrX");
		assertEquals(originalGene.getStart(), 5354467);
		assertEquals(originalGene.getEnd(), 5351467);
		assertEquals(originalGene.getVariantCount(), 1);
		assertEquals(originalGene.getVariantPointers()[0], 0);
		assertEquals(originalGene.getVariantPointers().length, 1);

		assertEquals(loadedGene.getName(), "Gene2");
		assertEquals(loadedGene.getChr(), "ChrX");
		assertEquals(loadedGene.getStart(), 5354467);
		assertEquals(loadedGene.getEnd(), 5351467);
		assertEquals(loadedGene.getVariantCount(), 1);
		assertEquals(loadedGene.getVariantPointers()[0], 0);
		assertEquals(loadedGene.getVariantPointers().length, 1);

		assertEquals(createdInteractions.getCovariates().get(0), "CellCount");
		assertEquals(loadedInteractions.getCovariates().get(0), "CellCount");
		
		assertEquals(createdInteractions.getCovariates().get(1), "Age");
		assertEquals(loadedInteractions.getCovariates().get(1), "Age");

		assertEquals(createdInteractions.getTotalNumberInteractions(), 4);
		assertEquals(loadedInteractions.getTotalNumberInteractions(), 4);
		
		Iterator<BinaryInteractionQueryResult> iterator = loadedInteractions.readVariantGeneResults("Var1", "Gene1");
		BinaryInteractionQueryResult interaction = iterator.next();
		assertEquals(interaction.getGeneName(), "Gene1");
		assertEquals(interaction.getVariantName(), "Var1");
		assertEquals(interaction.getCovariateName(), "CellCount");
		assertEqualsInteractionZscores(interaction.getInteractionZscores(), interactionZscores1);
		assertEqualsQtlZscores(interaction.getQtlZscores(), qtlZscore1);
		
		interaction = iterator.next();
		assertEquals(interaction.getGeneName(), "Gene1");
		assertEquals(interaction.getVariantName(), "Var1");
		assertEquals(interaction.getCovariateName(), "Age");
		assertEqualsInteractionZscores(interaction.getInteractionZscores(), interactionZscores2);
		assertEqualsQtlZscores(interaction.getQtlZscores(), qtlZscore1);
			
		iterator = loadedInteractions.readVariantGeneResults("Var2", "Gene1");
		interaction = iterator.next();
		assertEquals(interaction.getGeneName(), "Gene1");
		assertEquals(interaction.getVariantName(), "Var2");
		assertEquals(interaction.getCovariateName(), "CellCount");
		assertEqualsInteractionZscores(interaction.getInteractionZscores(), interactionZscores4);
		assertEqualsQtlZscores(interaction.getQtlZscores(), qtlZscore2);
		
		assertFalse(iterator.hasNext());
		
		iterator = loadedInteractions.readVariantGeneResults("Var1", "Gene2");
		interaction = iterator.next();
		assertEquals(interaction.getGeneName(), "Gene2");
		assertEquals(interaction.getVariantName(), "Var1");
		assertEquals(interaction.getCovariateName(), "Age");
		assertEqualsInteractionZscores(interaction.getInteractionZscores(), interactionZscores3);
		assertEqualsQtlZscores(interaction.getQtlZscores(), qtlZscore3);
		
		assertFalse(iterator.hasNext());
		
		createdInteractions.close();
		loadedInteractions.close();

	}
	
	@Test
	public void test3() throws BinaryInteractionFileException, FileNotFoundException, IOException {

		File file = new File(tmpOutputFolder, "testInt3.bin");
		BinaryInteractionCohort[] cohorts = new BinaryInteractionCohort[2];
		BinaryInteractionGeneCreator[] genes = new BinaryInteractionGeneCreator[2];
		BinaryInteractionVariantCreator[] variants = new BinaryInteractionVariantCreator[2];
		String[] covariates = new String[2];

		cohorts[0] = new BinaryInteractionCohort("cohort1", 12);
		cohorts[1] = new BinaryInteractionCohort("cohort2", 60000);
		genes[0] = new BinaryInteractionGeneCreator("Gene1", "Chr1", 10, 100);
		genes[1] = new BinaryInteractionGeneCreator("Gene2", "ChrX", 5354467, 5351467);
		variants[0] = new BinaryInteractionVariantCreator("Var1", "Chr2", 4, Allele.A, Allele.T);
		variants[1] = new BinaryInteractionVariantCreator("Var2", "Chr2", 5, Allele.create("AT"), Allele.T);
		covariates[0] = "CellCount";
		covariates[1] = "Age";
		
		String[] covVar1Gene1 = {"CellCount", "Age"};
		String[] covVar2Gene1 = {"CellCount"};
		String[] covVar1Gene2 = {"Age"};

		BinaryInteractionFileCreator creator = new BinaryInteractionFileCreator(file, variants, genes, cohorts, covariates, false, false, true, false);

		creator.addTestedVariantGene("Var1", "Gene1");
		creator.addTestedVariantGene("Var2", "Gene1");
		creator.addTestedVariantGene("Var1", "Gene2");
		
		creator.setDescription("Test file 3");
		
		creator.addTestedInteraction("Var1", "Gene1", covVar1Gene1);
		creator.addTestedInteraction("Var2", "Gene1", covVar2Gene1);
		creator.addTestedInteraction("Var1", "Gene2", covVar1Gene2);
		
		BinaryInteractionFile createdInteractions = creator.create();

		assertFalse(createdInteractions.isReadOnly());

		double[] zscores = {1d,1000d};
		int[] samples = {11,50000};
		BinaryInteractionQtlZscores qtlZscore1 = new BinaryInteractionQtlZscores(zscores, samples);
		
		double[] zscores2 = {2d, 1001d};
		int[] samples2 = {12, 50001};
		BinaryInteractionQtlZscores qtlZscore2 = new BinaryInteractionQtlZscores(zscores2, samples2);
		
		double[] zscores3 = {3d,1002d};
		int[] samples3 = {12,12434};
		BinaryInteractionQtlZscores qtlZscore3 = new BinaryInteractionQtlZscores(zscores3, samples3);

		createdInteractions.setQtlResults("Var1", "Gene1", qtlZscore1);
		createdInteractions.setQtlResults("Var2", "Gene1", qtlZscore2);
		createdInteractions.setQtlResults("Var1", "Gene2", qtlZscore3);
		assertEqualsQtlZscores(createdInteractions.readQtlResults("Var1", "Gene1"), qtlZscore1);
		assertEqualsQtlZscores(createdInteractions.readQtlResults("Var2", "Gene1"), qtlZscore2);
		assertEqualsQtlZscores(createdInteractions.readQtlResults("Var1", "Gene2"), qtlZscore3);
		

		final int[] samplesInteractionCohort = {10,10000};
		final double[] zscoreSnpCohort = {3d, 5000d};
		final double[] zscoreCovariateCohort = {4d, 1000d};
		final double[] zscoreInteractionCohort = {5d, 2000d};
		final double[] rSquaredCohort = {0.987654321d, 0.999999999};
		BinaryInteractionZscores interactionZscores1 = new BinaryInteractionZscores(samplesInteractionCohort, zscoreSnpCohort, zscoreCovariateCohort, zscoreInteractionCohort, rSquaredCohort);
		
		createdInteractions.setInteractionResults("Var1", "Gene1", "CellCount", interactionZscores1);
		
		final int[] samplesInteractionCohort2 = {11,10001};
		final double[] zscoreSnpCohort2 = {4d, 6000d};
		final double[] zscoreCovariateCohort2 = {5d, 1001d};
		final double[] zscoreInteractionCohort2 = {6d, 2001d};
		final double[] rSquaredCohort2 = {0.98765432d, 0.999999998};
		BinaryInteractionZscores interactionZscores2 = new BinaryInteractionZscores(samplesInteractionCohort2, zscoreSnpCohort2, zscoreCovariateCohort2, zscoreInteractionCohort2, rSquaredCohort2);
		
		createdInteractions.setInteractionResults("Var1", "Gene1", "Age", interactionZscores2);
		
		final int[] samplesInteractionCohort4 = {8,8};
		final double[] zscoreSnpCohort4 = {Double.MIN_NORMAL, Double.NEGATIVE_INFINITY};
		final double[] zscoreCovariateCohort4 = {Double.MAX_VALUE, Double.NaN};
		final double[] zscoreInteractionCohort4 = {Double.MIN_VALUE, Double.POSITIVE_INFINITY};
		final double[] rSquaredCohort4 = {0.987651d, 0.9999};
		BinaryInteractionZscores interactionZscores4 = new BinaryInteractionZscores(samplesInteractionCohort4, zscoreSnpCohort4, zscoreCovariateCohort4, zscoreInteractionCohort4, rSquaredCohort4);
		
		createdInteractions.setInteractionResults("Var2", "Gene1", "CellCount", interactionZscores4);
		
		final int[] samplesInteractionCohort3 = {9,9999};
		final double[] zscoreSnpCohort3 = {5d, 5003d};
		final double[] zscoreCovariateCohort3 = {7d, 1002d};
		final double[] zscoreInteractionCohort3 = {8d, 2002d};
		final double[] rSquaredCohort3 = {0.987654321d, 0.999999999};
		BinaryInteractionZscores interactionZscores3 = new BinaryInteractionZscores(samplesInteractionCohort3, zscoreSnpCohort3, zscoreCovariateCohort3, zscoreInteractionCohort3, rSquaredCohort3);
		
		createdInteractions.setInteractionResults("Var1", "Gene2", "Age", interactionZscores3);
		
		
		
		assertEqualsInteractionZscores(createdInteractions.readInteractionResults("Var1", "Gene1", "CellCount"), interactionZscores1);
		assertEqualsInteractionZscores(createdInteractions.readInteractionResults("Var1", "Gene2", "Age"), interactionZscores3);
		assertEqualsInteractionZscores(createdInteractions.readInteractionResults("Var1", "Gene1", "Age"), interactionZscores2);
		assertEqualsInteractionZscores(createdInteractions.readInteractionResults("Var2", "Gene1", "CellCount"), interactionZscores4);
		
		assertEqualsInteractionZscores(createdInteractions.readInteractionResults("Var1", "Gene2", "Age"), interactionZscores3);
		assertEqualsInteractionZscores(createdInteractions.readInteractionResults("Var1", "Gene1", "Age"), interactionZscores2);
		assertEqualsInteractionZscores(createdInteractions.readInteractionResults("Var1", "Gene2", "Age"), interactionZscores3);
		assertEqualsInteractionZscores(createdInteractions.readInteractionResults("Var1", "Gene1", "CellCount"), interactionZscores1);
		
		createdInteractions.finalizeWriting();
		
		
		BinaryInteractionFile loadedInteractions = BinaryInteractionFile.load(file);	
		
		assertEquals(loadedInteractions.getFileDescription(), "Test file 3");
		assertEquals(createdInteractions.getFileDescription(), "Test file 3");

		assertEquals(createdInteractions.getCreationDataEpoch(), loadedInteractions.getCreationDataEpoch());

		assertFalse(createdInteractions.areAllCovariatesTestedForAllVariantGenes());
		assertFalse(loadedInteractions.areAllCovariatesTestedForAllVariantGenes());

		assertFalse(createdInteractions.isFlippedZscoreStored());
		assertFalse(loadedInteractions.isFlippedZscoreStored());

		assertFalse(createdInteractions.isMetaAnalysis());
		assertFalse(loadedInteractions.isMetaAnalysis());

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
		assertEquals(originalVariant.getGeneCount(), 2);
		assertEquals(originalVariant.getGenePointers()[0], 0);
		assertEquals(originalVariant.getGenePointers()[1], 1);
		assertEquals(originalVariant.getGenePointers().length, 2);

		assertEquals(loadedVariant.getName(), "Var1");
		assertEquals(loadedVariant.getChr(), "Chr2");
		assertEquals(loadedVariant.getPos(), 4);
		assertEquals(loadedVariant.getRefAllele(), Allele.A);
		assertEquals(loadedVariant.getAltAllele(), Allele.T);
		assertEquals(loadedVariant.getGeneCount(), 2);
		assertEquals(loadedVariant.getGenePointers()[0], 0);
		assertEquals(loadedVariant.getGenePointers()[1], 1);
		assertEquals(loadedVariant.getGenePointers().length, 2);

		//variants[1] = new BinaryInteractionVariantCreator("Var2", "Chr2", 5, Allele.create("AT"), Allele.T);
		originalVariant = createdInteractions.getVariants().get(1);
		loadedVariant = loadedInteractions.getVariants().get(1);

		assertEquals(originalVariant.getName(), "Var2");
		assertEquals(originalVariant.getChr(), "Chr2");
		assertEquals(originalVariant.getPos(), 5);
		assertEquals(originalVariant.getRefAllele(), Allele.create("AT"));
		assertEquals(originalVariant.getAltAllele(), Allele.T);
		assertEquals(originalVariant.getGeneCount(), 1);
		assertEquals(originalVariant.getGenePointers()[0], 0);
		assertEquals(originalVariant.getGenePointers().length, 1);

		assertEquals(loadedVariant.getName(), "Var2");
		assertEquals(loadedVariant.getChr(), "Chr2");
		assertEquals(loadedVariant.getPos(), 5);
		assertEquals(loadedVariant.getRefAllele(), Allele.create("AT"));
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
		assertEquals(originalGene.getVariantCount(), 2);
		assertEquals(originalGene.getVariantPointers()[0], 0);
		assertEquals(originalGene.getVariantPointers()[1], 1);
		assertEquals(originalGene.getVariantPointers().length, 2);

		assertEquals(loadedGene.getName(), "Gene1");
		assertEquals(loadedGene.getChr(), "Chr1");
		assertEquals(loadedGene.getStart(), 10);
		assertEquals(loadedGene.getEnd(), 100);
		assertEquals(loadedGene.getVariantCount(), 2);
		assertEquals(loadedGene.getVariantPointers()[0], 0);
		assertEquals(loadedGene.getVariantPointers()[1], 1);
		assertEquals(loadedGene.getVariantPointers().length, 2);
		
		originalGene = createdInteractions.getGenes().get(1);
		loadedGene = loadedInteractions.getGenes().get(1);
		
		assertEquals(originalGene.getName(), "Gene2");
		assertEquals(originalGene.getChr(), "ChrX");
		assertEquals(originalGene.getStart(), 5354467);
		assertEquals(originalGene.getEnd(), 5351467);
		assertEquals(originalGene.getVariantCount(), 1);
		assertEquals(originalGene.getVariantPointers()[0], 0);
		assertEquals(originalGene.getVariantPointers().length, 1);

		assertEquals(loadedGene.getName(), "Gene2");
		assertEquals(loadedGene.getChr(), "ChrX");
		assertEquals(loadedGene.getStart(), 5354467);
		assertEquals(loadedGene.getEnd(), 5351467);
		assertEquals(loadedGene.getVariantCount(), 1);
		assertEquals(loadedGene.getVariantPointers()[0], 0);
		assertEquals(loadedGene.getVariantPointers().length, 1);

		assertEquals(createdInteractions.getCovariates().get(0), "CellCount");
		assertEquals(loadedInteractions.getCovariates().get(0), "CellCount");
		
		assertEquals(createdInteractions.getCovariates().get(1), "Age");
		assertEquals(loadedInteractions.getCovariates().get(1), "Age");

		assertEquals(createdInteractions.getTotalNumberInteractions(), 4);
		assertEquals(loadedInteractions.getTotalNumberInteractions(), 4);

		Iterator<BinaryInteractionQueryResult> iterator = loadedInteractions.readVariantGeneResults("Var1", "Gene1");
		BinaryInteractionQueryResult interaction = iterator.next();
		assertEquals(interaction.getGeneName(), "Gene1");
		assertEquals(interaction.getVariantName(), "Var1");
		assertEquals(interaction.getCovariateName(), "CellCount");
		assertEqualsInteractionZscores(interaction.getInteractionZscores(), interactionZscores1);
		assertEqualsQtlZscores(interaction.getQtlZscores(), qtlZscore1);
		
		interaction = iterator.next();
		assertEquals(interaction.getGeneName(), "Gene1");
		assertEquals(interaction.getVariantName(), "Var1");
		assertEquals(interaction.getCovariateName(), "Age");
		assertEqualsInteractionZscores(interaction.getInteractionZscores(), interactionZscores2);
		assertEqualsQtlZscores(interaction.getQtlZscores(), qtlZscore1);
			
		iterator = loadedInteractions.readVariantGeneResults("Var2", "Gene1");
		interaction = iterator.next();
		assertEquals(interaction.getGeneName(), "Gene1");
		assertEquals(interaction.getVariantName(), "Var2");
		assertEquals(interaction.getCovariateName(), "CellCount");
		assertEqualsInteractionZscores(interaction.getInteractionZscores(), interactionZscores4);
		assertEqualsQtlZscores(interaction.getQtlZscores(), qtlZscore2);
		
		assertFalse(iterator.hasNext());
		
		iterator = loadedInteractions.readVariantGeneResults("Var1", "Gene2");
		interaction = iterator.next();
		assertEquals(interaction.getGeneName(), "Gene2");
		assertEquals(interaction.getVariantName(), "Var1");
		assertEquals(interaction.getCovariateName(), "Age");
		assertEqualsInteractionZscores(interaction.getInteractionZscores(), interactionZscores3);
		assertEqualsQtlZscores(interaction.getQtlZscores(), qtlZscore3);
		
		assertFalse(iterator.hasNext());
		
		createdInteractions.close();
		loadedInteractions.close();

	}

	public void assertEqualsDoubleArray(double[] actual, double[] expected, double delta){
		
		assertEquals(actual.length, expected.length);
		
		for(int i = 0 ; i < actual.length ; ++i){
			assertEqualsNaN(actual[i], expected[i], delta);
		}
		
	}
	
	public void assertEqualsNaN(double actual, double expected, double delta){
		if(!Double.isNaN(actual) && Double.isNaN(expected)){
			assertEquals(actual, expected, delta);
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
		
		assertEqualsDoubleArray(actual.getZscores(), expected.getZscores(), 1e-10);
		assertEqualsIntArray(actual.getSampleCounts(), expected.getSampleCounts());
		assertEqualsNaN(actual.getMetaZscore(), expected.getMetaZscore(), 1e-10);
		
	}
	
	public void assertEqualsInteractionZscores(BinaryInteractionZscores actual, BinaryInteractionZscores expected){
		
		assertEquals(actual.getCohortCount(), expected.getCohortCount());
		
		assertEqualsDoubleArray(actual.getZscoreCovariateCohort(), expected.getZscoreCovariateCohort(), 1e-10);
		assertEqualsDoubleArray(actual.getZscoreInteractionCohort(), expected.getZscoreInteractionCohort(), 1e-10);
		assertEqualsDoubleArray(actual.getZscoreInteractionFlippedCohort(), expected.getZscoreInteractionFlippedCohort(), 1e-10);
		assertEqualsDoubleArray(actual.getZscoreSnpCohort(), expected.getZscoreSnpCohort(), 1e-10);
		assertEqualsDoubleArray(actual.getrSquaredCohort(), expected.getrSquaredCohort(), 1e-10);
		
		assertEqualsIntArray(actual.getSamplesInteractionCohort(), expected.getSamplesInteractionCohort());
		
		assertEqualsNaN(actual.getZscoreCovariateMeta(), expected.getZscoreCovariateMeta(), 1e-10);
		assertEqualsNaN(actual.getZscoreInteractionFlippedMeta(), expected.getZscoreInteractionFlippedMeta(), 1e-10);
		assertEqualsNaN(actual.getZscoreInteractionMeta(), expected.getZscoreInteractionMeta(), 1e-10);
		assertEqualsNaN(actual.getZscoreSnpMeta(), expected.getZscoreSnpMeta(), 1e-10);
		
	}
	
}