/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package umcg.genetica.io.binInteraction;

import java.io.File;
import java.io.FileNotFoundException;
import java.io.IOException;
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
	public void test1() throws BinaryInteractionFileException, FileNotFoundException, IOException{
		
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
		
		createdInteractions.makeReadOnly();
		
		BinaryInteractionFile loadedInteractions = BinaryInteractionFile.load(file);
		
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
		
		
		
	}


}