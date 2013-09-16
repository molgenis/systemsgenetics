/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package nl.umcg.deelenp.genotypealigner;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileInputStream;
import java.io.InputStreamReader;
import java.net.URISyntaxException;
import java.nio.charset.Charset;
import java.text.DateFormat;
import java.text.SimpleDateFormat;
import java.util.Date;
import java.util.HashSet;
import java.util.Iterator;
import org.molgenis.genotype.Allele;
import org.molgenis.genotype.Alleles;
import org.molgenis.genotype.RandomAccessGenotypeData;
import org.molgenis.genotype.plink.BedBimFamGenotypeData;
import org.molgenis.genotype.variant.GeneticVariant;
import static org.testng.Assert.*;
import org.testng.annotations.BeforeTest;
import org.testng.annotations.Test;

/**
 *
 * @author Patrick Deelen
 */
public class GenotypeAlignerNGTest {

	private File tmpOutputFolder;
	private String fileSep = System.getProperty("file.separator");
	private File testFilesFolder;
	private static final Charset FILE_ENCODING = Charset.forName("UTF-8");
	private static final Alleles MISSING_ALLELES = Alleles.createAlleles(Allele.ZERO, Allele.ZERO);

	public GenotypeAlignerNGTest() throws URISyntaxException {
		testFilesFolder = new File(this.getClass().getResource("/").toURI());
	}

	@BeforeTest
	public void setUpMethod() throws Exception {
		File tmpDir = new File(System.getProperty("java.io.tmpdir"));

		DateFormat dateFormat = new SimpleDateFormat("yyyy-MM-dd_HH-mm-ss");
		Date date = new Date();

		tmpOutputFolder = new File(tmpDir, "GenotyperAligerTest_" + dateFormat.format(date));

		Runtime.getRuntime().addShutdownHook(new Thread() {
			@Override
			public void run() {
				System.out.println("Removing tmp dir and files");
				for (File file : tmpOutputFolder.listFiles()) {
					System.out.println(" - Deleting: " + file.getAbsolutePath());
					file.delete();
				}
				System.out.println(" - Deleting: " + tmpOutputFolder.getAbsolutePath());
				tmpOutputFolder.delete();
			}
		});

		tmpOutputFolder.mkdir();


		System.out.println("Temp folder with output of this test: " + tmpOutputFolder.getAbsolutePath());
	}

	/**
	 * Test of main method, of class GenotypeAligner.
	 */
	@Test
	public void testMain() throws Exception {
		System.out.println("main");

		String studyDataBasePath = testFilesFolder + fileSep + "hapmap3CeuChr20B37Mb6RandomStrand";
		System.out.println(studyDataBasePath);
		String refData = testFilesFolder + fileSep + "1000gCeuChr20Mb6";

		GenotypeAligner.main("--debug", "--inputType", "PLINK_BED", "--input", studyDataBasePath, "--update-id", "--outputType", "PLINK_BED", "--output", tmpOutputFolder.getAbsolutePath() + fileSep + "test", "--refType", "VCF", "-ref", refData);

		System.out.println("Alignement complete now going to check using the real forward data");

		RandomAccessGenotypeData aligenedHapmap3Data = new BedBimFamGenotypeData(tmpOutputFolder.getAbsolutePath() + fileSep + "test");
		RandomAccessGenotypeData forwardHapmap3Data = new BedBimFamGenotypeData(testFilesFolder + fileSep + "hapmap3CeuChr20B37Mb6");

		//Check if the alles ar as expected acourding to real hapmap3 in forward strand
		int variantCounter = 0;
		for (GeneticVariant aligendVariant : aligenedHapmap3Data) {

			++variantCounter;

			GeneticVariant originalVariant = forwardHapmap3Data.getSnpVariantByPos(aligendVariant.getSequenceName(), aligendVariant.getStartPos());

			//Do not test these SNPs it is not on forward strand in hapmap3
			if (originalVariant.getPrimaryVariantId().equals("rs1047527") || originalVariant.getPrimaryVariantId().equals("rs2076553") || originalVariant.getPrimaryVariantId().equals("rs3761248")) {
				continue;
			}

			Iterator<Alleles> aligendVariantSampleAllelesIterator = aligendVariant.getSampleVariants().iterator();
			Iterator<Alleles> orignalVariantSampleAllelesIterator = originalVariant.getSampleVariants().iterator();

			while (aligendVariantSampleAllelesIterator.hasNext()) {
				Alleles aligned = aligendVariantSampleAllelesIterator.next();
				Alleles original = orignalVariantSampleAllelesIterator.next();
				assertEquals(aligned, original, "Inconsistant for SNP: " + aligendVariant.getPrimaryVariantId() + " " + aligned.getAllelesAsString() + " vs " + original.getAllelesAsString());
			}
			assertEquals(orignalVariantSampleAllelesIterator.hasNext(), false);

		}

		assertEquals(variantCounter, 3743);

		//Check if ID is updated based on 1000G
		assertEquals(aligenedHapmap3Data.getSnpVariantByPos("20", 809930).getPrimaryVariantId(), "rs78472400");

		//Check if number of samples is correct
		assertEquals(aligenedHapmap3Data.getSamples().size(), 165);
		assertEquals(aligenedHapmap3Data.getSnpVariantByPos("20", 809930).getSampleVariants().size(), 165);


	}

	/**
	 * Test of main method, of class GenotypeAligner.
	 */
	@Test
	public void testMain2() throws Exception {
		System.out.println("main");

		String studyDataBasePath = testFilesFolder + fileSep + "hapmap3CeuChr20B37Mb6RandomStrand";
		System.out.println(studyDataBasePath);
		String refData = testFilesFolder + fileSep + "1000gCeuChr20Mb6";

		GenotypeAligner.main("--debug", "--inputType", "PLINK_BED", "--input", studyDataBasePath, "--update-id", "--outputType", "PLINK_BED", "--output", tmpOutputFolder.getAbsolutePath() + fileSep + "test2", "--refType", "VCF", "-ref", refData, "--keep");

		System.out.println("Alignement complete now going to check using the real forward data");

		RandomAccessGenotypeData aligenedHapmap3Data = new BedBimFamGenotypeData(tmpOutputFolder.getAbsolutePath() + fileSep + "test2");
		RandomAccessGenotypeData forwardHapmap3Data = new BedBimFamGenotypeData(testFilesFolder + fileSep + "hapmap3CeuChr20B37Mb6");



		BufferedReader keepFileReader = new BufferedReader(new InputStreamReader(new FileInputStream(new File(testFilesFolder, "IncludedByKeep.txt")), FILE_ENCODING));

		HashSet<String> snpsKeptByKeepOption = new HashSet<String>();
		String snp;

		while ((snp = keepFileReader.readLine()) != null) {
			snpsKeptByKeepOption.add(snp);
		}

		//Check if the alles ar as expected acourding to real hapmap3 in forward strand
		int variantCounter = 0;
		for (GeneticVariant aligendVariant : aligenedHapmap3Data) {

			++variantCounter;

			GeneticVariant originalVariant = forwardHapmap3Data.getSnpVariantByPos(aligendVariant.getSequenceName(), aligendVariant.getStartPos());

			//Do not test these SNPs it is not on forward strand in hapmap3
			if (snpsKeptByKeepOption.contains(originalVariant.getPrimaryVariantId()) || originalVariant.getPrimaryVariantId().equals("rs1047527") || originalVariant.getPrimaryVariantId().equals("rs2076553") || originalVariant.getPrimaryVariantId().equals("rs3761248")) {
				continue;
			}

			Iterator<Alleles> aligendVariantSampleAllelesIterator = aligendVariant.getSampleVariants().iterator();
			Iterator<Alleles> orignalVariantSampleAllelesIterator = originalVariant.getSampleVariants().iterator();

			while (aligendVariantSampleAllelesIterator.hasNext()) {
				Alleles aligned = aligendVariantSampleAllelesIterator.next();
				Alleles original = orignalVariantSampleAllelesIterator.next();
				assertEquals(aligned, original, "Inconsistant for SNP: " + aligendVariant.getPrimaryVariantId() + " " + aligned.getAllelesAsString() + " vs " + original.getAllelesAsString());
			}
			assertEquals(orignalVariantSampleAllelesIterator.hasNext(), false);

		}

		assertEquals(variantCounter, 4084);

		//Check if ID is updated based on 1000G
		assertEquals(aligenedHapmap3Data.getSnpVariantByPos("20", 809930).getPrimaryVariantId(), "rs78472400");

		//Check if number of samples is correct
		assertEquals(aligenedHapmap3Data.getSamples().size(), 165);
		assertEquals(aligenedHapmap3Data.getSnpVariantByPos("20", 809930).getSampleVariants().size(), 165);


	}
	
	/**
	 * Test of main method, of class GenotypeAligner.
	 */
	@Test
	public void testMain3() throws Exception {
		System.out.println("main");

		String studyDataBasePath = testFilesFolder + fileSep + "hapmap3CeuChr20B37Mb6RandomStrand";
		System.out.println(studyDataBasePath);
		String refData = testFilesFolder + fileSep + "1000gCeuChr20Mb6";

		GenotypeAligner.main("--debug", "--inputType", "SHAPEIT2", "--input", studyDataBasePath, "--update-id", "--outputType", "PLINK_BED", "--output", tmpOutputFolder.getAbsolutePath() + fileSep + "test3", "--refType", "VCF", "-ref", refData, "--keep", "--forceChr", "20");

		System.out.println("Alignement complete now going to check using the real forward data");

		RandomAccessGenotypeData aligenedHapmap3Data = new BedBimFamGenotypeData(tmpOutputFolder.getAbsolutePath() + fileSep + "test3");
		RandomAccessGenotypeData forwardHapmap3Data = new BedBimFamGenotypeData(testFilesFolder + fileSep + "hapmap3CeuChr20B37Mb6");



		BufferedReader keepFileReader = new BufferedReader(new InputStreamReader(new FileInputStream(new File(testFilesFolder, "IncludedByKeep.txt")), FILE_ENCODING));

		HashSet<String> snpsKeptByKeepOption = new HashSet<String>();
		String snp;

		while ((snp = keepFileReader.readLine()) != null) {
			snpsKeptByKeepOption.add(snp);
		}

		//Check if the alles ar as expected acourding to real hapmap3 in forward strand
		int variantCounter = 0;
		for (GeneticVariant aligendVariant : aligenedHapmap3Data) {

			++variantCounter;

			GeneticVariant originalVariant = forwardHapmap3Data.getSnpVariantByPos(aligendVariant.getSequenceName(), aligendVariant.getStartPos());

			//Do not test these SNPs it is not on forward strand in hapmap3
			if (snpsKeptByKeepOption.contains(originalVariant.getPrimaryVariantId()) || originalVariant.getPrimaryVariantId().equals("rs1047527") || originalVariant.getPrimaryVariantId().equals("rs2076553") || originalVariant.getPrimaryVariantId().equals("rs3761248")) {
				continue;
			}

			Iterator<Alleles> aligendVariantSampleAllelesIterator = aligendVariant.getSampleVariants().iterator();
			Iterator<Alleles> orignalVariantSampleAllelesIterator = originalVariant.getSampleVariants().iterator();

			while (aligendVariantSampleAllelesIterator.hasNext()) {
				Alleles aligned = aligendVariantSampleAllelesIterator.next();
				Alleles original = orignalVariantSampleAllelesIterator.next();
				
				//Shapeits imputs sporadic missing genotypes. Ignore this
				if(original == MISSING_ALLELES){
					continue;
				}
				
				assertEquals(aligned, original, "Inconsistant for SNP: " + aligendVariant.getPrimaryVariantId() + " " + aligned.getAllelesAsString() + " vs " + original.getAllelesAsString());
			}
			assertEquals(orignalVariantSampleAllelesIterator.hasNext(), false);

		}

		assertEquals(variantCounter, 4084);

		//Check if ID is updated based on 1000G
		assertEquals(aligenedHapmap3Data.getSnpVariantByPos("20", 809930).getPrimaryVariantId(), "rs78472400");

		//Check if number of samples is correct
		assertEquals(aligenedHapmap3Data.getSamples().size(), 165);
		assertEquals(aligenedHapmap3Data.getSnpVariantByPos("20", 809930).getSampleVariants().size(), 165);


	}
}
