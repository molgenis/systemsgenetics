package umcg.genetica.io.pileup;

import gnu.trove.map.TIntObjectMap;
import gnu.trove.map.hash.TIntObjectHashMap;
import java.io.BufferedWriter;
import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.text.DateFormat;
import java.text.SimpleDateFormat;
import java.util.Date;
import java.util.HashMap;
import org.apache.commons.cli.CommandLine;
import org.apache.commons.cli.HelpFormatter;
import org.apache.commons.cli.OptionBuilder;
import org.apache.commons.cli.Options;
import org.apache.commons.cli.ParseException;
import org.apache.commons.cli.PosixParser;
import org.molgenis.genotype.Allele;
import org.molgenis.genotype.Alleles;
import org.molgenis.genotype.RandomAccessGenotypeData;
import org.molgenis.genotype.RandomAccessGenotypeDataReaderFormats;
import org.molgenis.genotype.variant.GeneticVariant;

/**
 *
 * @author Patrick Deelen
 */
public class PileupToVcf {

	private static final String VERSION = "1.3";
	private static final String HEADER =
			"  /---------------------------------------\\\n"
			+ "  |             Pileup to VCF             |\n"
			+ "  |                                       |\n"
			+ "  |             Patrick Deelen            |\n"
			+ "  |        patrickdeelen@gmail.com        |\n"
			+ "  |                                       |\n"
			+ "  |           Dasha Zhernakova,           |\n"
			+ "  |      Lude Franke, Morris Swertz       |\n"
			+ "  |                                       |\n"
			+ "  |     Genomics Coordination Center      |\n"
			+ "  |        Department of Genetics         |\n"
			+ "  |  University Medical Center Groningen  |\n"
			+ "  \\---------------------------------------/";
	private static final Options OPTIONS;
	private static final DateFormat DATE_TIME_FORMAT = new SimpleDateFormat("yyyy-MM-dd HH:mm:ss");
	private static final Date currentDataTime = new Date();
	
	static {
		OPTIONS = new Options();
	
		OptionBuilder.withArgName("path");
		OptionBuilder.hasArg();
		OptionBuilder.isRequired();
		OptionBuilder.withDescription("Pileup file max 1 sample");
		OPTIONS.addOption(OptionBuilder.create('p'));
		
		OptionBuilder.withArgName("string");
		OptionBuilder.hasArg();
		OptionBuilder.isRequired();
		OptionBuilder.withDescription("Sample name");
		OPTIONS.addOption(OptionBuilder.create('s'));
		
		OptionBuilder.withArgName("path");
		OptionBuilder.hasArg();
		OptionBuilder.isRequired();
		OptionBuilder.withDescription("Output VCF");
		OPTIONS.addOption(OptionBuilder.create('v'));
		
		OptionBuilder.withArgName("path");
		OptionBuilder.hasArgs();
		OptionBuilder.isRequired();
		OptionBuilder.withDescription("Reference vcf. Output only contains SNPs present in this file with same ref and alt alleles");
		OPTIONS.addOption(OptionBuilder.create('r'));
		
		OptionBuilder.withArgName("int");
		OptionBuilder.hasArgs();
		OptionBuilder.withDescription("Minimum base quality to include read in ref / alt count");
		OPTIONS.addOption(OptionBuilder.create('b'));
	
	}
	
	public static void main(String[] args) throws ParseException, IOException {
		
		System.out.println(HEADER);
		System.out.println();
		System.out.println("          --- Version: " + VERSION + " ---");
		System.out.println();
		System.out.println("More information: http://molgenis.org/systemsgenetics");
		System.out.println();

		System.out.println("Current date and time: " + DATE_TIME_FORMAT.format(currentDataTime));
		System.out.println();

		System.out.flush(); //flush to make sure header is before errors
		try {
			Thread.sleep(25); //Allows flush to complete
		} catch (InterruptedException ex) {
		}
		
		if (args.length == 0) {
			new HelpFormatter().printHelp(" ", OPTIONS);
		}
			
		final CommandLine commandLine = new PosixParser().parse(OPTIONS, args, true);

		final int minBaseQuality = commandLine.hasOption('b') ? Integer.parseInt(commandLine.getOptionValue('b')) : 0;
		final RandomAccessGenotypeData referenceVcf = RandomAccessGenotypeDataReaderFormats.VCF.createGenotypeData(commandLine.getOptionValue('r'));
		final PileupFile pileupFile = new PileupFile(commandLine.getOptionValue('p'), minBaseQuality);
		final File outputVcf = new File(commandLine.getOptionValue('v'));
		final String sample = commandLine.getOptionValue('s');
				
		System.out.println("Ref VCF: " + commandLine.getOptionValue('r'));
		System.out.println("Pileup file: " + pileupFile.getAbsolutePath());
		System.out.println("Output VCF: " + outputVcf.getAbsolutePath());
		System.out.println("Sample name: " + sample);
		System.out.println("Min base quality: " + minBaseQuality);
		
		convertPileupToVcf2(pileupFile, sample, outputVcf, referenceVcf, minBaseQuality);
		
		System.out.println("Conversion complete");
		
	}
	
		
   /**
	* Parse pileup file and create vcf files with only SNPs that are in referenceGenotype data. VCF only contains the counts for the 2 alleles in the reference genotype data. Only parses the first samples in the pileup file
	* 
	* @param pileupFile
	* @param sample
	* @param outputVcf
	* @param referenceGenotypes
	* @throws IOException 
	*/
    public static void convertPileupToVcf(final PileupFile pileupFile, final String sample, final File outputVcf, final RandomAccessGenotypeData referenceGenotypes, final int minBaseQuality) throws IOException {
        
		final BufferedWriter vcfWriter = new BufferedWriter(new FileWriter(outputVcf));
			
		vcfWriter.append("##fileformat=VCFv4.2\n");
		vcfWriter.append("##source=ConvertPileupToVcf version: " + VERSION + "\n");
		vcfWriter.append("##pileupFile=" + pileupFile.getAbsolutePath() + "\n");
		vcfWriter.append("##minBaseQuality=" + minBaseQuality + "\n");
		vcfWriter.append("##FORMAT=<ID=AD,Number=R,Type=Integer,Description=\"Allelic depths for the ref and alt alleles in the order listed\">\n");
		vcfWriter.append("##FORMAT=<ID=RQ,Number=R,Type=Float,Description=\"Average read quality of alleles\">\n");
		vcfWriter.append("##FORMAT=<ID=TD,Number=1,Type=Integer,Description=\"Total read depth including non ref and alt reads and read with low quality bases\">\n");
		vcfWriter.append("##FORMAT=<ID=ADQ,Number=1,Type=Float,Description=\"Minimum base quality to include in ref and alt counts\">\n");
		vcfWriter.append("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t" + sample + "\n");
		
		for(PileupEntry pileupEntry : pileupFile){
			GeneticVariant refVariant = referenceGenotypes.getSnpVariantByPos(pileupEntry.getChr(), pileupEntry.getPos());
			if(refVariant != null){
				
				Allele refAllele = refVariant.getVariantAlleles().get(0);
				Allele altAllele = refVariant.getVariantAlleles().get(1);
				
				int refCount = pileupEntry.getAlleleCount(refAllele);
				int altCount = pileupEntry.getAlleleCount(altAllele);
				
				if(refCount == 0 && altCount == 0){
					continue;
				}
			
				vcfWriter.append(refVariant.getSequenceName());
				vcfWriter.append('\t');
				vcfWriter.append(String.valueOf(refVariant.getStartPos()));
				vcfWriter.append('\t');
				vcfWriter.append(refVariant.getPrimaryVariantId());
				vcfWriter.append('\t');
				vcfWriter.append(refAllele.getAlleleAsString());
				vcfWriter.append('\t');
				vcfWriter.append(altAllele.getAlleleAsString());
				vcfWriter.append("\t.\tPASS\t.\tAD:RQ:TD:ADQ\t");
				vcfWriter.append(String.valueOf(refCount));
				vcfWriter.append(',');
				vcfWriter.append(String.valueOf(altCount));
				vcfWriter.append(':');
				vcfWriter.append(String.valueOf(pileupEntry.getAlleleAverageQuality(refAllele)));
				vcfWriter.append(',');
				vcfWriter.append(String.valueOf(pileupEntry.getAlleleAverageQuality(altAllele)));
				vcfWriter.append(':');
				vcfWriter.append(String.valueOf(pileupEntry.getReadDepth()));
				vcfWriter.append(':');
				vcfWriter.append(String.valueOf(pileupEntry.getMinimumBaseQuality()));
				vcfWriter.append('\n');
				
			}
		}
		
		vcfWriter.close();
	
    }

	/**
	* Parse pileup file and create vcf files with only SNPs that are in referenceGenotype data. VCF only contains the counts for the 2 alleles in the reference genotype data. Only parses the first samples in the pileup file
	* 
	* @param pileupFile
	* @param sample
	* @param outputVcf
	* @param referenceGenotypes
	* @throws IOException 
	*/
    public static void convertPileupToVcf2(final PileupFile pileupFile, final String sample, final File outputVcf, final RandomAccessGenotypeData referenceGenotypes, final int minBaseQuality) throws IOException {
        
		HashMap<String, TIntObjectMap<Alleles>> variantsAlleles = new HashMap<String, TIntObjectMap<Alleles>>();
		
		for(GeneticVariant variant : referenceGenotypes){
			
			TIntObjectMap<Alleles> chrVariantAlleles = variantsAlleles.get(variant.getSequenceName());
			if(chrVariantAlleles == null){
				chrVariantAlleles = new TIntObjectHashMap<Alleles>();
				variantsAlleles.put(variant.getSequenceName(), chrVariantAlleles);
			}
			chrVariantAlleles.put(variant.getStartPos(), variant.getVariantAlleles());
			
		}
		
		
		final BufferedWriter vcfWriter = new BufferedWriter(new FileWriter(outputVcf));
			
		vcfWriter.append("##fileformat=VCFv4.2\n");
		vcfWriter.append("##source=ConvertPileupToVcf version: " + VERSION + "\n");
		vcfWriter.append("##pileupFile=" + pileupFile.getAbsolutePath() + "\n");
		vcfWriter.append("##minBaseQuality=" + minBaseQuality + "\n");
		vcfWriter.append("##FORMAT=<ID=AD,Number=R,Type=Integer,Description=\"Allelic depths for the ref and alt alleles in the order listed\">\n");
		vcfWriter.append("##FORMAT=<ID=RQ,Number=R,Type=Float,Description=\"Average read quality of alleles\">\n");
		vcfWriter.append("##FORMAT=<ID=TD,Number=1,Type=Integer,Description=\"Total read depth including non ref and alt reads and read with low quality bases\">\n");
		vcfWriter.append("##FORMAT=<ID=ADQ,Number=1,Type=Float,Description=\"Minimum base quality to include in ref and alt counts\">\n");
		vcfWriter.append("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t" + sample + "\n");
		
		for(PileupEntry pileupEntry : pileupFile){
			
			Alleles variantAlleles = null;
			TIntObjectMap<Alleles> chrVariantAlleles = variantsAlleles.get(pileupEntry.getChr());
			if(chrVariantAlleles != null){
				variantAlleles = chrVariantAlleles.get(pileupEntry.getPos());
			}
			
			
			if(variantAlleles != null){
				
				Allele refAllele = variantAlleles.get(0);
				Allele altAllele = variantAlleles.get(1);
				
				int refCount = pileupEntry.getAlleleCount(refAllele);
				int altCount = pileupEntry.getAlleleCount(altAllele);
				
				if(refCount == 0 && altCount == 0){
					continue;
				}
			
				vcfWriter.append(pileupEntry.getChr());
				vcfWriter.append('\t');
				vcfWriter.append(String.valueOf(pileupEntry.getPos()));
				vcfWriter.append('\t');
				vcfWriter.append('.');
				vcfWriter.append('\t');
				vcfWriter.append(refAllele.getAlleleAsString());
				vcfWriter.append('\t');
				vcfWriter.append(altAllele.getAlleleAsString());
				vcfWriter.append("\t.\tPASS\t.\tAD:RQ:TD:ADQ\t");
				vcfWriter.append(String.valueOf(refCount));
				vcfWriter.append(',');
				vcfWriter.append(String.valueOf(altCount));
				vcfWriter.append(':');
				vcfWriter.append(String.valueOf(pileupEntry.getAlleleAverageQuality(refAllele)));
				vcfWriter.append(',');
				vcfWriter.append(String.valueOf(pileupEntry.getAlleleAverageQuality(altAllele)));
				vcfWriter.append(':');
				vcfWriter.append(String.valueOf(pileupEntry.getReadDepth()));
				vcfWriter.append(':');
				vcfWriter.append(String.valueOf(pileupEntry.getMinimumBaseQuality()));
				vcfWriter.append('\n');
				
			}
		}
		
		vcfWriter.close();
	
    }

	
}
