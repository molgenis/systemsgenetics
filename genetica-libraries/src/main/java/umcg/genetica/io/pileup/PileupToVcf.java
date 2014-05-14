package umcg.genetica.io.pileup;

import java.io.BufferedWriter;
import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import org.apache.commons.cli.CommandLine;
import org.apache.commons.cli.HelpFormatter;
import org.apache.commons.cli.OptionBuilder;
import org.apache.commons.cli.Options;
import org.apache.commons.cli.ParseException;
import org.apache.commons.cli.PosixParser;
import org.molgenis.genotype.Allele;
import org.molgenis.genotype.RandomAccessGenotypeData;
import org.molgenis.genotype.RandomAccessGenotypeDataReaderFormats;
import org.molgenis.genotype.variant.GeneticVariant;

/**
 *
 * @author Patrick Deelen
 */
public class PileupToVcf {

	private static final Options OPTIONS;
	
	static {
		OPTIONS = new Options();
	
		OptionBuilder.withArgName("path");
		OptionBuilder.hasArgs();
		OptionBuilder.isRequired();
		OptionBuilder.withDescription("Pileup file max 1 sample");
		OPTIONS.addOption(OptionBuilder.create('p'));
		
		OptionBuilder.withArgName("string");
		OptionBuilder.hasArgs();
		OptionBuilder.isRequired();
		OptionBuilder.withDescription("Sample name");
		OPTIONS.addOption(OptionBuilder.create('s'));
		
		OptionBuilder.withArgName("path");
		OptionBuilder.hasArgs();
		OptionBuilder.isRequired();
		OptionBuilder.withDescription("Output VCF");
		OPTIONS.addOption(OptionBuilder.create('v'));
		
		OptionBuilder.withArgName("path");
		OptionBuilder.hasArgs();
		OptionBuilder.isRequired();
		OptionBuilder.withDescription("Reference vcf. Output only contains SNPs present in this file with same ref and alt alleles");
		OPTIONS.addOption(OptionBuilder.create('r'));
	
	}
	
	public static void main(String[] args) throws ParseException, IOException {
		
		if (args.length == 0) {
			new HelpFormatter().printHelp(" ", OPTIONS);
		}
			
		final CommandLine commandLine = new PosixParser().parse(OPTIONS, args, true);

		
		RandomAccessGenotypeData referenceVcf = RandomAccessGenotypeDataReaderFormats.VCF.createGenotypeData(commandLine.getOptionValue('r'));
		
		PileupFile pileupFile = new PileupFile(commandLine.getOptionValue('p'));
		
		File outputVcf = new File(commandLine.getOptionValue('v'));
		
		String sample = commandLine.getOptionValue('s');
		
		convertPileupToVcf(pileupFile, sample, outputVcf, referenceVcf);
		
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
    public static void convertPileupToVcf(final PileupFile pileupFile, final String sample, final File outputVcf, final RandomAccessGenotypeData referenceGenotypes) throws IOException {
        
		final BufferedWriter vcfWriter = new BufferedWriter(new FileWriter(outputVcf));
			
		vcfWriter.append("##fileformat=VCFv4.2\n");
		vcfWriter.append("##source=ConvertPileupToVcf\n");
		vcfWriter.append("##FORMAT=<ID=AD,Number=R,Type=Integer,Description=\"Allelic depths for the ref and alt alleles in the order listed\">\n");
		vcfWriter.append("##FORMAT=<ID=RQ,Number=R,Type=Float,Description=\"Average read quality of alleles\">\n");
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
				vcfWriter.append("\t.\tPASS\t.\tAD:RQ\t");
				vcfWriter.append(String.valueOf(refCount));
				vcfWriter.append(',');
				vcfWriter.append(String.valueOf(altCount));
				vcfWriter.append(':');
				vcfWriter.append(String.valueOf(pileupEntry.getAlleleAverageQuality(refAllele)));
				vcfWriter.append(',');
				vcfWriter.append(String.valueOf(pileupEntry.getAlleleAverageQuality(altAllele)));
				vcfWriter.append('\n');
				
			}
		}
		
		vcfWriter.close();
	
    }

}
