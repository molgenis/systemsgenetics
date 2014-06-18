package nl.umcg.deelenp.genotypeharmonizer;

import java.io.BufferedWriter;
import java.io.Closeable;
import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import org.molgenis.genotype.variant.GeneticVariant;

/**
 *
 * @author Patrick Deelen
 */
public class SnpLogWriter implements Closeable{

	private final BufferedWriter snpLogWriter;

	public SnpLogWriter(File snpLogFile) throws IOException {
		 snpLogWriter = new BufferedWriter(new FileWriter(snpLogFile));
		 snpLogWriter.append("chr\tpos\tid\talleles\taction\tmessage\n");
	}
	
	public void addToLog(GeneticVariant variant,String action, String message) throws IOException{
		snpLogWriter.append(variant.getSequenceName());
		snpLogWriter.append('\t');
		snpLogWriter.append(String.valueOf(variant.getStartPos()));
		snpLogWriter.append('\t');
		snpLogWriter.append(variant.getPrimaryVariantId());
		snpLogWriter.append('\t');
		snpLogWriter.append(variant.getVariantAlleles().toString());
		snpLogWriter.append('\t');
		snpLogWriter.append(action);
		snpLogWriter.append('\t');
		snpLogWriter.append(message);
		snpLogWriter.append('\n');
	}

	@Override
	public void close() throws IOException {
		snpLogWriter.close();
	}
	
	
	
}
