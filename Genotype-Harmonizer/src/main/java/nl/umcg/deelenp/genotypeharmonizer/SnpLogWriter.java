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

	public static enum Actions{
		
		EXCLUDED("Excluded"),
		SWAPPED("Swapped"),
		MAINTAINED("Maintained");
		
		private String actionString;
		
		private Actions(String actionString) {
			this.actionString = actionString;
		}
		
		public String getActionString() {
			return actionString;
		}
			
	}
	
	private final BufferedWriter snpLogWriter;

	public SnpLogWriter(File snpLogFile) throws IOException {
		 snpLogWriter = new BufferedWriter(new FileWriter(snpLogFile));
		 snpLogWriter.append("chr\tpos\tid\talleles\taction\tmessage\n");
	}
	
	public void addToLog(GeneticVariant variant, Actions action, String message) throws IOException{
		snpLogWriter.append(variant.getSequenceName());
		snpLogWriter.append('\t');
		snpLogWriter.append(String.valueOf(variant.getStartPos()));
		snpLogWriter.append('\t');
		snpLogWriter.append(variant.getPrimaryVariantId());
		snpLogWriter.append('\t');
		snpLogWriter.append(variant.getVariantAlleles().toString());
		snpLogWriter.append('\t');
		snpLogWriter.append(action.getActionString());
		snpLogWriter.append('\t');
		snpLogWriter.append(message);
		snpLogWriter.append('\n');
	}

	@Override
	public void close() throws IOException {
		snpLogWriter.close();
	}
	
	
	
}
