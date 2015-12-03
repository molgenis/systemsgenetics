package nl.systemsgenetics.eqtlinteractionanalyser.eqtlinteractionanalyser;

import au.com.bytecode.opencsv.CSVWriter;
import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.util.ArrayList;

/**
 *
 * @author Patrick Deelen
 */
public class SkippedInteractionWriter {
	
	private final CSVWriter writer;
	private final String[] row = new String[5];
	private int c;
	private StringBuilder tmp;

	public SkippedInteractionWriter(File skippedInteractionsFile) throws IOException {
		writer = new CSVWriter(new FileWriter(skippedInteractionsFile), '\t', CSVWriter.NO_QUOTE_CHARACTER);
		
		c = 0;
		row[c++] = "Covariate";
		row[c++] = "CountSingular";
		row[c++] = "CountSharedQtl";
		row[c++] = "SingularQtls";
		row[c++] = "SharedQtls";
		
		writer.writeNext(row);
	}
	
	public void close() throws IOException{
		writer.close();
	}

	synchronized void add(SkippedInteractionTracker skipped){
		
		ArrayList<String> singular = skipped.getSkipped().get(SkippedInteractionTracker.Reason.SINGULAR);
		ArrayList<String> sharedQtl = skipped.getSkipped().get(SkippedInteractionTracker.Reason.SHARED_QTL);
		
		c = 0;
		row[c++] = skipped.getCovariate();
		row[c++] = String.valueOf(singular.size());
		row[c++] = String.valueOf(sharedQtl.size());
		
		tmp = new StringBuilder();
		for(String qtl : singular){
			tmp.append(qtl);
			tmp.append(';');
		}
		row[c++] = tmp.toString();
		
		tmp = new StringBuilder();
		for(String qtl : sharedQtl){
			tmp.append(qtl);
			tmp.append(';');
		}
		row[c++] = tmp.toString();
		
		writer.writeNext(row);
		
	}
	
}
