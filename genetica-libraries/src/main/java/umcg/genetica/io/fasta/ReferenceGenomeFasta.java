package umcg.genetica.io.fasta;

import cern.colt.matrix.tdouble.impl.DenseDoubleMatrix1D;
import com.lowagie.text.LargeElement;
import java.io.BufferedReader;
import java.io.File;
import java.io.FileReader;
import java.io.IOException;
import java.util.HashMap;
import java.util.regex.Matcher;
import java.util.regex.Pattern;

/**
 *
 * @author Patrick Deelen
 */
public class ReferenceGenomeFasta {

	private final HashMap<String, LargeByteArray> chromosomes;
	private static final Pattern FASTA_HEADER_PATTERN = Pattern.compile("^\\>(\\S+).*(\\d+):\\d+$");

	public ReferenceGenomeFasta(File fastaFile) throws IOException, Exception {

		chromosomes = new HashMap<String, LargeByteArray>(64);

		BufferedReader fasteReader = new BufferedReader(new FileReader(fastaFile));

		String line;
		DenseDoubleMatrix1D currentSeq = null;

		while ((line = fasteReader.readLine()) != null) {

			if (line.charAt(0) == '>') {

				Matcher headerMatcher = FASTA_HEADER_PATTERN.matcher(line);

				if (!headerMatcher.matches()) {
					throw new Exception("Error parsing refernece genome fasta header: " + line);
				}

				String chr = headerMatcher.group(1);
				long lenght = Long.parseLong(headerMatcher.group(2));
				
				//currentSeq = new LargeByteArray
				
			}

		}

	}
}
