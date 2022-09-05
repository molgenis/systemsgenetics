package mbqtl.vcf;

import mbqtl.vcf.filter.genotypefilters.VCFGenotypeFilter;
import umcg.genetica.io.text.TextFile;

import java.io.IOException;
import java.util.ArrayList;
import java.util.HashSet;
import java.util.Iterator;

/**
 * Created by hwestra on 5/11/15.
 */
public class VCFGenotypeData implements Iterator<VCFVariant> {

	private ArrayList<VCFGenotypeFilter> genotypeFilters;
	TextFile tf = null;
	String next = null;

	public VCFGenotypeData(TextFile tf, HashSet<String> excludeTheseSamples, ArrayList<VCFGenotypeFilter> genotypeFilters) throws IOException {
		this.tf = tf;
		this.genotypeFilters = genotypeFilters;
		tf.close();
		tf.open();

		init();
	}

	public VCFGenotypeData(TextFile tf, HashSet<String> excludeTheseSamples) throws IOException {
		this.tf = tf;
		tf.close();
		tf.open();

		init();
	}

	public VCFGenotypeData(TextFile tf) throws IOException {
		this.tf = tf;
		tf.close();
		tf.open();

		init();
	}

	public VCFGenotypeData(String vcffile) throws IOException {
		tf = new TextFile(vcffile, TextFile.R);
		init();

	}

	public void open() throws IOException {
		tf.open();
		init();
	}

	private void init() throws IOException {
		String ln = tf.readLine();
		while (ln.startsWith("#")) {
			if (ln.startsWith("#CHROM")) {
				parseheader(ln);
			}
			ln = tf.readLine();
		}

		if (ln != null) {
			next = ln; // new VCFVariant(ln);
		} else {
			next = null;
		}
	}

	ArrayList<String> samples = new ArrayList<String>();

	public ArrayList<String> getSamples() {
		return samples;
	}

	private void parseheader(String ln) {
		String[] elems = ln.split("\t");
		// #CHROM  POS     ID      REF     ALT     QUAL    FILTER  INFO    FORMAT

		for (int e = 9; e < elems.length; e++) {
			String sample = elems[e];

			samples.add(sample);


		}


	}

	@Override
	public boolean hasNext() {
		return next != null;
	}

	@Override
	public VCFVariant next() {
		VCFVariant current = new VCFVariant(next, genotypeFilters, true);
		try {
			String ln = tf.readLine();
			if (ln != null) {
				next = ln; // new VCFVariant(ln, genotypeFilters, true);
			} else {
				next = null;
			}

		} catch (IOException ex) {
			throw new RuntimeException(ex);
		}
//		return current;
		return null;
	}

	@Override
	public void remove() {

	}

	public void close() throws IOException {
		this.tf.close();
	}

//	public VCFVariant nextLoadGenotypesOnly() {
//		VCFVariant current = next;
//		try {
//
//			String ln = tf.readLine();
//			if (ln != null) {
//				next = new VCFVariant(ln, VCFVariant.PARSE.GENOTYPES);
//			} else {
//				next = null;
//			}
//
//		} catch (IOException ex) {
//			throw new RuntimeException(ex);
//		}
//		return current;
//	}

	public VCFVariant nextLoadHeader() {
		VCFVariant current = new VCFVariant(next, VCFVariant.PARSE.HEADER);
		try {

			String ln = tf.readLine();
			if (ln != null) {
				next = ln; // new VCFVariant(ln, VCFVariant.PARSE.HEADER);
			} else {
				next = null;
			}

		} catch (IOException ex) {
			throw new RuntimeException(ex);
		}
//		return current;
		return null;
	}


}
