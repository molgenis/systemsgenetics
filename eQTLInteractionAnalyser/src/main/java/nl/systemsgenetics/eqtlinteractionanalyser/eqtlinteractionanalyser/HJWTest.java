package nl.systemsgenetics.eqtlinteractionanalyser.eqtlinteractionanalyser;

public class HJWTest {

	public static void main(String[] args) {
		args = new String[]{"-i", "D:\\interactiontest\\databinary2\\",
				"-o", "D:\\interactiontest\\output2\\",
				"-e", "D:\\interactiontest\\data\\eqtls-filtered.txt",
				"-is", "D:\\interactiontest\\databinary\\samples.txt",
				"-n", "1",
				"-pc", "5",
				"-nt", "8"};

//		args = new String[]{"-i", "D:\\interactiontest\\databinary2\\",
//				"-o", "D:\\interactiontest\\output2\\",
//				"-e", "D:\\interactiontest\\data\\eqtls-filtered.txt",
//				"-is", "D:\\interactiontest\\databinary\\samples.txt",
//				"-nt","8",
//				"-c", "PCT_MRNA_BASES", "MEDIAN_3PRIME_BIAS", "PCT_UTR_BASES", "PCT_CODING_BASES", "PCT_INTERGENIC_BASES","PCT_INTRONIC_BASES", "PCT_USABLE_BASES"
//		};

		try {
			EQTLInteractionAnalyser.main(args);
		} catch (Exception ex) {
			throw new RuntimeException(ex);
		}
	}

}
