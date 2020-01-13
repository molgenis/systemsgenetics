package umcg.genetica.features;

public enum FeatureType {
	
	STARTCODON,
	STOPCODON,
	CDS,
	EXON,
	UTR,
	OTHER;
	
	public static FeatureType parse(String str) {
		String lc = str.toLowerCase();
		if (lc.equals("start_codon")) {
			return STARTCODON;
		} else if (lc.equals("stop_codon")) {
			return STOPCODON;
		} else if (lc.equals("cds")) {
			return CDS;
		} else if (lc.equals("exon")) {
			return EXON;
		} else if (lc.equals("utr")) {
			return UTR;
		} else {
			return OTHER;
		}
	}
	
	
}
