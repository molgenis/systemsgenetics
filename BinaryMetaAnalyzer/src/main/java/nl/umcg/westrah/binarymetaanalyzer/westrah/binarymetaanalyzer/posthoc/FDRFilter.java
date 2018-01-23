package nl.umcg.westrah.binarymetaanalyzer.westrah.binarymetaanalyzer.posthoc;

import umcg.genetica.io.text.TextFile;
import umcg.genetica.text.Strings;

import java.io.IOException;

public class FDRFilter {
	
	public static void main(String[] args) {
		FDRFilter f = new FDRFilter();
		String[] filesin = new String[]{
				"C:\\Sync\\OneDrive\\Postdoc2\\2017-11-eQTLMeta\\Heterogeneity\\eqtl\\stegle\\trans_eQTL_replication_CD4_Stegle_20171227.txt.gz",
				"C:\\Sync\\OneDrive\\Postdoc2\\2017-11-eQTLMeta\\Heterogeneity\\eqtl\\stegle\\trans_eQTL_replication_CD8_Stegle_20171227.txt.gz",
				"C:\\Sync\\OneDrive\\Postdoc2\\2017-11-eQTLMeta\\Heterogeneity\\eqtl\\stegle\\trans_eQTL_replication_CD14_Stegle_20171227.txt.gz",
				"C:\\Sync\\OneDrive\\Postdoc2\\2017-11-eQTLMeta\\Heterogeneity\\eqtl\\stegle\\trans_eQTL_replication_CD15_Stegle_20171227.txt.gz",
				"C:\\Sync\\OneDrive\\Postdoc2\\2017-11-eQTLMeta\\Heterogeneity\\eqtl\\stegle\\trans_eQTL_replication_CD19_Stegle_20171227.txt.gz",
				"C:\\Sync\\OneDrive\\Postdoc2\\2017-11-eQTLMeta\\Heterogeneity\\eqtl\\stegle\\trans_eQTL_replication_Mac_Stegle_20171227.txt.gz",
				"C:\\Sync\\OneDrive\\Postdoc2\\2017-11-eQTLMeta\\Heterogeneity\\eqtl\\stegle\\trans_eQTL_replication_PLT_Stegle_20171227.txt.gz",
				"C:\\Sync\\OneDrive\\Postdoc2\\2017-11-eQTLMeta\\Heterogeneity\\eqtl\\2017-12-20-LCLMeta\\eQTLs-FDR0.05.txt"
			
		};
		
		String[] filesout = new String[]{
				"C:\\Sync\\OneDrive\\Postdoc2\\2017-11-eQTLMeta\\Heterogeneity\\eqtl\\stegle\\trans_eQTL_replication_CD4_Stegle_20171227-fdrfilter.txt.gz",
				"C:\\Sync\\OneDrive\\Postdoc2\\2017-11-eQTLMeta\\Heterogeneity\\eqtl\\stegle\\trans_eQTL_replication_CD8_Stegle_20171227-fdrfilter.txt.gz",
				"C:\\Sync\\OneDrive\\Postdoc2\\2017-11-eQTLMeta\\Heterogeneity\\eqtl\\stegle\\trans_eQTL_replication_CD14_Stegle_20171227-fdrfilter.txt.gz",
				"C:\\Sync\\OneDrive\\Postdoc2\\2017-11-eQTLMeta\\Heterogeneity\\eqtl\\stegle\\trans_eQTL_replication_CD15_Stegle_20171227-fdrfilter.txt.gz",
				"C:\\Sync\\OneDrive\\Postdoc2\\2017-11-eQTLMeta\\Heterogeneity\\eqtl\\stegle\\trans_eQTL_replication_CD19_Stegle_20171227-fdrfilter.txt.gz",
				"C:\\Sync\\OneDrive\\Postdoc2\\2017-11-eQTLMeta\\Heterogeneity\\eqtl\\stegle\\trans_eQTL_replication_Mac_Stegle_20171227-fdrfilter.txt.gz",
				"C:\\Sync\\OneDrive\\Postdoc2\\2017-11-eQTLMeta\\Heterogeneity\\eqtl\\stegle\\trans_eQTL_replication_PLT_Stegle_20171227-fdrfilter.txt.gz",
				"C:\\Sync\\OneDrive\\Postdoc2\\2017-11-eQTLMeta\\Heterogeneity\\eqtl\\2017-12-20-LCLMeta\\eQTLs-FDR0.05-fdrfilter.txt"
			
		};
		
		for (int q = 0; q < filesin.length; q++) {
			try {
				f.filter(filesin[q], filesout[q], 0.05);
			} catch (IOException e) {
				e.printStackTrace();
			}
		}
	}
	
	
	public void filter(String inf, String outf, double threshold) throws IOException {
		
		TextFile in = new TextFile(inf, TextFile.R);
		TextFile out = new TextFile(outf, TextFile.W);
		out.writeln(in.readLine());
		String[] elems = in.readLineElems(TextFile.tab);
		while (elems != null) {
			String fdr = elems[elems.length - 1];
			Double fdrd = Double.parseDouble(fdr);
			if (fdrd < threshold) {
				out.writeln(Strings.concat(elems, Strings.tab));
			}
			elems = in.readLineElems(TextFile.tab);
		}
		in.close();
		out.close();
	}
	
}


