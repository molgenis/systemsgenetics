/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package nl.systemsgenetics.eqtlpermutationtranscriptionfactoranalysis;
import java.io.IOException;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.HashSet;
import java.util.regex.Pattern;
//import org.bouncycastle.util.Arrays;
import umcg.genetica.genomicboundaries.GenomicBoundaries;
import umcg.genetica.graphics.ViolinBoxPlot;
import umcg.genetica.io.text.TextFile;
import umcg.genetica.io.ucsc.UCSCDataObject;
import umcg.genetica.io.ucsc.WigFile;
import umcg.genetica.math.stats.WilcoxonMannWhitney;
import com.google.common.primitives.Doubles;
import umcg.genetica.io.gwascatalog.GWASCatalog;
import umcg.genetica.math.stats.ZScores;

/**
 *
 * @author Matthieu
 */
public class Testing {
	private static final Pattern TAB_PATTERN = Pattern.compile("\t");
	public static void main(String[] args)throws IOException{
		/*
		WigFile wf = new WigFile("C:\\Users\\Matthieu\\Documents\\Afstudeerstage\\Pilot\\3.histones\\Gm12864Ctcf.wig", false);
		String fileLine;
		
		
		long totalN = wf.size();
		long n = 0;
		while(n < totalN){
			UCSCDataObject ucscdo = wf.parseLn();
			System.out.println( ucscdo.getChr() + "\t" + ucscdo.getPositionStart() + "\t" + ucscdo.getPositionEnd() + "\t" + ucscdo.getValue() );
		}
		*/
		
		
//		System.out.println("Testing the Genetica Library Wilcox test.");
//		double[] A = new double[]{35,23,58,74,70,67,34,76,34,73,40,76,34,80,7};
//		double[] B = new double[]{23,48,48,90,68,34,90,86,24,86,92,63,78,94,65,89,34,67,89,34,6};
//		System.out.println("Made two arrays with values: {A: " + A.length + "} :-: {B:" + B.length + "}");
//		WilcoxonMannWhitney wcmw = new WilcoxonMannWhitney();
//		double testPvalue = wcmw.returnWilcoxonMannWhitneyPValue(A, B);
//		double wilcoxAuc = wcmw.getAUC();
//		System.out.println("PValue is: " + testPvalue + " and the UAC value is: " + wilcoxAuc);
//		
//		
//		System.out.println("Iets met boxplots...");
//		ViolinBoxPlot vbp = new ViolinBoxPlot();
//		double[][][] setA = new double[1][1][15];
//		setA[0][0] = A;
//		System.out.println(setA[0][0]);
		
		
		
//		double[] A = new double[]{35,23,58,74,70,67,34,76,34,73,40,76,34,80,7};
//		double[] B = new double[]{23,48,48,90,68,34,90,86,24,86,92,63,78,94,65,89,34,67,89,34,6};
//		double[][][] setA = new double[1][2][0];
//		setA[0][0] = A;
//		setA[0][1] = B;
//		
//		String[][] plotLabels = new String[][]{{"Dingen", "Dingen2"}};
//		String[] dataset = new String[]{"Aap"}; 
//		
//		
//		//Test to see if the values are added.
//		for(int n=0;n<15;n++){
//			System.out.println(setA[0][0][n]);
//		}
//		
//		
//		ViolinBoxPlot vbp = new ViolinBoxPlot();
//		System.out.println("About to draw something;");
//		vbp.draw(setA, dataset, plotLabels, "Iets", ViolinBoxPlot.Output.PDF, "C:\\Users\\Matthieu\\Documents\\Afstudeerstage\\Pilot\\3.histones\\aap.pdf");
		
		
		/*
		WriteWigToText wwtt = new WriteWigToText();
		wwtt.writeWigToText("C:\\Users\\Matthieu\\Documents\\Afstudeerstage\\Pilot\\3.histones\\test\\wgEncodeUwHistoneGm06990H3k27me3StdRawRep2.wig",
				"C:\\Users\\Matthieu\\Documents\\Afstudeerstage\\Pilot\\3.histones\\test\\wgEncodeUwHistoneGm06990H3k27me3StdRawRep2.txt");
		
		
		SplitWigTextDataOnChromosome swtdoc = new SplitWigTextDataOnChromosome("C:\\Users\\Matthieu\\Documents\\Afstudeerstage\\Pilot\\3.histones\\test\\wgEncodeUwHistoneGm06990H3k27me3StdRawRep2.txt",
				"C:\\Users\\Matthieu\\Documents\\Afstudeerstage\\Pilot\\3.histones\\test\\Gm06990");
		
		
		Testing test = new Testing();
		test.compareWigSets("C:\\Users\\Matthieu\\Documents\\Afstudeerstage\\Pilot\\3.histones\\test\\Gm06990\\Gm06990wgEncodeUwHistoneGm06990H3k27me3StdRawRep1_chr5.txt",
				"C:\\Users\\Matthieu\\Documents\\Afstudeerstage\\Pilot\\3.histones\\test\\Gm06990\\Gm06990wgEncodeUwHistoneGm06990H3k27me3StdRawRep2_chr5.txt");
		
		
		
		ArrayList<Double> aapje1 = new ArrayList<Double>();
		aapje1.add(1.0);
		aapje1.add(2.0);
		aapje1.add(3.0);
		aapje1.add(4.0);
		aapje1.add(5.0);
		aapje1.add(6.0);
		ArrayList<Double> aapje2 = new ArrayList<Double>();
		aapje2.add(3.0);
		aapje2.add(4.0);
		aapje2.add(5.0);
		aapje2.add(6.0);
		aapje2.add(7.0);
		aapje2.add(8.0);
		
		double[] arr_aapje1 = Doubles.toArray(aapje1);
		double[] arr_aapje2 = Doubles.toArray(aapje2);
		
		for(int i=0;i<=5;i++){
			System.out.println(arr_aapje1[i] + "\t" + arr_aapje2[i]);
		}
		*/
		
	}
	
	
	
	public GenomicBoundaries readWigData(String pruikBestand)throws IOException{
		GenomicBoundaries<Object> histoneRegions = new GenomicBoundaries();
		
		WigFile wf = new WigFile(pruikBestand, false);
		long totalN = wf.size();
		long n = 0;
		while(n < totalN){
			UCSCDataObject ucscdo = wf.parseLn();
			String chr = Byte.toString(ucscdo.getChr());
			System.out.println( ucscdo.getChr() + "\t" + ucscdo.getPositionStart() + "\t" + ucscdo.getPositionEnd() + "\t" + ucscdo.getValue() );
			//histoneRegions.addBoundary(chr, ucscdo.getPositionStart(), ucscdo.getPositionEnd(), ucscdo.getValue());
		}
		
		return histoneRegions;
	}
	
	
	public void compareWigSets(String wigFilePath1, String wigFilePath2)throws IOException{
		HashSet<Integer> wigFile1 = readData(wigFilePath1);
		System.out.println("Read first wigfile.");
		
		HashSet<Integer> wigFile2 = readData(wigFilePath2);
		System.out.println("Read second wigfile.");
		
		System.out.println(wigFile1.size() + " :-: " + wigFile2.size());
		//wigFile1.addAll(wigFile2);
		wigFile1.retainAll(wigFile2);
		System.out.println(wigFile1.size());
	}
	
	
	
	
	public HashSet<Integer> readHistoneSiteDataFromWig(String wigFileLocation) throws IOException{
		HashSet<Integer> chrAndStart = new HashSet<Integer>();
		WigFile wf = new WigFile(wigFileLocation, false);
		long totalN = wf.countLines();
		long n = 0;
		while(n < totalN){
			UCSCDataObject ucscdo = wf.parseLn();
			
			if(ucscdo != null){
				String chr = Byte.toString(ucscdo.getChr());
				chrAndStart.add(ucscdo.getPositionStart());

				if(n%10000==0){
					System.out.println("Read " + n + " elements.");
				}
			}
			n++;
		}
		wf.close();
		return chrAndStart;
	}
	
	
	public HashSet<Integer> readData(String inFilePath)throws IOException{
		HashSet<Integer>  starts = new HashSet<Integer>();
		
		String fileLine;
		String[] fileLineData;
		TextFile tf = new TextFile(inFilePath, false);
		while( (fileLine=tf.readLine())!=null ){
			fileLineData = TAB_PATTERN.split(fileLine);
			starts.add(Integer.parseInt(fileLineData[1]));
		}
		tf.close();
		return starts;
	}
}
