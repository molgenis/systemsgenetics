/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package nl.systemsgenetics.eqtlpermutationtranscriptionfactoranalysis;
import java.io.IOException;
import org.bouncycastle.util.Arrays;
import umcg.genetica.genomicboundaries.GenomicBoundaries;
import umcg.genetica.graphics.ViolinBoxPlot;
import umcg.genetica.io.ucsc.UCSCDataObject;
import umcg.genetica.io.ucsc.WigFile;
import umcg.genetica.math.stats.WilcoxonMannWhitney;

/**
 *
 * @author Matthieu
 */
public class Testing {
	
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
		
		/*
		System.out.println("Testing the Genetica Library Wilcox test.");
		double[] A = new double[]{35,23,58,74,70,67,34,76,34,73,40,76,34,80,7};
		double[] B = new double[]{23,48,48,90,68,34,90,86,24,86,92,63,78,94,65,89,34,67,89,34,6};
		System.out.println("Made two arrays with values: {A: " + A.length + "} :-: {B:" + B.length + "}");
		WilcoxonMannWhitney wcmw = new WilcoxonMannWhitney();
		double testPvalue = wcmw.returnWilcoxonMannWhitneyPValue(A, B);
		double wilcoxAuc = wcmw.getAUC();
		System.out.println("PValue is: " + testPvalue + " and the UAC value is: " + wilcoxAuc);
		
		
		System.out.println("Iets met boxplots...");
		ViolinBoxPlot vbp = new ViolinBoxPlot();
		double[][][] setA = new double[1][1][15];
		setA[0][0] = A;
		System.out.println(setA[0][0]);
		
		
		
		double[] A = new double[]{35,23,58,74,70,67,34,76,34,73,40,76,34,80,7};
		double[] B = new double[]{23,48,48,90,68,34,90,86,24,86,92,63,78,94,65,89,34,67,89,34,6};
		double[][][] setA = new double[1][2][0];
		setA[0][0] = A;
		setA[0][1] = B;
		
		String[][] plotLabels = new String[][]{{"Dingen", "Dingen2"}};
		String[] dataset = new String[]{"Aap"}; 
		
		/*
		//Test to see if the values are added.
		for(int n=0;n<15;n++){
			System.out.println(setA[0][0][n]);
		}
		
		
		ViolinBoxPlot vbp = new ViolinBoxPlot();
		System.out.println("About to draw something;");
		vbp.draw(setA, dataset, plotLabels, "Iets", ViolinBoxPlot.Output.PDF, "C:\\Users\\Matthieu\\Documents\\Afstudeerstage\\Pilot\\3.histones\\test.pdf");
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
	
}
