package nl.umcg.westrah.binarymetaanalyzer;

import org.apache.commons.io.FilenameUtils;
import umcg.genetica.io.bin.BinaryFile;
import umcg.genetica.io.text.TextFile;

import java.io.*;
import java.util.Arrays;
import java.util.zip.GZIPInputStream;

/**
 * Created by harmbrugge on 30-06-17.
 */
public class MatrixConverter {
	
	
	private static final String NR_CALLED = "5075";
	private static final String MAF = "0.5";
	private static final String HWE = "1.0";
	private static final String CALL_RATE = "1.0";
	
	private File inputFile;
	private String directory;
	private String baseName;
	private int geneCount;
	
	private BinaryFile binaryOutputFile;
	private TextFile snpFile;
	private TextFile columnFile;
	
	public MatrixConverter(File inputMatrix) throws IOException {
		this.inputFile = inputMatrix;
		this.directory = inputMatrix.getParent();
		this.baseName = FilenameUtils.getBaseName(inputMatrix.getName());
		this.binaryOutputFile = new BinaryFile(directory + "/" + baseName + ".dat", BinaryFile.W);
		this.snpFile = new TextFile(directory + "/" + baseName + "-RowNames.txt.gz", TextFile.W);
		this.columnFile = new TextFile(directory + "/" + baseName + "-ColNames.txt.gz", TextFile.W);
		
		snpFile.writeln("SNP\tAlleles\tMinorAllele\tAlleleAssessed\tNrCalled\tMaf\tHWE\tCallRate");
		
	}
	
	public void convert() throws IOException {
		// TODO: Called 'Magic number' in comments, what is this number? if (m_cisOnly) 1 else 0
		binaryOutputFile.writeInt(0);
		
		InputStream fileStream = new FileInputStream(inputFile);
		InputStream gzipStream = new GZIPInputStream(fileStream);
		Reader decoder = new InputStreamReader(gzipStream);
		
		try (BufferedReader br = new BufferedReader(new BufferedReader(decoder))) {
			
			String line = br.readLine();
			String[] header = line.split("\\s");
			String[] genes = Arrays.copyOfRange(header, 3, header.length);
			writeColumnFile(genes);
			
			line = br.readLine();
			
			while (line != null) {
				String[] splitLine = line.split("\\s");
				writeSnpLine(splitLine);
				writeDatFile(splitLine);
				line = br.readLine();
			}
		} catch (IOException e) {
			e.printStackTrace();
		} finally {
			snpFile.close();
			columnFile.close();
			binaryOutputFile.close();
		}
	}
	
	private void writeDatFile(String[] splitLine) throws IOException {
		for (int i = 3; i < splitLine.length; i++) {
			float zScore = Float.parseFloat(splitLine[i]);
			binaryOutputFile.writeFloat(zScore);
		}
	}
	
	private void writeColumnFile(String[] genes) throws IOException {
		geneCount = genes.length;
		for (String gene :genes) columnFile.writeln(gene);
		columnFile.close();
	}
	
	private void writeSnpLine(String[] splitLine) throws IOException {
		
		String snpName = splitLine[0];
		String alleles = splitLine[1] + "/" + splitLine [2];
		String minorAllele = splitLine[2];
		String geneCount = Integer.toString(this.geneCount);
		
		snpFile.writelnTabDelimited(new String[]{snpName, alleles, minorAllele, minorAllele, NR_CALLED, MAF, HWE, CALL_RATE, geneCount, "-"});
		
	}
	
}