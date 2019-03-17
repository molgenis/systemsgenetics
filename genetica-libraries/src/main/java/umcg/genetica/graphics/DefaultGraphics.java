package umcg.genetica.graphics;

import com.itextpdf.text.Document;
import com.itextpdf.text.DocumentException;
import com.itextpdf.text.pdf.PdfContentByte;
import com.itextpdf.text.pdf.PdfWriter;
import umcg.genetica.graphics.themes.DefaultTheme;

import javax.imageio.ImageIO;
import java.awt.*;
import java.awt.image.BufferedImage;
import java.io.File;
import java.io.FileNotFoundException;
import java.io.IOException;
import java.util.Locale;


/**
 * @author hwestra
 */
public class DefaultGraphics {
	
	protected Output output = Output.PDF;
	protected BufferedImage bi;
	protected Graphics2D g2d;
	protected Document document;
	protected PdfWriter writer;
	protected int figureWidth;
	protected int figureHeight;
	protected String outputFileName;
	private Locale defaultLocale;
	private PdfContentByte cb;
	
	
	public enum Output {
		
		PDF, PNG
	}
	
	protected DefaultGraphics() {
	}
	
	public DefaultGraphics(String outputFileName, int width, int height) throws FileNotFoundException, DocumentException {
		this.initializePlot(outputFileName, width, height);
	}
	
	public Output getOutput() {
		return output;
	}
	
	public void setOutput(Output output) {
		this.output = output;
	}
	
	public BufferedImage getBi() {
		return bi;
	}
	
	public void setBi(BufferedImage bi) {
		this.bi = bi;
	}
	
	public void setG2d(Graphics2D g2d) {
		this.g2d = g2d;
	}
	
	public Document getDocument() {
		return document;
	}
	
	public void setDocument(Document document) {
		this.document = document;
	}
	
	public PdfWriter getWriter() {
		return writer;
	}
	
	public void setWriter(PdfWriter writer) {
		this.writer = writer;
	}
	
	public int getFigureWidth() {
		return figureWidth;
	}
	
	public void setFigureWidth(int figureWidth) {
		this.figureWidth = figureWidth;
	}
	
	public int getFigureHeight() {
		return figureHeight;
	}
	
	public void setFigureHeight(int figureHeight) {
		this.figureHeight = figureHeight;
	}
	
	public String getOutputFileName() {
		return outputFileName;
	}
	
	public void setOutputFileName(String outputFileName) {
		this.outputFileName = outputFileName;
	}
	
	public Locale getDefaultLocale() {
		return defaultLocale;
	}
	
	public void setDefaultLocale(Locale defaultLocale) {
		this.defaultLocale = defaultLocale;
	}
	
	public PdfContentByte getCb() {
		return cb;
	}
	
	public void setCb(PdfContentByte cb) {
		this.cb = cb;
	}
	
	protected void initializePlot(String outputFileName, int width, int height) throws DocumentException, FileNotFoundException {
		if (outputFileName.toLowerCase().endsWith("png")) {
			output = Output.PNG;
		}
		defaultLocale = Locale.getDefault();
		Locale.setDefault(Locale.US);
		// set up Graphics2D depending on required format using iText in case PDF
		g2d = null;
		document = null;
		writer = null;
		bi = null;
		this.outputFileName = outputFileName;
		bi = new BufferedImage(1, 1, BufferedImage.TYPE_INT_RGB);
		g2d = bi.createGraphics();
		
		figureWidth = width;
		figureHeight = height;
		
		// initialize plot
		if (output == Output.PDF) {
			com.itextpdf.text.Rectangle rectangle = new com.itextpdf.text.Rectangle(width, height);
			document = new Document(rectangle);
			writer = com.itextpdf.text.pdf.PdfWriter.getInstance(document, new java.io.FileOutputStream(outputFileName));
			
			document.open();
			cb = writer.getDirectContent();
			cb.saveState();
			
			g2d = cb.createGraphics(width, height);
		} else {
			bi = new BufferedImage(width, height, BufferedImage.TYPE_INT_RGB);
			g2d = bi.createGraphics();
		}
		
		g2d.setRenderingHint(RenderingHints.KEY_ANTIALIASING, RenderingHints.VALUE_ANTIALIAS_ON);
		g2d.setColor(Color.white);
		g2d.fillRect(0, 0, width, height);
		
		g2d.setStroke(new DefaultTheme().getStroke());
	}
	
	public void close() throws IOException {
		// dispose
		g2d.dispose();
		if (output == Output.PDF) {
			
			g2d.dispose();
			cb.restoreState();
			document.close();
			writer.close();
		} else {
			bi.flush();
			ImageIO.write(bi, "PNG", new File(outputFileName));
		}
		
		Locale.setDefault(defaultLocale);
	}
	
	public Graphics2D getG2d() {
		return g2d;
	}
	
	public Color interpolateColor(Color color1, Color color2, double blending) {
		double inverse_blending = 1 - blending;
		int red = (int) Math.floor(color1.getRed() * blending + color2.getRed() * inverse_blending);
		int green = (int) Math.floor(color1.getGreen() * blending + color2.getGreen() * inverse_blending);
		int blue = (int) Math.floor(color1.getBlue() * blending + color2.getBlue() * inverse_blending);
		int alpha = color1.getAlpha();

//note that if i pass float values they have to be in the range of 0.0-1.0
//and not in 0-255 like the ones i get returned by the getters.
		Color blended = new Color(red, green, blue, alpha);
		return blended;
	}
	
	public double determineUnit(double range) {
		
		double divisor = Math.log10(range);
		divisor = Math.floor(divisor);
		divisor = Math.pow(10, divisor);
		return divisor;
	}
	
	
}
