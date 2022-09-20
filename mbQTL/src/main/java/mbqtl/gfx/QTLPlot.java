package mbqtl.gfx;

import com.itextpdf.text.DocumentException;
import umcg.genetica.graphics.DefaultGraphics;
import umcg.genetica.graphics.Grid;

import java.io.IOException;

public class QTLPlot {

	public static void main(String[] args) {

		// generate some test data
		double[] x = new double[1000];
		double[] y = new double[1000];
		for (int i = 0; i < 1000; i++) {
			if (i < 300) {
				x[i] = 0 + Math.random() * 0.4;
				x[i] = 0;
				y[i] = Math.random() * 500;
			} else if (i < 800) {
				if (Math.random() > 0.5) {
					x[i] = 1 - (Math.random() * 0.9);
				} else {
					x[i] = 1 + (Math.random() * 0.9);
				}
				x[i] = 1;
				y[i] = Math.random() * 500 + 400;
			} else {
				x[i] = 2 - Math.random() * 0.4;
				x[i] = 2;
				y[i] = Math.random() * 500 + 800;
			}
		}

		Grid g = new Grid(200, 200, 1, 3, 50, 50);
		QTLPanel panel = new QTLPanel(1, 1);
		panel.setData(x, y);
		double z = 5.57855;
		double p = 1e-10;
		double r = 0.35797;
		panel.setAlleles(new String[]{"C", "T"});
		panel.setDatasetDetails("Dataset1", "ENSG001", "rs0123", z, p, r);
		g.addPanel(panel);

		panel = new QTLPanel(1, 1);
		panel.setData(x, y);
		panel.setAlleles(new String[]{"C", "T"});
		z = 10.0587;
		p = 1e-232;
		r = 0.546874;
		panel.setDatasetDetails("Dataset2", "ENSG001", "rs0123", z, p, r);
		g.addPanel(panel);

		panel = new QTLPanel(1, 1);
		z = 10.0587;
		p = 1e-232;
		r = 0.546874;
		panel.setDatasetDetails("Dataset3", "ENSG001", "rs0123", z, p, r);
		panel.setNotTested();
		g.addPanel(panel);

		try {
			g.draw("d:\\test.pdf", DefaultGraphics.Output.PDF);
		} catch (IOException e) {
			throw new RuntimeException(e);
		} catch (DocumentException e) {
			throw new RuntimeException(e);
		}

	}


}
