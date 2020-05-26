package umcg.genetica.graphics.panels;


import umcg.genetica.graphics.DefaultGraphics;
import umcg.genetica.graphics.themes.DefaultTheme;
import umcg.genetica.graphics.themes.Theme;

import java.awt.*;

/**
 * Created by hwestra on 7/15/15.
 */
public abstract class Panel {
	protected int width = 0;
	protected int height = 0;
	protected int marginX = 0;
	protected int marginY = 0;
	protected int x0 = 0;
	protected int y0 = 0;
	protected int nrRows = 0;
	protected int nrCols = 0;
	protected Theme theme = new DefaultTheme();
	protected String title;

	public abstract void draw(DefaultGraphics g);

	public Panel(int nrRows, int nrCols) {
		this.nrRows = nrRows;
		this.nrCols = nrCols;
	}

	public void setDimensions(int panelWidth, int panelHeight) {
		this.width = panelWidth;
		this.height = panelHeight;
	}

	public void drawRotate(Graphics2D g2d, double x, double y, double angle, String text) {
		g2d.translate((float) x, (float) y);
		g2d.rotate(Math.toRadians(angle));
		g2d.drawString(text, 0, 0);
		g2d.rotate(-Math.toRadians(angle));
		g2d.translate(-(float) x, -(float) y);
	}

	public int getStringWidth(Graphics2D g2d, String str) {
		FontMetrics metrics = g2d.getFontMetrics(g2d.getFont());
		return metrics.stringWidth(str);

	}

	public void setTitle(String title) {
		this.title = title;
	}

	public void setX0(int x0) {
		this.x0 = x0;
	}

	public void setY0(int y0) {
		this.y0 = y0;
	}

	public int getNrRows() {
		return nrRows;
	}

	public int getNrCols() {
		return nrCols;
	}

	public int getMarginX() {
		return marginX;
	}

	public void setMarginX(int marginX) {
		this.marginX = marginX;
	}

	public int getMarginY() {
		return marginY;
	}

	public void setMarginY(int marginY) {
		this.marginY = marginY;
	}

	public void setTheme(Theme theme) {
		this.theme = theme;
	}


}
