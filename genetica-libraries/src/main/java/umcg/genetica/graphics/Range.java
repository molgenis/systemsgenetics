package umcg.genetica.graphics;

import JSci.maths.ArrayMath;

/**
 * Created by hwestra on 7/16/15.
 */
public class Range {
	private double maxX;
	private double maxY;
	private double minX;
	private double minY;
	private Double unitY;
	
	@Override
	public String toString() {
		return "Range{" +
				"maxX=" + maxX +
				", maxY=" + maxY +
				", minX=" + minX +
				", minY=" + minY +
				", unitY=" + unitY +
				", unitX=" + unitX +
				'}';
	}
	
	private Double unitX;
	
	public Range(double minX, double minY, double maxX, double maxY) {
		this.minX = minX;
		this.minY = minY;
		this.maxX = maxX;
		this.maxY = maxY;
	}
	
	public Range(double[] x, double[] y) {
		
		minX = ArrayMath.min(x);
		maxX = ArrayMath.max(x);
		minY = ArrayMath.min(y);
		maxY = ArrayMath.max(y);
		
		
	}
	
	public Range(double[][] dataY) {
		
		minX = 0;
		maxX = dataY[0].length;
		
		for (int i = 0; i < dataY.length; i++) {
			double min = ArrayMath.min(dataY[i]);
			if (min < minY) {
				minY = min;
			}
			double max = ArrayMath.max(dataY[i]);
			if (max > maxY) {
				maxY = max;
			}
		}
	}
	
	public Range(double[][][] dataY) {
		maxY = -Double.MAX_VALUE;
		minY = Double.MAX_VALUE;
		
		for (int r = 0; r < dataY.length; r++) {
			for (int c = 0; c < dataY[r].length; c++) {
				for (int q = 0; q < dataY[r][c].length; q++) {
					double v = dataY[r][c][q];
					if (v > maxY) {
						maxY = v;
					}
					if (v < minY) {
						minY = v;
					}
				}
			}
		}
	}
	
	public double getMaxX() {
		return maxX;
	}
	
	public double getMaxY() {
		return maxY;
	}
	
	public double getMinX() {
		return minX;
	}
	
	public double getMinY() {
		return minY;
	}
	
	public void roundX() {
		// round up max X
		double rangeX = Math.abs(maxX - minX);
		unitX = determineUnit(rangeX);
		double remainder = Math.abs(maxX) % unitX;
		if (remainder != 0) {
			maxX += (unitX - remainder);
		}
		
		// round down min X
		remainder = Math.abs(minX) % unitX;
		System.out.println(remainder+" remainder minx");
		if (remainder != 0) {
			minX -= (unitX - remainder);
		}
	}
	
	public void roundY() {
		// round up max Y
		double rangeY = Math.abs(maxY - minY);
		unitY = determineUnit(rangeY);
		double remainder = Math.abs(maxY) % unitY;
		if (remainder != 0) {
			maxY += (unitY - remainder);
		}
		
		// round down min Y
		remainder = Math.abs(minY) % unitY;
		if (remainder != 0) {
			minY -= (unitY - remainder);
		}
	}
	
	public void round() {
		roundX();
		roundY();
	}
	
	public double getRangeX() {
		double rangeX = Math.abs(maxX - minX);
		return rangeX;
	}
	
	public double getRangeY() {
		double rangeY = Math.abs(maxY - minY);
		return rangeY;
	}
	
	public double determineUnit(double range) {
		
		double divisor = Math.log10(range);
		divisor = Math.floor(divisor);
		divisor = Math.pow(10, divisor);
		return divisor;
	}
	
	public double getRelativePositionX(double x) {
		double tmpmaxX = maxX;
		
		// move to 0..maxX
		double tmpX = x;
		if (minX < 0) {
			tmpmaxX += (-1 * minX);
			tmpX += (-1 * minX);
		} else {
			tmpmaxX -= minX;
			tmpX -= minX;
		}
		return tmpX / tmpmaxX;
	}
	
	public double getRelativePositionY(double y) {
		double tmpmaxY = maxY;
		
		// move to 0..maxX
		double tmpY = y;
		if (minY < 0) {
			tmpmaxY += (-1 * minY);
			tmpY += (-1 * minY);
		} else {
			tmpmaxY -= minY;
			tmpY -= minY;
		}
		return tmpY / tmpmaxY;
	}
	
	public double getUnitX() {
		if (unitX == null) {
			unitX = determineUnit(getRangeX());
		}
		return unitX;
	}
	
	public double getUnitY() {
		if (unitY == null) {
			unitY = determineUnit(getRangeY());
		}
		return unitY;
	}
}
