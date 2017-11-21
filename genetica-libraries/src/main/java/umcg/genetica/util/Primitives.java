/*
 * as taken from http://pastebin.com/f57ae8d47
 */
package umcg.genetica.util;

import java.util.ArrayList;
import java.util.List;

/**
 *
 *
 */
public class Primitives {
	
	public static double max(double[] values) {
		double max = -Double.MAX_VALUE;
		for (double value : values) {
			if (!Double.isNaN(value) && value > max) {
				max = value;
			}
		}
		return max;
	}
	
	public static double max(double[][] values) {
		double max = -Double.MAX_VALUE;
		for (double[] x : values) {
			for (double value : x) {
				if (!Double.isNaN(value) && value > max) {
					max = value;
				}
			}
		}
		return max;
	}
	
	public static char max(char[] values) {
		char max = Character.MIN_VALUE;
		for (char value : values) {
			if (value > max) {
				max = value;
			}
		}
		return max;
	}
	
	public static byte max(byte[] values) {
		byte max = Byte.MIN_VALUE;
		for (byte value : values) {
			if (value > max) {
				max = value;
			}
		}
		return max;
	}
	
	public static short max(short[] values) {
		short max = Short.MIN_VALUE;
		for (short value : values) {
			if (value > max) {
				max = value;
			}
		}
		return max;
	}
	
	public static int max(int[] values) {
		int max = Integer.MIN_VALUE;
		for (int value : values) {
			if (value > max) {
				max = value;
			}
		}
		return max;
	}
	
	public static long max(long[] values) {
		long max = Long.MIN_VALUE;
		for (long value : values) {
			if (value > max) {
				max = value;
			}
		}
		return max;
	}
	
	public static char maxChar(char... values) {
		return max(values);
	}
	
	public static byte maxByte(byte... values) {
		return max(values);
	}
	
	public static short maxShort(short... values) {
		return max(values);
	}
	
	public static int maxInt(int... values) {
		return max(values);
	}
	
	public static long maxLong(long... values) {
		return max(values);
	}
	
	public static double maxDouble(double... values) {
		return max(values);
	}
	
	public static double min(double[] values) {
		double min = Double.MAX_VALUE;
		for (double value : values) {
			if (!Double.isNaN(value) && value < min) {
				min = value;
			}
		}
		return min;
	}
	
	public static double min(double[][] values) {
		double min = Double.MAX_VALUE;
		for (double[] x : values) {
			for (double value : x) {
				if (!Double.isNaN(value) && value < min) {
					min = value;
				}
			}
		}
		return min;
	}
	
	public static char min(char[] values) {
		char min = Character.MAX_VALUE;
		for (char value : values) {
			if (value < min) {
				min = value;
			}
		}
		return min;
	}
	
	public static byte min(byte[] values) {
		byte min = Byte.MAX_VALUE;
		for (byte value : values) {
			if (value < min) {
				min = value;
			}
		}
		return min;
	}
	
	public static short min(short[] values) {
		short min = Short.MAX_VALUE;
		for (short value : values) {
			if (value < min) {
				min = value;
			}
		}
		return min;
	}
	
	public static int min(int[] values) {
		int min = Integer.MAX_VALUE;
		for (int value : values) {
			if (value < min) {
				min = value;
			}
		}
		return min;
	}
	
	public static long min(long[] values) {
		long min = Long.MAX_VALUE;
		for (long value : values) {
			if (value < min) {
				min = value;
			}
		}
		return min;
	}
	
	public static double minDouble(double... values) {
		return min(values);
	}
	
	public static char minChar(char... values) {
		return min(values);
	}
	
	public static byte minByte(byte... values) {
		return min(values);
	}
	
	public static short minShort(short... values) {
		return min(values);
	}
	
	public static int minInt(int... values) {
		return min(values);
	}
	
	public static long minLong(long... values) {
		return min(values);
	}
	
	public static int getIndexOfMin(double[] mse) {
		
		double min = mse[0];
		int output = 0;
		for (int i = 0; i < mse.length; i++) {
			if (mse[i] < min) {
				output = i;
			}
		}
		return output;
	}
	
	public static double[] toPrimitiveArr(Double[] toArray) {
		double[] arr = new double[toArray.length];
		for (int i = 0; i < toArray.length; i++) {
			arr[i] = toArray[i];
		}
		return arr;
	}
	
	public static int[] toPrimitiveArr(Integer[] toArray) {
		int[] arr = new int[toArray.length];
		for (int i = 0; i < toArray.length; i++) {
			arr[i] = toArray[i];
		}
		return arr;
	}
	
	public static int[] toPrimitiveArr(List<Integer> selectN) {
		return toPrimitiveArr(selectN.toArray(new Integer[0]));
	}
	
	public static double[] toPrimitiveArr(ArrayList<Double> selectZ) {
		return toPrimitiveArr(selectZ.toArray(new Double[0]));
	}
}
