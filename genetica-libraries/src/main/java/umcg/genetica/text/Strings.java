/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package umcg.genetica.text;

import gnu.trove.map.hash.THashMap;

import java.io.UnsupportedEncodingException;
import java.text.DecimalFormat;
import java.util.ArrayList;
import java.util.List;
import java.util.Set;
import java.util.function.BiConsumer;
import java.util.function.BinaryOperator;
import java.util.function.Function;
import java.util.function.Supplier;
import java.util.regex.Matcher;
import java.util.regex.Pattern;
import java.util.stream.Collector;

/**
 * @author harmjan
 */
public class Strings {
	
	public static final Pattern tab = Pattern.compile("\t");
	public static final Pattern whitespace = Pattern.compile("\\s");
	public static final Pattern comma = Pattern.compile(",");
	public static final Pattern semicolon = Pattern.compile(";");
	public static final Pattern colon = Pattern.compile(":");
	public static final Pattern pipe = Pattern.compile("\\|");
	public static final Pattern forwardslash = Pattern.compile("/");
	public static final Pattern backwardslash = Pattern.compile("\\\\");
	public static final Pattern dot = Pattern.compile("\\.");
	public static final Pattern space = Pattern.compile(" ");
	public static final Pattern dash = Pattern.compile("-");
	
	public static String concat(String[] s, Pattern t) {
		String sepstr = t.toString();
		int len = 0;
		for (int i = 0; i < s.length; i++) {
			if (s[i] != null) {
				len += s[i].length();
			} else {
				len += 4;
			}
		}
		len += ((s.length - 1) * sepstr.length());
		
		StringBuilder output = new StringBuilder(len);
		for (int i = 0; i < s.length; i++) {
			if (i == 0) {
				output.append(s[i]);
			} else {
				output.append(sepstr).append(s[i]);
			}
		}
		return output.toString();
	}
	
	public static String concat(Object[] s, Pattern t) {
		String sepstr = t.toString();
		int len = 0;
		for (int i = 0; i < s.length; i++) {
			if (s[i] != null) {
				len += s[i].toString().length();
			} else {
				len += 4;
			}
		}
		len += ((s.length - 1) * t.toString().length());
		
		
		StringBuilder output = new StringBuilder(len);
		for (int i = 0; i < s.length; i++) {
			if (i == 0) {
				output.append(s[i].toString());
			} else {
				output.append(sepstr).append(s[i].toString());
			}
		}
		return output.toString();
	}
	
	public static String concat(double[] s, Pattern t) {
		String[] str = new String[s.length];
		for (int i = 0; i < s.length; i++) {
			str[i] = "" + s[i];
		}
		return concat(str, t);
	}
	
	public static String concat(double[] s, DecimalFormat f, Pattern t) {
		String[] str = new String[s.length];
		for (int i = 0; i < s.length; i++) {
			if (Double.isNaN(s[i])) {
				str[i] = "" + Double.NaN;
			} else {
				str[i] = "" + f.format(s[i]);
			}
			
		}
		return concat(str, t);
	}
	
	public static String concat(float[] s, DecimalFormat f, Pattern t) {
		String[] str = new String[s.length];
		for (int i = 0; i < s.length; i++) {
			str[i] = "" + f.format(s[i]);
		}
		return concat(str, t);
	}
	
	public static String concat(int[] s, Pattern t) {
		String[] str = new String[s.length];
		for (int i = 0; i < s.length; i++) {
			str[i] = "" + s[i];
		}
		return concat(str, t);
	}
	
	public static String concat(List<String> s, Pattern t) {
		String[] data = s.toArray(new String[0]);
		return concat(data, t);
	}
	
	public static String concat(String[] s, Pattern t, int start, int end) {
		String[] data = new String[end - start];
		for (int i = start; i < end; i++) {
			data[i - start] = s[i];
		}
		return concat(data, t);
	}
	
	public static String[] split(String in) {
		List<String> list = new ArrayList<String>();
		StringBuilder sb = new StringBuilder();
		int len = in.length();
		int i = 0;
		char c;
		while (i < len) {
			c = in.charAt(i);
			if (c == '\t' || i == len) {
				list.add(sb.toString());
				sb.delete(0, len - 1);
			} else {
				sb.append(c);
			}
			i++;
		}
		
		return list.toArray(new String[0]);
	}
	
	public static String reverse(String str) {
		return new StringBuffer(str).reverse().toString();
	}
	
	// split a string, and match against a certain other string
	public static int[] countOccurrences(String ln, Pattern tab, String match) {
		int totalElems = 0;
		int matchingElems = 0;
		Matcher m = tab.matcher(ln);
		int prevsta = 0;
		while (m.find()) {
			String substr = ln.substring(prevsta, m.start());
			if (match.equals(substr)) {
				matchingElems++;
			}
			totalElems++;
			prevsta = m.end();
		}
		String substr = ln.substring(prevsta, ln.length());
		if (match.equals(substr)) {
			matchingElems++;
		}
		totalElems++;
		
		return new int[]{totalElems, matchingElems};
	}
	
	// split a string, but only return a part of the output
	public static String[] subsplit(String ln, Pattern tab, int lowerBound, int upperBound) {
		
		String[] output = new String[upperBound - lowerBound];
		Matcher m = tab.matcher(ln);
		
		int ctr1 = 0;
		int ctr2 = 0;
		
		int prevsta = 0;
		
		while (m.find()) {
//			System.out.println(m.group() + "\t" + m.start() + "\t" + m.end() + "\t" + ln.substring(prevsta, m.start()));
			
			if (ctr1 >= lowerBound && ctr1 != upperBound) {
				String substr = ln.substring(prevsta, m.start());
//				System.out.println("inner loop " + ctr2 + "\t" + substr);
				output[ctr2] = substr;
				ctr2++;
			}
			if (ctr1 >= upperBound) {
//				System.out.println("breaking");
				break;
			}
//			System.out.println(ctr1);
			prevsta = m.end();
			ctr1++;
		}
		
		if (ctr1 < upperBound && prevsta < ln.length()) {
			// last element
			String substr = ln.substring(prevsta, ln.length());
//			System.out.println("outer loop " + ctr2 + "\t" + substr);
			output[ctr2] = substr;
			ctr2++;
		}
//		System.out.println("end ctr2: " + ctr2);
//		System.out.println(upperBound);
//		System.out.println("enc ctr1: " + ctr1);
		if (ctr2 == 0) {
			return null;
		} else {
			return output;
		}
		
	}
	
	public static int countSeparatorOccurrences(String permln, Pattern tab) {
		Matcher m = tab.matcher(permln);
		int ctr = 0;
		while (m.find()) {
			ctr++;
		}
		return ctr;
	}
	
	private static THashMap<String, String> cache;
	
	public synchronized static String cache(String s) {
		if (cache == null) {
			cache = new THashMap<>(1000000);
		}
		if (!cache.containsKey(s)) {
			try {
				String str = new String(s.getBytes("UTF-8"));
				cache.put(str, str);
				return str;
			} catch (UnsupportedEncodingException e) {
				return s;
			}
		} else {
			return cache.get(s);
		}
	}
	
	
}
