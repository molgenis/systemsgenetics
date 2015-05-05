/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package umcg.genetica.text;

import java.text.DecimalFormat;
import java.util.ArrayList;
import java.util.List;
import java.util.regex.Pattern;

/**
 *
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

        StringBuilder output = new StringBuilder();
        for (int i = 0; i < s.length; i++) {
            if (i == 0) {
                output.append(s[i]);
            } else {
                output.append(t.toString()).append(s[i]);
            }
        }
        return output.toString();
    }

    public static String concat(Object[] s, Pattern t) {

        StringBuilder output = new StringBuilder();
        for (int i = 0; i < s.length; i++) {
            if (i == 0) {
                output.append(s[i].toString());
            } else {
                output.append(t.toString()).append(s[i].toString());
            }
        }
        return output.toString();
    }

    public static String concat(double[] s, Pattern t) {

        StringBuilder output = new StringBuilder();
        for (int i = 0; i < s.length; i++) {
            if (i == 0) {
                output.append(s[i]);
            } else {
                output.append(t.toString()).append(s[i]);
            }
        }
        return output.toString();
    }

    public static String concat(double[] s, DecimalFormat f, Pattern t) {

        StringBuilder output = new StringBuilder();
        for (int i = 0; i < s.length; i++) {
            if (i == 0) {
                output.append(f.format(s[i]));
            } else {
                output.append(t.toString()).append(f.format(s[i]));
            }
        }
        return output.toString();
    }

    public static String concat(float[] s, DecimalFormat f, Pattern t) {

        StringBuilder output = new StringBuilder();
        for (int i = 0; i < s.length; i++) {
            float floatVal = s[i];
            String str = f.format(floatVal);
            if (Float.isNaN(floatVal)) {
                str = "NaN";
            }

            if (i == 0) {
                output.append(str);
            } else {
                output.append(t.toString()).append(str);
            }
        }
        return output.toString();
    }

    public static String concat(int[] s, Pattern t) {

        StringBuilder output = new StringBuilder();
        for (int i = 0; i < s.length; i++) {
            if (i == 0) {
                output.append(s[i]);
            } else {
                output.append(t.toString()).append(s[i]);
            }
        }
        return output.toString();
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
}
