/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package umcg.genetica.io.trityper.probeannotation;

import java.util.ArrayList;
import java.util.HashMap;

/**
 *
 * @author harm-jan
 */
public class ProbeAnnotation {
    HashMap<Integer, Integer> probeTranslation = new HashMap<Integer, Integer>();
    String name;
    HashMap<String, Integer> seqToProbe = new HashMap<String, Integer>();
    ArrayList<String> ilmn = new ArrayList<String>();
    ArrayList<Integer> arrayId = new ArrayList<Integer>();
    ArrayList<String> seq = new ArrayList<String>();
    ArrayList<String> chr = new ArrayList<String>();
    ArrayList<String> chrpos = new ArrayList<String>();
    ArrayList<String> symbol = new ArrayList<String>();
    String[] annotation;
}
