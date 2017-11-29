package eqtlmappingpipeline.util;

import umcg.genetica.io.text.TextFile;
import umcg.genetica.io.trityper.ProbeAnnotation;

import java.io.IOException;
import java.util.HashMap;
import java.util.HashSet;

public class ExpressionSampleFilter {
	
	public void filter(String efileIn, String efileOut, String sampleList) throws IOException {
		
		
		TextFile tf = new TextFile(sampleList, TextFile.R);
		String[] elems = tf.readLineElems(TextFile.tab);
		HashMap<String, String> sampleToCohort = new HashMap<>();
		HashSet<String> cohorts = new HashSet<String>();
		while (elems != null) {
			
			if (elems.length > 1) {
				// sample\tcohort
				sampleToCohort.put(elems[0], elems[1]);
				cohorts.add(elems[1]);
			} else {
				sampleToCohort.put(elems[0], null);
			}
			
			elems = tf.readLineElems(TextFile.tab);
		}
		tf.close();
		
		System.out.println(sampleToCohort.size() + " sample annotations loaded");
		if (cohorts != null) {
			System.out.println(cohorts.size() + " cohorts");
		}
		int[] include = null;
		TextFile[] outarr = null;
		
		if (cohorts == null) {
			TextFile tf2 = new TextFile(efileIn, TextFile.R);
			String[] header = tf2.readLineElems(TextFile.tab);
			include = new int[header.length];
			
			String headerOut = header[0];
			
			for (int i = 1; i < header.length; i++) {
				if (sampleToCohort.containsKey(header[i])) {
					include[i] = 0;
					headerOut += "\t" + header[i];
				} else {
					include[i] = -1;
				}
			}
			tf2.close();
			outarr = new TextFile[1];
			System.out.println("Creating: " + efileOut);
			outarr[0] = new TextFile(efileOut, TextFile.W);
			outarr[0].writeln(headerOut);
		} else {
			HashMap<String, Integer> cohortToTxt = new HashMap<String, Integer>();
			int ctr = 0;
			outarr = new TextFile[cohorts.size()];
			for (String cohort : cohorts) {
				cohortToTxt.put(cohort, ctr);
				System.out.println("Creating: " + efileOut + cohort + ".txt.gz");
				outarr[ctr] = new TextFile(efileOut + cohort + ".txt.gz", TextFile.W);
				ctr++;
			}
			
			TextFile tf2 = new TextFile(efileIn, TextFile.R);
			String[] header = tf2.readLineElems(TextFile.tab);
			include = new int[header.length];
			String[] outputHeaders = new String[cohorts.size()];
			for (int i = 0; i < outputHeaders.length; i++) {
				outputHeaders[i] = "-";
			}
			int[] samplesPerCohort = new int[cohorts.size()];
			for (int i = 1; i < header.length; i++) {
				if (sampleToCohort.containsKey(header[i])) {
					String cohort = sampleToCohort.get(header[i]);
					Integer id = cohortToTxt.get(cohort);
					include[i] = id;
					samplesPerCohort[id]++;
					outputHeaders[id] += "\t" + header[i];
				} else {
					include[i] = -1;
				}
			}
			tf2.close();
			for (int i = 0; i < outputHeaders.length; i++) {
				outarr[i].writeln(outputHeaders[i]);
			}
			System.out.println("Final samples per cohort: ");
			int sum = 0;
			for (String cohort : cohorts) {
				Integer id = cohortToTxt.get(cohort);
				
				System.out.println(cohort + "\t" + samplesPerCohort[id]);
				sum += samplesPerCohort[id];
			}
			System.out.println("Total:\t" + sum);
		}
		
		// now parse the matrix
		System.out.println("Parsing: " + efileIn);
		TextFile tf2 = new TextFile(efileIn, TextFile.R);
		tf2.readLineElems(TextFile.tab); // skip header
		elems = tf2.readLineElems(TextFile.tab);
		int ctr = 0;
		while (elems != null) {
			for (int i = 0; i < outarr.length; i++) {
				outarr[i].append(elems[0]);
			}
			for (int i = 1; i < elems.length; i++) {
				int id = include[i];
				if (id > -1) {
					outarr[id].append("\t");
					outarr[id].append(elems[i]);
				}
			}
			for (int i = 0; i < outarr.length; i++) {
				outarr[i].append("\n");
			}
			ctr++;
			if (ctr % 1000 == 0) {
				System.out.println(ctr + " lines parsed.");
			}
			elems = tf2.readLineElems(TextFile.tab);
		}
		tf2.close();
		for (int i = 0; i < outarr.length; i++) {
			outarr[i].close();
		}
		System.out.println("Done.");
		System.out.println();
	}
}
