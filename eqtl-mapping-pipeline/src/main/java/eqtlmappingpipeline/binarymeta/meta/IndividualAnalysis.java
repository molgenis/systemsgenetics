/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package eqtlmappingpipeline.binarymeta.meta;

import umcg.genetica.io.trityper.probeannotation.ProbeTranslation;
import eqtlmappingpipeline.binarymeta.meta.graphics.MultiVenn;
import java.io.IOException;
import java.util.HashSet;
import java.util.zip.DataFormatException;
import umcg.genetica.io.trityper.EQTL;
import umcg.genetica.io.trityper.QTLTextFile;

/**
 *
 * @author harm-jan
 */
public class IndividualAnalysis extends MetaAnalyze {

    @Override
    public void init(String args, String texttoreplace, String replacetextwith) throws IOException {
	m_settings = new MetaSettings();
	m_settings.parse(args, texttoreplace, replacetextwith);
	probeTranslation = new ProbeTranslation();
	probeTranslation.load(m_settings.getProbetranslationfile());

    }

    @Override
    public void analyze() throws IOException, DataFormatException, Exception {
	System.out.println("");
	System.out.println("Starting analysis!");

	if (!m_settings.getOutput().endsWith("/")) {
	    m_settings.setOutput(m_settings.getOutput() + "/");
	}

	String origOuputDir = m_settings.getOutput();


	System.out.println(m_settings.getDatasetnames().size() + " datasets to analyze");

	int numTotalDatasets = m_settings.getDatasetnames().size();
//	for (int d = 0; d < numTotalDatasets; d++) {
//	    String[] locations = new String[1];
//	    locations[0] = m_settings.getDatasetlocations().get(d);
//	    String[] datasets = new String[1];
//	    datasets[0] = m_settings.getDatasetnames().get(d);
//
//	    m_settings.setOutput(origOuputDir + "/" + datasets[0] + "/");
//
//	    if (!Gpio.exists(m_settings.getOutput())) {
//		Gpio.createDir(m_settings.getOutput());
//	    }
//
//	    for (int perm = 0; perm < m_settings.getNrPermutations() + 1; perm++) {
//		ds = new Dataset[1];
//		runCalculationRound(perm, locations, datasets, d);
//	    }
//
//	    if (m_settings.getNrPermutations() > 0) {
//		FDR.calculateFDR(m_settings.getOutput(), m_settings.getNrPermutations(), m_settings.getFinalEQTLBufferMaxLength(), m_settings.getFdrthreshold());
//		EQTLDotPlot edp = new EQTLDotPlot();
//		edp.draw(m_settings.getOutput() + "/eQTLsFDR" + m_settings.getFdrthreshold() + ".txt", m_settings.getOutput() + "/DotPlot-FDR" + m_settings.getFdrthreshold() + ".png"); // "/eQTLsFDR" + fdrCutOff + ".txt", outputReportsDir + "/eQTLsFDR" + fdrCutOff + "DotPlot.png"
//		edp = null;
//	    }
//	}

	String[] setnames = new String[m_settings.getDatasetnames().size() + 1];
	for (int i = 0; i < setnames.length - 1; i++) {
	    setnames[i] = m_settings.getDatasetnames().get(i);
	}

	setnames[setnames.length - 1] = "Meta-Analysis";

	EQTL[][] eqtls = new EQTL[setnames.length][0];
	double[] setsizes = new double[setnames.length];
	for (int i = 0; i < setnames.length - 1; i++) {
	    QTLTextFile etf = new QTLTextFile(origOuputDir + "/" + setnames[i] + "/eQTLsFDR" + m_settings.getFdrthreshold() + ".txt", QTLTextFile.R);
	    eqtls[i] = etf.read();
	    HashSet<String> uniqueprobes = new HashSet<String>();
	    for (EQTL u : eqtls[i]) {
		uniqueprobes.add(u.getProbe());
	    }
	    setsizes[i] = uniqueprobes.size();
	    etf.close();
	}

	QTLTextFile etf = new QTLTextFile(origOuputDir + "/eQTLsFDR" + m_settings.getFdrthreshold() + ".txt", QTLTextFile.R);
	eqtls[setnames.length - 1] = etf.read();
	HashSet<String> uniqueprobes = new HashSet<String>();
	for (EQTL u : eqtls[setnames.length - 1]) {
	    uniqueprobes.add(u.getProbe());
	}
	setsizes[setnames.length - 1] = uniqueprobes.size();
	etf.close();



	double[][] overlaps = new double[setnames.length][setnames.length];
	for (int i = 0; i < setnames.length - 1; i++) {
	    EQTL[] a = eqtls[i];

	    HashSet<String> uniqueProbesA = new HashSet<String>();
	    for (EQTL e : a) {
		uniqueProbesA.add(e.getProbe());
	    }

	    String[] probeList = uniqueProbesA.toArray(new String[0]);

	    for (int j = i + 1; j < setnames.length; j++) {
		EQTL[] b = eqtls[j];
		HashSet<String> uniqueProbesB = new HashSet<String>();
		for (EQTL e : b) {
		    uniqueProbesB.add(e.getProbe());
		}


		int overlap = 0;
		for(String s: probeList){
		    if(uniqueProbesB.contains(s)){
			overlap++;
		    }
		}

		overlaps[i][j] = overlap;
	    }
	}


	System.out.println("");
	System.out.println("");
	System.out.println("");
	System.out.println("");

	MultiVenn m = new MultiVenn();
	m.plot(setnames, setsizes, overlaps);
	System.out.println("drawing");
	m.draw(origOuputDir + "/multivenn.png");

    }
}
