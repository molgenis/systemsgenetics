/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package umcg.genetica.io.trityper.bin;

import java.io.DataInputStream;
import java.io.DataOutputStream;
import java.io.EOFException;
import java.io.File;
import java.io.FileInputStream;
import java.io.FileOutputStream;
import java.io.IOException;
import java.util.ArrayList;
import umcg.genetica.io.trityper.TriTyperGeneticalGenomicsDataset;

/**
 *
 * @author harm-jan
 */
public class BinaryResultProbeSummary {

    private DataOutputStream out;
    private int counter = 0;
    private DataInputStream in;
    public static boolean W = true;
    public static boolean R = false;

    public BinaryResultProbeSummary(String filename, boolean W) throws IOException {
	if (W) {
	    out = new DataOutputStream(new FileOutputStream(new File(filename + ".ProbeSummary.dat")));
	} else {
	    in = new DataInputStream(new FileInputStream(new File(filename)));
	}
    }

    public void close() throws IOException {
	if (out != null) {
	    out.flush();
	    out.close();
	}
    }

    public void write(String[] probeList, Integer[][] probetranslation, TriTyperGeneticalGenomicsDataset[] gg) throws IOException {
	for (int p = 0; p < probeList.length; p++) {
	    for (int d = 0; d < gg.length; d++) {
		if (probetranslation[d][p] != null) {
		    int probeId = probetranslation[d][p];
		    int start = gg[d].getExpressionData().getChrStart()[probeId];
		    int stop = gg[d].getExpressionData().getChrStart()[probeId];
		    int midpoint = (int) Math.floor((double) (start + stop) / 2);

		    out.writeInt(p);
		    out.writeUTF(gg[d].getExpressionData().getProbes()[probeId]);

		    out.writeByte(gg[d].getExpressionData().getChr()[probeId]);
		    out.writeInt(midpoint);

		    if (gg[d].getExpressionData().getAnnotation()[probeId] == null) {
			out.writeUTF("-");
		    } else {
			out.writeUTF(gg[d].getExpressionData().getAnnotation()[probeId]);
		    }


//                    System.out.println(p+"\t"+gg[d].getExpressionData().getProbes()[probeId]+"\t"+gg[d].getExpressionData().getChr()[probeId]+"\t"+midpoint+"\t"+gg[d].getExpressionData().getAnnotation()[probeId]);
		    // index probe chr midpoint annotation
		    break;
		}
	    }
	}
    }

    public void write(String[] probeNames) throws IOException {
	for (int p = 0; p < probeNames.length; p++) {
	    out.writeInt(p);
	    out.writeUTF(probeNames[p]);
	    out.writeByte(Byte.MAX_VALUE);
	    out.writeInt(Integer.MAX_VALUE);
	    out.writeUTF("-");
	}
    }

    public BinaryResultProbe[] readAllProbes() throws IOException {
	ArrayList<BinaryResultProbe> probes = new ArrayList<BinaryResultProbe>();
	BinaryResultProbe probe = readNextProbe();

	int ct = 0;
	while (probe != null) {
	    probes.add(probe);
	    probe = readNextProbe();
	}

	BinaryResultProbe[] probelist = new BinaryResultProbe[probes.size()];
	for (int p = 0; p < probelist.length; p++) {
	    probelist[p] = probes.get(p);
	}
	return probelist;
    }

    public BinaryResultProbe readNextProbe() throws IOException {
	BinaryResultProbe p = null;
	try {
	    p = new BinaryResultProbe();
	    p.setId(in.readInt());
	    p.setName(in.readUTF());
	    p.setChr(in.readByte());
	    p.setMidpoint(in.readInt());
	    p.setAnnotation(in.readUTF());
	} catch (EOFException e) {
	    return null;
	}
	return p;
    }
}
