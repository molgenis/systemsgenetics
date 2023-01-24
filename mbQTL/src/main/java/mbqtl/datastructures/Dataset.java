package mbqtl.datastructures;

import java.util.ArrayList;

public class Dataset {
	public String name;
	public ArrayList<Integer> RNAIds;
	public ArrayList<Integer> DNAIds;
	public int[] expressionIds;
	public int[] genotypeIds;

	public void append(Integer dnaSampleId, Integer rnaSampleId) {
		DNAIds.add(dnaSampleId);
		RNAIds.add(rnaSampleId);
	}

	public void toStr() {
		System.out.println(this.name + "\tSize: " + DNAIds.size());
	}

	public double[] select(double[] data, int[] ids) {
		double[] output = new double[ids.length];
		for (int i = 0; i < ids.length; i++) {
			output[i] = data[ids[i]];
		}
		return output;
	}

	public double[] select(double[][] data, int[] ids) {
		double[] output = new double[ids.length];
		for (int i = 0; i < ids.length; i++) {
			output[i] = data[0][ids[i]];
		}
		return output;
	}

	public double[] select(byte[] data, int[] ids) {
		double[] output = new double[expressionIds.length];
		for (int i = 0; i < expressionIds.length; i++) {
			output[i] = data[expressionIds[i]];
		}
		return output;
	}
}
