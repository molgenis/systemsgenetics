//
// Source code recreated from a .class file by IntelliJ IDEA
// (powered by FernFlower decompiler)
//

package mbqtl.datastructures;

import java.util.ArrayList;
import java.util.Objects;

import umcg.genetica.util.Primitives;

public class Dataset implements Comparable<Dataset> {
	private String name;
	private int[] genotypeIds;
	private int[] expressionIds;
	private ArrayList<Integer> DNAIds = new ArrayList();
	private ArrayList<Integer> RNAIds = new ArrayList();

	public Dataset() {
	}

	public String getName() {
		return name;
	}

	public int[] getGenotypeIds() {
		if (genotypeIds == null) {
			this.genotypeIds = Primitives.toPrimitiveArr(this.DNAIds);
		}
		return genotypeIds;
	}

	public int[] getExpressionIds() {
		if (expressionIds == null) {
			this.expressionIds = Primitives.toPrimitiveArr(this.RNAIds);
		}
		return expressionIds;
	}

	public ArrayList<Integer> getDNAIds() {
		return DNAIds;
	}

	public ArrayList<Integer> getRNAIds() {
		return RNAIds;
	}

	public void append(Integer dnaId, Integer rnaId) {
		this.DNAIds.add(dnaId);
		this.RNAIds.add(rnaId);
	}

	public void toArr() {
		this.genotypeIds = Primitives.toPrimitiveArr(this.DNAIds);
		this.expressionIds = Primitives.toPrimitiveArr(this.RNAIds);
		this.DNAIds = null;
		this.RNAIds = null;
	}

	public double[] select(double[] vals, int[] ids) {
		double[] data = new double[ids.length];

		for (int g = 0; g < ids.length; ++g) {
			data[g] = vals[ids[g]];
		}

		return data;
	}

	public double[] select(byte[] vals, int[] ids) {
		double[] data = new double[ids.length];

		for (int g = 0; g < ids.length; ++g) {
			data[g] = (double) vals[ids[g]];
		}

		return data;
	}

	public double[] select(double[][] vals, int[] ids) {
		double[] data = new double[ids.length];

		for (int g = 0; g < ids.length; ++g) {
			data[g] = vals[ids[g]][0];
		}

		return data;
	}

	public boolean equals(Object o) {
		if (this == o) {
			return true;
		} else if (o != null && this.getClass() == o.getClass()) {
			Dataset dataset = (Dataset) o;
			return Objects.equals(this.name, dataset.name);
		} else {
			return false;
		}
	}

	public int hashCode() {
		return Objects.hash(new Object[]{this.name});
	}

	public int compareTo(Dataset dataset) {
		return this.equals(dataset) ? 0 : this.name.compareTo(dataset.name);
	}

	public int size() {
		return RNAIds.size();
	}

	public boolean isEmpty() {
		return RNAIds.isEmpty();
	}

	public void setName(String datasetName) {
		this.name = datasetName;
	}
}
