package betaqtl;

public class VariantQCObj {
	double maf;
	double cr;
	double hwep;
	boolean passqc;

	@Override
	public String toString() {
		return "QCObj{" +
				"maf=" + maf +
				", cr=" + cr +
				", hwep=" + hwep +
				", passqc=" + passqc +
				'}';
	}
}
