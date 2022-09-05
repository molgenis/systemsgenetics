package mbqtl;

public class VariantQCObj {
	double maf = 0;
	double cr = 0;
	double hwep = 0;
	boolean passqc = false;

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
