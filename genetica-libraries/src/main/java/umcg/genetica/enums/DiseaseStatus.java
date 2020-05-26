package umcg.genetica.enums;

/**
 * Created by hwestra on 5/27/16.
 */
public enum DiseaseStatus {

	CONTROL(0, "Control", 1),
	CASE(1, "Case", 2),
	UNKNOWN(-1, "Unknown", -1);

	private final int number;
	private final String name;
	private final int famnumber;

	private DiseaseStatus(int num, String name, int famnumber) {
		this.number = num;
		this.name = name;
		this.famnumber = famnumber;
	}

	public static DiseaseStatus parseStatus(String statusStr) {
		statusStr = statusStr.toLowerCase().trim();

		if (statusStr.equals("control") || statusStr.equals("" + CONTROL.famnumber)) {
			return DiseaseStatus.CONTROL;
		} else if (statusStr.equals("case") || statusStr.equals("" + CASE.famnumber)) {
			return DiseaseStatus.CASE;
		} else {
			return DiseaseStatus.UNKNOWN;
		}
	}

	public String getName() {
		return name;
	}

	public int getNumber() {
		return number;
	}

	public int compare(DiseaseStatus other) {
		if (this.equals(other)) {
			return 0;
		} else if (this.number > other.number) {
			return 1;
		} else {
			return 0;
		}
	}

	public boolean equals(DiseaseStatus other) {
		return this.number == other.number;
	}

	public String toString() {
		return name;
	}
}
