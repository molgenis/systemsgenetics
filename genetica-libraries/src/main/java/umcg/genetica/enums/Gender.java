package umcg.genetica.enums;

/**
 * Created by hwestra on 5/27/16.
 */
public enum Gender {

	MALE(0, "Male"),
	FEMALE(1, "Female"),
	UNKNOWN(-1, "Unknown");

	private final int number;
	private final String name;

	private Gender(int num, String name) {
		this.number = num;
		this.name = name;
	}

	public static Gender parseStatus(String statusStr) {
		statusStr = statusStr.toLowerCase().trim();

		if (statusStr.equals("male") || statusStr.equals("1")) {
			return Gender.MALE;
		} else if (statusStr.equals("female") || statusStr.equals("2")) {
			return Gender.FEMALE;
		} else {
			return Gender.UNKNOWN;
		}
	}

	public String getName() {
		return name;
	}

	public int getNumber() {
		return number;
	}

	public int compare(Gender other) {
		if (this.equals(other)) {
			return 0;
		} else if (this.number > other.number) {
			return 1;
		} else {
			return 0;
		}
	}

	public boolean equals(Gender other) {
		return this.number == other.number;
	}

	public String toString() {
		return name;
	}
}
