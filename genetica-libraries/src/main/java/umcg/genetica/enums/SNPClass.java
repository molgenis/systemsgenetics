package umcg.genetica.enums;

/**
 * Created by hwestra on 10/24/16.
 */
public enum SNPClass {

	NONCODING("Noncoding", 0),
	EXONIC("Exonic", 1),
	INTRONIC("Intronic", 2),
	PROMOTER("Promoter", 3),
	INDEL("Indel", 4),
	UTR("UTR", 5);

	private final String name;
	private final int number;


	SNPClass(String name, int num) {
		this.name = name;
		this.number = num;
	}

	public int getNumber() {
		return number;
	}

	public String getName() {
		return name;
	}
}
