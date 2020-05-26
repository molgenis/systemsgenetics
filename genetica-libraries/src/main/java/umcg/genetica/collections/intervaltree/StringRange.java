package umcg.genetica.collections.intervaltree;

/**
 *
 * @author Patrick Deelen
 */
public class StringRange implements Range{

	private final int start;
	private final int end;
	private final String value;

	public StringRange(int start, int end, String value) {
		this.start = start;
		this.end = end;
		this.value = value;
	}
	
	@Override
	public int getStart() {
		return start;
	}

	@Override
	public int getEnd() {
		return end;
	}

	public String getValue() {
		return value;
	}

	@Override
	public String toString() {
		return value;
	}

	@Override
	public int hashCode() {
		int hash = 3;
		hash = 97 * hash + this.start;
		hash = 97 * hash + this.end;
		hash = 97 * hash + (this.value != null ? this.value.hashCode() : 0);
		return hash;
	}

	@Override
	public boolean equals(Object obj) {
		if (obj == null) {
			return false;
		}
		if (getClass() != obj.getClass()) {
			return false;
		}
		final StringRange other = (StringRange) obj;
		if (this.start != other.start) {
			return false;
		}
		if (this.end != other.end) {
			return false;
		}
		if ((this.value == null) ? (other.value != null) : !this.value.equals(other.value)) {
			return false;
		}
		return true;
	}

	
}
