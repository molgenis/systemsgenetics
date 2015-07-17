package nl.systemsgenetics.eqtlinteractionanalyser.eqtlinteractionanalyser;

import java.io.IOException;
import java.util.ArrayList;
import java.util.EnumMap;

/**
 *
 * @author Patrick Deelen
 */
public class SkippedInteractionTracker {

	public static enum Reason {
		SINGULAR, SHARED_QTL 
	}
	
	private final String covariate;
	EnumMap<Reason, ArrayList<String>> skipped;

	public SkippedInteractionTracker(String covariate) {
		this.covariate = covariate;
		skipped = new EnumMap<Reason, ArrayList<String>>(Reason.class);
		for(Reason r : Reason.values()){
			skipped.put(r, new ArrayList<String>());
		}
	}
	
	public void addSkipped(Reason r, String qtl){
		skipped.get(r).add(qtl);
	}

	public String getCovariate() {
		return covariate;
	}

	public EnumMap<Reason, ArrayList<String>> getSkipped() {
		return skipped;
	}
		
}
