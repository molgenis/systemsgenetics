package nl.systemsgenetics.functionenrichmentoftransqtls;

import java.io.IOException;
import org.apache.commons.lang3.ArrayUtils;

public class FunctionEnrichMain {
	public static void main(String[] args) throws IOException {
		
		
		if (args.length < 1) {
			System.out.println("Usage: (convertchisq|correlatesum)");
		} else {
			if (args[0].equals("convertchisq")) {
				if (args.length < 5) {
					System.out.println("Usage: convertchisq traitfile zscoreloc outputloc nrperm");
				} else {
					try {
						ConvertToChiSq c = new ConvertToChiSq();
						c.run(args[1], args[2], args[3], Integer.parseInt(args[4]));
					} catch (IOException e) {
						e.printStackTrace();
					}
				}
			} else if (args[0].equals("correlatesum")) {
				if (args.length < 5) {
					System.out.println("Usage: correlatesum pathwayMatrix listSignifcantPathways sumChi2Matrix outputMatrix");
				} else {
					CorrelateSumChi2ToPathways.main(ArrayUtils.remove(args,0));
				}
			}
		}
	}
}
