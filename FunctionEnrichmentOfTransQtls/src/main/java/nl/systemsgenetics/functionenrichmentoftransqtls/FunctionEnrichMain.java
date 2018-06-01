package nl.systemsgenetics.functionenrichmentoftransqtls;

import java.io.IOException;

public class FunctionEnrichMain {
	public static void main(String[] args) {
		ConvertToChiSq c = new ConvertToChiSq();
		
		if (args.length < 1) {
			System.out.println("Usage: (convertchisq|correlatesum)");
		} else {
			if (args[0].equals("convertchisq")) {
				if (args.length < 5) {
					System.out.println("Usage: convertchisq traitfile zscoreloc outputloc nrperm");
				} else {
					try {
						c.run(args[1], args[2], args[3], Integer.parseInt(args[4]));
					} catch (IOException e) {
						e.printStackTrace();
					}
				}
			} else if (args[0].equals("correlatesum")) {
			
			}
		}
	}
}
