package eqtlmappingpipeline.metaqtl4;

import org.apache.log4j.Logger;

/**
 * An awesome MetaQTL4
 *
 */
public class MetaQTL4ConsoleUI {

    private static final Logger logger = Logger.getLogger(MetaQTL4ConsoleUI.class);

    public MetaQTL4ConsoleUI(String[] args) {
        if (args.length < 1) {
            System.err.println("Use settings file");
            System.exit(-1);
        }

        String settings = null;

        for (int i = 0; i < args.length; i++) {
            String arg = args[i];
            String val = null;

            if (i + 1 < args.length) {
                val = args[i + 1];
            }

            if (arg.equals("--settings")) {
                settings = val;
            }
        }

        if (settings == null) {
            System.err.println("use --settings settings.xml");
        } else {
            try {
                MetaQTL4Settings mqtl4settings = new MetaQTL4Settings(settings, null, null);
                SingleDatasetAnalysis q = new SingleDatasetAnalysis(mqtl4settings);
            } catch (Exception e) {
                e.printStackTrace();
            }
        }
    }

}
