package eqtlmappingpipeline.metaqtl4;


import org.apache.log4j.Logger;

/**
 * An awesome MetaQTL4
 *
 */
public class Main {

    private static final Logger logger = Logger.getLogger(Main.class);
    
    public static void main(final String[] args) throws Exception {
    
        if(args.length < 1){
            System.err.println("Use settings file");
            System.exit(-1);
        }
        try{
            MetaQTL4Settings settings = new MetaQTL4Settings(args[0], null, null);
            MetaQTL4 q = new MetaQTL4(settings);
        } catch (Exception e){
            e.printStackTrace();
        }
    }

    
}
