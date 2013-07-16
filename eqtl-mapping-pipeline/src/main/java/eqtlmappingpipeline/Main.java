/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package eqtlmappingpipeline;

import eqtlmappingpipeline.gui.EQTLMappingPipelineConsole;

/**
 *
 * @author harmjan
 */
public class Main {

	public static final String VERSION = Main.class.getPackage().getImplementationVersion();

    /**
     * @param args the command line arguments
     */
    public static void main(String[] args) {

        EQTLMappingPipelineConsole app = new EQTLMappingPipelineConsole();
        app.main(args);


        System.out.println("Have a nice day :)");
        System.exit(0);
        
    }
}
