/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package eqtlmappingpipeline;

import eqtlmappingpipeline.gui.EQTLMappingPipelineConsole;

import java.text.DateFormat;
import java.text.SimpleDateFormat;
import java.util.ResourceBundle;

/**
 * @author harmjan
 */
public class Main {
	
	public static final String VERSION = ResourceBundle.getBundle("version").getString("application.version");
	public static final DateFormat DATE_TIME_FORMAT = new SimpleDateFormat("yyyy-MM-dd HH:mm:ss");
	
	/**
	 * @param args the command line arguments
	 */
	public static void main(String... args) {
		try {
			EQTLMappingPipelineConsole app = new EQTLMappingPipelineConsole();
			app.main(args);
		} catch (Exception e) {
			e.printStackTrace();
		}
		
	}
}
