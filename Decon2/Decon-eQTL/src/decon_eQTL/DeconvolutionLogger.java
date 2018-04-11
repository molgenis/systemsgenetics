package decon_eQTL;

import java.io.File;
import java.io.IOException;
import java.nio.file.Files;
import java.text.DateFormat;
import java.text.SimpleDateFormat;
import java.util.Date;
import java.util.logging.ConsoleHandler;
import java.util.logging.FileHandler;
import java.util.logging.Level;
import java.util.logging.LogManager;
import java.util.logging.Logger;

public class DeconvolutionLogger {
	static private FileHandler outfilePath;
	protected final static Logger log = Logger.getLogger(Logger.GLOBAL_LOGGER_NAME);
	static public void setup(String outputDir, Boolean noConsole) throws IOException {
		// get the global logger to configure it
		LogManager.getLogManager().reset();
		log.setLevel(Level.INFO);
		DateFormat dateFormat = new SimpleDateFormat("yyyyMMdd");
		Date date = new Date();
		File file = new File(outputDir+"/DeconvolutionLog_"+dateFormat.format(date)+".txt");
		Files.deleteIfExists(file.toPath());
		setOutfilePath(new FileHandler(outputDir+"/DeconvolutionLog_"+dateFormat.format(date)+".txt"));
		CustomRecordFormatter customFormatter = new CustomRecordFormatter();
		ConsoleHandler consoleHandler = new ConsoleHandler();
		consoleHandler.setFormatter(customFormatter);
		outfilePath.setFormatter(customFormatter);
		log.setUseParentHandlers(false);
		if (!noConsole){
			log.addHandler(consoleHandler);
		}
		log.addHandler(outfilePath);

	}
	public static FileHandler getOutfilePath() {
		return outfilePath;
	}
	public static void setOutfilePath(FileHandler outfilePath) {
		DeconvolutionLogger.outfilePath = outfilePath;
	}
}


