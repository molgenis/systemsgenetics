package nl.umcg.westrah.binarymetaanalyzer.westrah.binarymetaanalyzer.posthoc;

import umcg.genetica.io.text.TextFile;

import java.io.IOException;
import java.util.ArrayList;

public class SettingsFileCreator {
	
	
	public static void main(String[] args) {
		
		String directoryfile = "C:\\Data\\tmp\\2017-12-27-Standardization\\directoryFile.txt";
		String outputdir = "/groups/umcg-wijmenga/tmp03/projects/eQTLGen/analysis/trans-eqtl/standardization/pcCorrected/v1012/output/";
		String averagingmethod = "mean";
		String nrpermutations = "" + 10;
		String settingloclocal = "C:\\Data\\tmp\\2017-12-27-Standardization\\settings\\";
		String settinglocserver = "/groups/umcg-wijmenga/tmp03/projects/eQTLGen/analysis/trans-eqtl/standardization/pcCorrected/v1012/settings/";
		String tool = "/groups/umcg-wijmenga/tmp03/projects/eQTLGen/tools/BinaryMetaAnalyzer-1.0.12-SNAPSHOT-jar-with-dependencies.jar";
		SettingsFileCreator s = new SettingsFileCreator();
		try {
//			s.createInternalMeta(directoryfile, outputdir, averagingmethod, nrpermutations, settingloclocal, settinglocserver, tool);
			
			String snpAnnotationfile = "/groups/umcg-bios/tmp03/users/umcg-hwestra/eqtl/HelpFiles/GiantSNPMappings_Filtered_on_Maf0.001_Indels.txt.gz";
			String probetranslation = "/groups/umcg-bios/tmp03/users/umcg-hwestra/eqtl/HelpFiles/ProbeAnnotation2_CorrectlyMapped_mapping.txt";
			directoryfile = "";
			String shellout = "";
			
			
			s.createBinaryMeta(10, snpAnnotationfile, 200000000, probetranslation, outputdir, 32, directoryfile, shellout);
			
		} catch (IOException e) {
			e.printStackTrace();
		}
		
	}
	
	public void createBinaryMeta(int nrpermutations, String snpannotation, int nrFinalEQTLs, String probetranslation, String outputdirectory, int nrthreads, String directoryfile, String shelloutfile) throws IOException {
		
		TextFile input = new TextFile(directoryfile, TextFile.R);
		String[] elems = input.readLineElems(TextFile.tab); // datasetname prefix location
		
		class Dataset {
			String dsname;
			String prefix;
			String location;
			String useannotation;
		}
		
		ArrayList<Dataset> combos = new ArrayList<>();
		while (elems != null) {
			
			String dsname = elems[0];
			String prefix = elems[1];
			String location = elems[2];
			String useannotation = elems[3];
			Dataset combo = new Dataset();
			combo.dsname = dsname;
			combo.prefix = prefix;
			combo.location = location;
			combo.useannotation = useannotation;
			combos.add(combo);
			elems = input.readLineElems(TextFile.tab);
		}
		input.close();
		
		String output = "<settings>\n" +
				"\t<defaults>\n" +
				"\t<permutations>" + nrpermutations + "</permutations>\n" +
				"\t\t<startpermutation>0</startpermutation>\n" +
				"\t\t<snpannotation>" + snpannotation + "</snpannotation>\n" +
				"\t\t<absolutezscore>false</absolutezscore>\n" +
				"\t\t<finalnreqtls>" + nrFinalEQTLs + "</finalnreqtls>\n" +
				"\t\t<fdrthreshold>0.05</fdrthreshold>\n" +
				"\t\t<cisprobedistance>1000000</cisprobedistance>\n" +
				"\t\t<transprobedistance>5000000</transprobedistance>\n" +
				"\t\t<includeprobeswithoutmapping>true</includeprobeswithoutmapping>\n" +
				"\t\t<includesnpswithoutmapping>true</includesnpswithoutmapping>\n" +
				"\t\t<makezscoreplot>false</makezscoreplot>\n" +
				"\t\t<makezscoretable>true</makezscoretable>\n" +
				"\t\t<probetranslationfile>" + probetranslation + "</probetranslationfile>\n" +
				"\t\t<output>" + outputdirectory + "</output>\n" +
				"\t\t<minimalnumberofdatasetsthatcontainprobe>1</minimalnumberofdatasetsthatcontainprobe>\n" +
				"\t\t<minimalnumberofdatasetsthatcontainsnp>1</minimalnumberofdatasetsthatcontainsnp>\n" +
				"\t\t<snpprobeselectsamplesizethreshold>-1</snpprobeselectsamplesizethreshold>\n" +
				"\t\t<runonlypermutation>-1</runonlypermutation>\n" +
				"\t\t<threads>" + nrthreads + "</threads>\n" +
				"\t\t<cis>false</cis>\n" +
				"\t\t<trans>true</trans>\n" +
				"\t\t<probeselection></probeselection>\n" +
				"\t\t<snpprobeselection/>\n" +
				"\t</defaults>";
		
		output += "\t<datasets>";
		
		for (Dataset d : combos) {
			output += "\t\t<dataset>\n" +
					"\t\t\t<name>" + d.dsname + "</name>\n" +
					"\t\t\t<prefix>" + d.prefix + "</prefix>\n" +
					"\t\t\t<location>" + d.location + "</location>\n" +
					"\t\t\t<expressionplatform>" + d.useannotation + "</expressionplatform>\n" +
					"\t\t</dataset>";
		}
		
		output += "\t</datasets>";
		
		output += "</settings>";
		
		TextFile outshell = new TextFile(shelloutfile, TextFile.W);
		outshell.writeln(output);
		outshell.close();
	}
	
	public void createInternalMeta(String directoryfile, String outputdir, String averagingMethod, String nrPermutations, String settingsloclocal, String settingslocserver, String tool) throws IOException {
		
		
		/*
		 
		 
		 */
		
		TextFile input = new TextFile(directoryfile, TextFile.R);
		String[] elems = input.readLineElems(TextFile.tab); // datasetname prefix location
		
		class Dataset {
			String dsname;
			String prefix;
			String location;
			String probetofeature;
		}
		
		ArrayList<Dataset> combos = new ArrayList<>();
		while (elems != null) {
			
			String dsname = elems[0];
			String prefix = elems[1];
			String location = elems[2];
			String probetofeature = elems[3];
			Dataset combo = new Dataset();
			combo.dsname = dsname;
			combo.prefix = prefix;
			combo.location = location;
			combo.probetofeature = probetofeature;
			combos.add(combo);
			elems = input.readLineElems(TextFile.tab);
		}
		input.close();
		
		if (averagingMethod == null) {
			averagingMethod = "mean";
		}
		
		int ctr = 1;
		for (Dataset combo : combos) {
			
			String outputloc = outputdir + combo.dsname + "/";
			String output = "<?xml version=\"1.0\" encoding=\"utf-8\" standalone=\"no\"?> \n" +
					"<settings>\n" +
					"\t<defaults>\n" +
					"\t\t<permutations>" + nrPermutations + "</permutations>\n" +
					"\t\t<startpermutations>0</startpermutations>\n" +
					"\t\t<output>" + outputloc + "</output>\n" +
					"\t\t<probeToFeature>" + combo.probetofeature + "</probeToFeature>\n" +
					"\t\t<zScoreMerge>" + averagingMethod + "</zScoreMerge>\n" +
					"\t\t<dataset>\n" +
					"\t\t\t<name>" + combo.dsname + "</name>\n" +
					"\t\t\t<prefix>" + combo.prefix + "</prefix>\n" +
					"\t\t\t<location>" + combo.location + "</location>\n" +
					"\t\t\t</dataset>\n" +
					"\t</defaults>\n" +
					"</settings>";
			TextFile tf = new TextFile(settingsloclocal + ctr + ".xml", TextFile.W);
			tf.writeln(output);
			tf.close();
			
			tf = new TextFile(settingsloclocal + ctr + ".sh", TextFile.W);
			String shellout = "#!/bin/bash\n" +
					"#SBATCH --job-name=ST-" + combo.dsname + "\n" +
					"#SBATCH --output=" + outputdir + "/logs/ST-" + combo.dsname + ".out\n" +
					"#SBATCH --error=" + outputdir + "/logs/ST-" + combo.dsname + ".err\n" +
					"#SBATCH --time=24:00:00\n" +
					"#SBATCH --constraint=tmp03\n" +
					"#SBATCH --cpus-per-task=12\n" +
					"#SBATCH --mem=6gb\n" +
					"#SBATCH --nodes=1\n" +
					"#SBATCH --open-mode=append\n" +
					"#SBATCH --export=NONE\n" +
					"#SBATCH --get-user-env=L\n" +
					"\n" +
					"" +
					"java -Xmx5g -jar  " + tool + " --internalmeta \\\n--settings " + settingslocserver + ctr + ".xml";
			tf.writeln(shellout);
			tf.close();
			
			ctr++;
		}
		
		TextFile tf = new TextFile(settingsloclocal + "arrayjob.sh", TextFile.W);
		
		String shellout = "#!/bin/bash\n" +
				"#SBATCH --job-name=ST\n" +
				"#SBATCH --array=1-" + (ctr - 1) + "\n" +
				"#SBATCH --output=" + outputdir + "/logs/ST-job%a.out\n" +
				"#SBATCH --error=" + outputdir + "/logs/ST-job%a.err\n" +
				"#SBATCH --time=96:00:00\n" +
				"#SBATCH --constraint=tmp03\n" +
				"#SBATCH --cpus-per-task=12\n" +
				"#SBATCH --mem=6gb\n" +
				"#SBATCH --nodes=1\n" +
				"#SBATCH --open-mode=append\n" +
				"#SBATCH --export=NONE\n" +
				"#SBATCH --get-user-env=L\n" +
				"\n" +
				"java -Xmx5g -jar  " + tool + " --internalmeta \\\n--settings " + settingslocserver + "$SLURM_ARRAY_TASK_ID.xml";
		tf.writeln(shellout);
		tf.close();
		
	}
}