/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package umcg.genetica.io.trityper.converters;

import java.io.IOException;
import java.util.HashSet;
import umcg.genetica.io.Gpio;
import umcg.genetica.io.text.TextFile;
import umcg.genetica.io.trityper.SNP;
import umcg.genetica.io.trityper.SNPLoader;
import umcg.genetica.io.trityper.TriTyperGenotypeData;
import umcg.genetica.io.trityper.util.BaseAnnot;

/**
 *
 * @author harmjan
 */
public class TriTyperToTabSeparated {

    public void convert(String trityperDir, String outDir, String snpfile) throws IOException {
        TextFile queryFile = new TextFile(snpfile, TextFile.R);
        HashSet<String> query = new HashSet<String>();
        query.addAll(queryFile.readAsArrayList());
        queryFile.close();

        TriTyperGenotypeData ds = new TriTyperGenotypeData();
        ds.load(trityperDir);

        if (outDir.endsWith("/")) {
            outDir += "/";
        }
        if (!Gpio.exists(outDir)) {
            Gpio.createDir(outDir);
        }

        SNPLoader loader = ds.createSNPLoader();

        for (String s : query) {
            TextFile out = new TextFile(outDir + s + ".txt", TextFile.W);

            int SNPId = ds.getSnpToSNPId().get(s);
            if (SNPId == -9) {
                out.writeln("SNP " + s + " not present in dataset");
            } else {
                SNP snpObj = ds.getSNPObject(SNPId);

                loader.loadGenotypes(snpObj);
                if (loader.hasDosageInformation()) {
                    loader.loadDosage(snpObj);
                }

                out.writeln("SNP: " + s + "\tHWE: " + snpObj.getHWEP() + "\tMAF: " + snpObj.getMAF() + "\tCR: " + snpObj.getCR() + "\tAlleles: " + BaseAnnot.toString(snpObj.getAlleles()[0]) + "/" + BaseAnnot.toString(snpObj.getAlleles()[1]) + "\tMinor: " + BaseAnnot.toString(snpObj.getMinorAllele()));

                String header = "Ind\tIsIncluded\tIsCase\tIsFemale\tGenotype";
                if(loader.hasDosageInformation()){
                    header+="\tDosage";
                }
                out.writeln(header);
                String[] inds = ds.getIndividuals();
                int nrInds = inds.length;

                for (int i = 0; i < nrInds; i++) {
                    String output = inds[i] +"\tisIncluded: "+ds.getIsIncluded()[i]+"\tisCase: "+ds.getIsCase()[i]+"\tisFemale:"+ds.getIsFemale()[i]+ "\t" + BaseAnnot.toString(snpObj.getAllele1()[i]) + BaseAnnot.toString(snpObj.getAllele2()[i]);
                    if (loader.hasDosageInformation()) {
                        output += "\t" + snpObj.getDosageValues()[i];
                    }
                    out.writeln(output);
                }
				snpObj.clearGenotypes();
            }

            out.close();
        }

        loader.close();
    }
}
