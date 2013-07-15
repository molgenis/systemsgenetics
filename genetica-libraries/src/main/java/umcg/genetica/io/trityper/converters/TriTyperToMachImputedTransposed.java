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
public class TriTyperToMachImputedTransposed {

    public static void main(String[] args){
    
        try{
            TriTyperToMachImputedTransposed.convert("/Data/GeneticalGenomicsDatasets/BloodH8v2ImputeTriTyper/", "/Data/GeneticalGenomicsDatasets/BloodH8v2ImputeTriTyper/MachImputedTransposed.txt", "/Users/harmjan/Downloads/cis_SNPs.txt");
        } catch (IOException e){
            
        }
        
    }
    
    public static void convert(String trityperdataset, String output, String snplists) throws IOException {

        if (trityperdataset == null) {
            throw new IllegalArgumentException("Trityper folder should be set");
        }

        if (output == null) {
            throw new IllegalArgumentException("Output folder should be set");
        }
        
        String parentDir = Gpio.getParentDir(output);
        
        Gpio.createDir(parentDir);


        // Convert TriTyper Data to MACH format...
        // SNP chr pos al1 al2 freq1 qual dosage1 dosage2
        TriTyperGenotypeData tt = new TriTyperGenotypeData();
        tt.load(trityperdataset);

        HashSet<String> selectedSNPs = null;

        if (snplists != null) {
            selectedSNPs = new HashSet<String>();
            TextFile tf = new TextFile(snplists, TextFile.R);
            selectedSNPs.addAll(tf.readAsArrayList());
            tf.close();
            System.out.println("About to export "+selectedSNPs.size()+" SNPs from "+snplists);
        }

        output+="ExportedSNPs-MachImputedTransposed.txt";

        TextFile outfile = new TextFile(output, TextFile.W);

        SNPLoader l = tt.createSNPLoader();
        String[] snps = tt.getSNPs();
        String[] individuals = tt.getIndividuals();

        // write the header

        String header = "SNP\tCHR\tPOS\tAL1\tAL2\tFreq\tQual";
        for (int ind = 0; ind < individuals.length; ind++) {
            if (tt.getIsIncluded()[ind]) {
                header += "\t" + individuals[ind];
            }
        }

        outfile.writeln(header);
        // iterate the sNPs.
        System.out.println("Exporting: ");
        int ctr = 0;
        for (int s = 0; s < snps.length; s++) {
            if (selectedSNPs == null || selectedSNPs.contains(snps[s])) {
                StringBuilder sb = new StringBuilder();
                sb.append(snps[s]).append("\t").append(tt.getChr(s)).append("\t").append(tt.getChrPos(s)).append("\t");
                SNP snpObj = tt.getSNPObject(s);
                l.loadGenotypes(snpObj);
                l.loadDosage(snpObj);

                sb.append(BaseAnnot.toString(snpObj.getAlleles()[0])).append("\t");
                sb.append(BaseAnnot.toString(snpObj.getAlleles()[1])).append("\t");

                if (snpObj.getMinorAllele() == (snpObj.getAlleles()[0])) {
                    sb.append(snpObj.getMAF()).append("\t");
                } else {
                    sb.append(1d - snpObj.getMAF()).append("\t");
                }


                sb.append("1.0");
                for (int ind = 0; ind < individuals.length; ind++) {
                    if (tt.getIsIncluded()[ind]) {
                        sb.append("\t").append(snpObj.getDosageValues()[ind]);
                    }
                }
                outfile.writeln(sb.toString());
                snpObj.clearGenotypes();
                ctr++;
                if(ctr%100 == 0){
                    System.out.println(ctr+"\texported");
                }
            }
        }
        outfile.close();

        System.out.println("Done exporting "+ctr+" SNPs ");


    }
}
