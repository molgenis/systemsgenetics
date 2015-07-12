/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package nl.systemsgenetics.cellTypeSpecificAlleleSpecificExpression;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileNotFoundException;
import java.io.FileReader;
import java.io.IOException;
import java.io.PrintWriter;
import java.io.UnsupportedEncodingException;
import java.util.ArrayList;

/**
 *
 * @author Adriaan van der Graaf 2015
 * 
 * This class will determine Maximum likelihood of Beta Binomial Over dispersion 
 * Then it will perform beta binomial likelihood ratio test per SNP.
 * 
 * Current uncertainty: 
 *                  I just take the likelihood of a certain a in null and compare 
 *                  this to freely estimated, don't know if this is correct.
 * 
 * Limitation to this is that it will read the as_files twice. 
 * Maybe there might be some better logic for this, but this would require quite 
 * some memory, so I'll read through them twice.
 * 
 */
public class BetaBinomEntry {
   
    
    public BetaBinomEntry(String asLocations, String outputLocation, int minReads, int minHets) throws FileNotFoundException, UnsupportedEncodingException, IOException{
        
        /*
         * PART 1: read all individual names
         */
        readFileIntoStringArrayList DetermineAllFiles = new readFileIntoStringArrayList(asLocations);
        ArrayList<String> allFiles = DetermineAllFiles.getLines();
        
        
        /*
         * PART 2: Determine all Overdispersion parameters for all samples.
         * TODO Implement some logic if there is already overdispersion data available.
         *      Perhaps implement this inside the betaBinomOverdispInSample side.
         *      Not necassary now, but later, yes.
         */
        
        ArrayList<betaBinomOverdispInSample>  dispersionParameters = new ArrayList<betaBinomOverdispInSample>();
        
        for(String sampleName : allFiles){
           
            dispersionParameters.add(new betaBinomOverdispInSample(sampleName));
            /*
             * TODO: determine some way to save the overdispersion sample.
             *       But we could also do it everytime you do this, perhaps based 
             *       on the relative time compared to the full calculation.
             */
        }
        double[] dispersionVals = new double[dispersionParameters.size()];  
        int i = 0;     
        
        for(betaBinomOverdispInSample sampleDispersion : dispersionParameters){
            dispersionVals[i] = sampleDispersion.getOverdispersion()[0];
            
            
            //do a check to make sure ordering is correct.
            if(!(sampleDispersion.getSampleName().equals(allFiles.get(i)))){
                System.out.println("ERROR! ordering is not correct filenames for overdispersion");
                System.out.println(sampleDispersion.getSampleName());
                System.out.println(allFiles.get(i));
            }
            
            i++;
        }
                
                
                
        /*
         * Part 3: Read the AS_files per line, and do the Binomial Test per SNP.
         *         
         */
        
        
        PrintWriter out_writer;
        out_writer = new PrintWriter(outputLocation, "UTF-8");
        
        //open all the files we want to open.
        ArrayList<File> all_files;
        all_files = initialize_all_files(asLocations);
        ArrayList<BufferedReader> all_file_readers;
        all_file_readers = initFileReaders(all_files);

        while (true) {

            //read one line from all the files.
            ArrayList<String> all_line_data = read_data_line(all_file_readers);

            //step out of the loop when there is no data
            //maybe there is a more elegant way to do this, but I'm not sure.
            if (all_line_data.isEmpty()) {
                break;
            }

            //parse the lines
            ArrayList<dispersedIndividualSnpData> all_snp_data;
            all_snp_data = parse_lines(all_line_data, all_files, dispersionVals);

            //do the beta binomial test:
            BetaBinomialTest results = new BetaBinomialTest(all_snp_data, minReads, minHets);

            // Write the results to the out_file.
            if (results.isTestPerformed()) {
                out_writer.println(results.writeTestStatistics(true));
                i++;
            }

        }
        
        out_writer.close();
        
    }


    private static ArrayList<File> initialize_all_files(String filenames_file) throws IOException {

        ArrayList<File> fileList = new ArrayList<File>();

        File file = new File(filenames_file);
        FileReader fr = new FileReader(file);
        BufferedReader br = new BufferedReader(fr);
        String line;

        while ((line = br.readLine()) != null) {

            File tempFile = new File(line);

            fileList.add(tempFile);

        }
        br.close();
        fr.close();
        return fileList;
    }

    private static ArrayList<BufferedReader> initFileReaders(ArrayList<File> all_files) throws FileNotFoundException {
        ArrayList<BufferedReader> fileList = new ArrayList<BufferedReader>();

        for (File i_file : all_files) {

            FileReader tempReader = new FileReader(i_file);
            BufferedReader tempBuff = new BufferedReader(tempReader);
            fileList.add(tempBuff);

        }
        return fileList;

    }

    private static ArrayList<String> read_data_line(ArrayList<BufferedReader> all_readers) throws IOException {

        ArrayList<String> all_lines = new ArrayList<String>();

        for (BufferedReader iFile : all_readers) {

            String line = iFile.readLine();

            if (line != null) {
                all_lines.add(line);
            } else {
                return new ArrayList<String>();

            }
        }

        return all_lines;
    }

    private static ArrayList<dispersedIndividualSnpData> parse_lines(ArrayList<String> all_line_data, ArrayList<File> all_files, double[] dispersion) {

        ArrayList<dispersedIndividualSnpData> all_individuals;
        all_individuals = new ArrayList<dispersedIndividualSnpData>();
        int i = 0;

        for (String iLine : all_line_data) {

            all_individuals.add(new dispersedIndividualSnpData(all_files.get(i).getAbsolutePath(), iLine, dispersion[i]));

            i++;
        }

        return all_individuals;
    }



}
