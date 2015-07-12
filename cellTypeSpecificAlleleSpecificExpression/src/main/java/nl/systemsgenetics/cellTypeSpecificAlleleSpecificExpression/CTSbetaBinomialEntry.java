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
import java.util.ArrayList;

/**
 *
 * @author adriaan
 */
class CTSbetaBinomialEntry {

    public CTSbetaBinomialEntry(String asLocations, String phenoTypeLocation, String outputLocation, int minReads, int minHets) throws IOException {
        
        
        // BETABINOMIAL TEST.
        // CURRENTLY DETERMINES THE OVERDISPERSION 
        // IN THE SAME WAY AS THE BETA BINOMIAL.
        
        //PART 1: read all individuals names from the files.
        readFileIntoStringArrayList DetermineAllFiles = new readFileIntoStringArrayList(asLocations);
        ArrayList<String> allFiles = DetermineAllFiles.getLines();
        
        //PART 2: determine the per sample overdispersion in the file.
        
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
        
        
        //PART 3: determine the cell proportions per sample, provided in the sample file.
        
        
        // filenames_file, is a file containing per line a filename, which to load.
        String filenames_file = asLocations;

        /*
            Determine cell proportion per sample
         */
         
        readFileIntoStringArrayList phenoReadObject = new readFileIntoStringArrayList(phenoTypeLocation);
        
        ArrayList<String> phenoString = phenoReadObject.getLines();
        
        /*
            Right now just assuming this is a file that is 
            ordered in the same way as the asLocation file.
            With per line the cell proportion that we can determine.
        */
        i = 0;
        double[] cellProp = new double[phenoString.size()];
        for(String samplePheno : phenoString){
            cellProp[i] = Double.parseDouble(samplePheno);
            i++;
        }
        
        
        //PART 4. Read the as files one line at a time and 
        //      determine the Cell type specific tests :).
        
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
            ArrayList<CTSdispersedIndividualSnpData> all_snp_data;
            all_snp_data = parse_lines(all_line_data, all_files,cellProp,  dispersionVals );

            //do the beta binomial test:
            CTSBetaBinomialTest results = new CTSBetaBinomialTest(all_snp_data, minReads, minHets);

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

    private static ArrayList<CTSdispersedIndividualSnpData> parse_lines(ArrayList<String> all_line_data, ArrayList<File> all_files, double[] cellProp ,double[] dispersion) {

        ArrayList<CTSdispersedIndividualSnpData> all_individuals;
        all_individuals = new ArrayList<CTSdispersedIndividualSnpData>();
        int i = 0;

        for (String iLine : all_line_data) {

            all_individuals.add(new CTSdispersedIndividualSnpData(all_files.get(i).getAbsolutePath(), iLine, cellProp[i], dispersion[i] ));

            i++;
        }

        return all_individuals;
    }
    
    
}
