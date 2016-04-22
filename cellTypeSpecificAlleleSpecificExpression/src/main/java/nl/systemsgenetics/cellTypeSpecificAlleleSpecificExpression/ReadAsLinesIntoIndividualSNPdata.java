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
import java.util.ArrayList;
import org.jdom.IllegalDataException;

/**
 *
 * @author adriaan
 * 
 * This class will be used to read all AS lines. 
 * It will be used in this way:
 * it will open up the files specified in some folder, and then it will 
 * have a method readNextlines that will output an arrayList of IndividualSnpData
 * 
 * WARNING! this class opens up as many files as there are lines in the files.
 * Make sure to not call this class too often.
 * 
 */
public class ReadAsLinesIntoIndividualSNPdata {
    
    //Contains the filename with locations.
    String originalFileName;
    //contains the locations as specified in originalFileName file;
    ArrayList<String> asFileNames;
    
    //For later reading.
    private final ArrayList<File> all_files;
    private final ArrayList<BufferedReader> all_file_readers;
    
    public ReadAsLinesIntoIndividualSNPdata(String allFiles) throws IOException{
        
        originalFileName = allFiles;
        
        all_files = initialize_all_files(originalFileName);
        
        //fill up asFilenames
        asFileNames = new ArrayList<String>();
        for(File iFile :all_files){
            asFileNames.add(iFile.getName());
        }
        
        //Start the fileReaders
        all_file_readers = initFileReaders(all_files);
    
    }
    
    public ArrayList<IndividualSnpData> getIndividualsFromNextLine() throws IOException, IllegalDataException, IllegalArgumentException{
        
        ArrayList<String> allData = read_data_line();
         
        ArrayList<IndividualSnpData> theseSnps = UtilityMethods.parse_AS_lines(allData, all_files);
        
        if(theseSnps.isEmpty()){
            return theseSnps;
        }
        
        
        //Now do some checks to make sure everything is correct.
        String refSnpName    = theseSnps.get(0).getSnpName();
        String refChromosome = theseSnps.get(0).getChromosome();
        String refPosition   = theseSnps.get(0).getPosition();    
        for(IndividualSnpData iSnp: theseSnps){
            
            if(!refSnpName.equals(iSnp.getSnpName())){
                throw new IllegalDataException("refSnp was not the same across AS files");
            }
            if(!refChromosome.equals(iSnp.getChromosome())){
                throw new IllegalDataException("snp chromosome was not the same across AS files");
            }
            if(!refPosition.equals(iSnp.getPosition())){
                throw new IllegalDataException("snp position was not the same across AS files");
            }
        
        }
        
        //if no exceptions are thrown then everything is correct and we move on.
        return theseSnps;
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
            //I tried to change buffersize, but this did not seem to work well for optimization purposes
            BufferedReader tempBuff = new BufferedReader(tempReader); 
            fileList.add(tempBuff);

        }
        return fileList;

    }

    private ArrayList<String> read_data_line() throws IOException, IllegalDataException {

        ArrayList<String> all_lines = new ArrayList<String>();
        
        //This for loop is also implicitly used to make sure there are no (partially) empty files. 
        for (BufferedReader iFile : all_file_readers) {

            String line = iFile.readLine();

            if (line != null) {
                all_lines.add(line);
            }
            
        }
        
        if(all_lines.isEmpty()){
            return new ArrayList<String>();
        } else if(all_lines.size() != all_file_readers.size()){
            throw new IllegalDataException("The AS files are of different length");
        }
        
        
        return all_lines;
    }


}
