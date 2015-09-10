/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package nl.systemsgenetics.cell_type_specific_ase;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileNotFoundException;
import java.io.FileReader;
import java.io.IOException;
import java.io.PrintWriter;
import java.io.UnsupportedEncodingException;
import java.util.ArrayList;

/*
 *
 * @author Adriaan
 */

public class BinomialEntry {

    public static void BinomialEntry(String asLocations) throws FileNotFoundException, UnsupportedEncodingException, IOException {

        PrintWriter out_writer;
        out_writer = new PrintWriter("/media/fast/GENETICA/implementInJava/CEUASCOUNTS/beta_binomial_output.txt", "UTF-8");

        // filenames_file, is a file containing per line a filename, which to load.
        String filenames_file = asLocations;

        //open all the files we want to open.
        ArrayList<File> all_files;
        all_files = initialize_all_files(filenames_file);
        ArrayList<BufferedReader> all_file_readers;
        all_file_readers = initFileReaders(all_files);

        int i = 0;

        while (true) {

            //read one line from all the files.
            ArrayList<String> all_line_data = read_data_line(all_file_readers);

            //step out of the loop when there is no data
            //maybe there is a more elegant way to do this, but I'm not sure.
            if (all_line_data.isEmpty()) {
                break;
            }

            //parse the lines
            ArrayList<IndividualSnpData> all_snp_data;
            all_snp_data = parse_lines(all_line_data, all_files);

            //do the binomial test:
            BinomialTest results = new BinomialTest(all_snp_data);

            // Write the results to the out_file.
            if (results.isTestPerformed()) {
                out_writer.println(results.writeTestStatistics(true));
                i++;
            }

        }

        out_writer.close();

        System.out.println("Finished " + Integer.toString(i) + " tests, now closing.");

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

    private static ArrayList<IndividualSnpData> parse_lines(ArrayList<String> all_line_data, ArrayList<File> all_files) {

        ArrayList<IndividualSnpData> all_individuals;
        all_individuals = new ArrayList<IndividualSnpData>();
        int i = 0;

        for (String iLine : all_line_data) {

            all_individuals.add(new IndividualSnpData(all_files.get(i).getAbsolutePath(), iLine));

            i++;
        }

        return all_individuals;
    }

}
