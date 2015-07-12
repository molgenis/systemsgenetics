/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */

package nl.systemsgenetics.cellTypeSpecificAlleleSpecificExpression;


import htsjdk.samtools.CigarElement;
import java.io.BufferedInputStream;
import java.io.DataInputStream;
import java.io.File;
import java.io.FileInputStream;
import java.io.FileNotFoundException;
import java.io.IOException;
import java.util.Arrays;

import java.util.HashMap;

import java.util.Set;
import org.molgenis.genotype.Alleles;
import org.molgenis.genotype.RandomAccessGenotypeData;
import org.molgenis.genotype.trityper.TriTyperGenotypeData;
import org.molgenis.genotype.variant.GeneticVariant;
import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.SAMRecordIterator;
import htsjdk.samtools.SamReader;
import htsjdk.samtools.SamReaderFactory;
import java.io.PrintWriter;
import java.util.ArrayList;
import org.apache.commons.io.FilenameUtils;



/**
 *
 * @author Adriaan van der Graaf
 */
public class readGenoAndAsFromIndividual {
    
    public static void readGenoAndAsFromIndividual(String loc_of_bam1, String genotype_loc, String coupling_location, String outputLocation) throws IOException, Exception{
        
        //Print ASREADS header

        System.out.println("---- Starting ASREADS for the following settings: ----");
        System.out.println("\t input bam:       " +  loc_of_bam1);
        System.out.println("\t genotype file:   " +  genotype_loc);
        System.out.println("\t coupling file:   " +  coupling_location);
        System.out.println("\t output location: " +  outputLocation);
        System.out.println("------------------------------------------------------");


        String path_of_trityper = genotype_loc;

        //parse command line arguments
        String loc_of_bam;
        loc_of_bam = loc_of_bam1; //"/media/fast/GENETICA/implementInJava/CEURNASEQ/ERR188047.MERGED.SORTED.bam";
        System.out.println("Location of bam file: ");
        System.out.println(loc_of_bam);
        
        if (!new File(loc_of_bam).exists())
        {
           System.out.println("Location of bam file is not an existing file. Exitting.");
           return;
        }else{
            System.out.println("Location of bam file is an existing file, will continue.");

        }
        
        
        RandomAccessGenotypeData dataset1;
        //open tri typer dataset
        try {
            dataset1 = new TriTyperGenotypeData(new File(path_of_trityper));
        } catch (IOException ex) {
            System.err.println("Error reading dataset 1: " + path_of_trityper);
            return;
        }
        
        //get the variants in the variantIdMAP
        HashMap<String, GeneticVariant> variantIdMap = dataset1.getVariantIdMap();
        Set<String> snpNames = variantIdMap.keySet();
        
        
        String[] individual_names;
        individual_names = dataset1.getSampleNames();
        
        //String path = "/gcc/groups/lld/tmp01/projects/bamFiles/";
        //sample_map contains all the individuals that are in the sample file.
        HashMap sample_map = convert_individual_names(individual_names, coupling_location);
        System.out.println("Sample names were loaded.");
        //System.out.println(sample_map.toString());

       
        //Twice beacuse my files have the .MERGED.sorted.bam suffix attached to them.
        String sample_name = FilenameUtils.getBaseName(FilenameUtils.getBaseName(FilenameUtils.getBaseName(loc_of_bam)));
        System.out.println("sample_name: " + sample_name);
        
        //sample name used for testing: "AC1C40ACXX-5-14";
        Object sample_idx = sample_map.get(sample_name);   
       
        if(sample_idx == null){
            System.out.println("Couldn't find the filename in the sample names. Quitting.");
            return;
        }

        int sample_index = Integer.parseInt(sample_idx.toString());
        
        System.out.println("sample_index: " + sample_index);
        
        //bam file path and filename
        String path_and_filename = loc_of_bam;
        File sample_file = new File(path_and_filename);

        SamReader bam_file = SamReaderFactory.makeDefault().open(sample_file);

        
        int number_of_samples;         
        number_of_samples = Arrays.asList(sample_map).size();
        
        System.out.println("Initialized for reading bam file");
        
        PrintWriter writer = new PrintWriter(outputLocation, "UTF-8");
        
        
        int i = 0;
        for(String i_snp : snpNames){
            
            //System.out.println(i_snp);
            
            GeneticVariant this_variant = variantIdMap.get(i_snp);
            
            String chromosome = this_variant.getSequenceName();
            String position = String.valueOf(this_variant.getStartPos());
            
//            System.out.println("chromosome: " + chromosome + " pos: " + position);
            
            // We only do analyses if we find a SNP and it is biallelic      
            if(this_variant.isSnp() & this_variant.isBiallelic()){
//                System.out.println(this_variant.getAllIds());
                                
                String row_of_table = get_allele_specific_overlap_at_snp(this_variant, 
                                                                        sample_index,
                                                                        chromosome,
                                                                        position,
                                                                        bam_file);
                 
                
                
                writer.println(chromosome + "\t" + position + "\t" + i_snp + "\t" + 
                               this_variant.getVariantAlleles().getAllelesAsChars()[0] + "\t" +
                               this_variant.getVariantAlleles().getAllelesAsChars()[1] + "\t" +
                               row_of_table + "\t" + 
                               Arrays.toString(this_variant.getSampleVariants().get(sample_index).getAllelesAsChars()) 
                            );
            }
         
            i++;
            if(i % 10000 == 0){

                System.out.println("Finished " + Integer.toString(i) + " SNPs");
            
            }
        
        }
        
        dataset1.close();
        writer.close();
    }
    
    
    public static String get_allele_specific_overlap_at_snp(GeneticVariant this_variant,
                                                            int sample_index, 
                                                            String chromosome,
                                                            String position, 
                                                            SamReader bam_file){
        
        int pos_int = Integer.parseInt(position);
        
        Alleles all_variants = this_variant.getVariantAlleles();
        Character ref_allele_char = all_variants.getAllelesAsChars()[0];
        String ref_allele = ref_allele_char.toString(); 
        //System.out.println("ref_allele: " + ref_allele);
        Character alt_allele_char = all_variants.getAllelesAsChars()[1];
        String alt_allele = alt_allele_char.toString(); 
        //System.out.println("alt_allele: " + alt_allele);
       
        
        
        int ref_overlap = 0;
        int alt_overlap = 0;
        int no__overlap = 0;
        
        

       
        // now determine per individual the sample variants.
        // I'm assuming the ordering is the same as the individual names created 
        // by the  getSampleNames() method. 
        // Otherwise the data will be nicely permuted, and I will have to convert some stuff.
        
        int i=0;

        String temp_string;
        
        int position_of_snp = Integer.parseInt(position);                       
        
        SAMRecordIterator all_reads_in_region;
        all_reads_in_region = bam_file.queryOverlapping(chromosome, position_of_snp, position_of_snp);
        
        // Right now assuming the above iterator provides me with reads.
        // Otherwise, I don't know.
        
        String bases = "";
        
        while(all_reads_in_region.hasNext()){
            
            SAMRecord read_in_region = all_reads_in_region.next();
            
            
            Character base_in_read = get_base_at_position(read_in_region, pos_int);
           //System.out.println("base_in_read: " + base_in_read);
            
            
            
            
            if( base_in_read == ref_allele.charAt(0) ){
                ref_overlap++;            
            } else if(base_in_read == alt_allele.charAt(0) ){
                alt_overlap++;
            }else if(base_in_read == '!'  || base_in_read == 'N'){
                continue;
            }else{
                no__overlap++;
            
            }
            
            
        }
        
        
        //Fuck this line below. cost me a day to figure out the error.
        all_reads_in_region.close();
        
        String string_for_output;
        string_for_output = Integer.toString(ref_overlap) + "\t" + 
                            Integer.toString(alt_overlap) + "\t" + 
                            Integer.toString(no__overlap);
        
        return(string_for_output);
    }
    
    
    public static Character get_base_at_position(SAMRecord sam_read,
                                             int pos_int
                                             ){
        
        int curr_pos = sam_read.getAlignmentStart();
        int end_pos = sam_read.getAlignmentEnd();
        
        
        
        String read = sam_read.getReadString();
        int idx_of_read = 0;
        int idx_of_cigar = 0;
        boolean debug = false;
        
        if(debug){
            System.out.println("positions");
            System.out.println("curr_pos: " + Integer.toString(curr_pos));
            System.out.println("end_pos: " +  Integer.toString( end_pos));
            System.out.println("SNP pos: " +  Integer.toString( pos_int));
        }
        
        while(curr_pos <= end_pos){
            
            CigarElement cigar = sam_read.getCigar().getCigarElement(idx_of_cigar);
            int final_pos;
            switch(cigar.getOperator()){
                case D:
                    curr_pos += cigar.getLength();
                    if(debug){
                        System.out.println("Found D in cigar: ");
                        System.out.println("current_pos + " + Integer.toString(cigar.getLength()));
                        System.out.println("final curr_pos: " + Integer.toString(curr_pos));

                    }
                    break;
                case EQ:
                    final_pos = curr_pos + cigar.getLength();
                    if(debug){
                        System.out.println("Found EQ in cigar: ");
                        System.out.println();
                        System.out.println("index_of_read + " + Integer.toString(cigar.getLength()));
                        System.out.println("final idx_of_read: " + Integer.toString(idx_of_read));
                        System.out.println();
                        System.out.println("current_pos + " + Integer.toString(cigar.getLength()));
                        System.out.println("final curr_pos: " + Integer.toString(curr_pos));
                        System.out.println("\nStarting iteration:");

                    }
                    while(curr_pos < final_pos){
                        if(debug){
                            System.out.println("index is now:" + Integer.toString(idx_of_read));
                            System.out.println("ref pos is now:" + Integer.toString(curr_pos));
                            System.out.println("The base at this pos is:" + read.charAt(idx_of_read));
                        }
                        
                        if(curr_pos  == pos_int){            
//                            System.out.println("index is now:" + Integer.toString(idx_of_read));
//                            System.out.println("ref pos is now:" + Integer.toString(curr_pos));
//                            System.out.println("The base at this pos is:" + read.charAt(idx_of_read));
                            return(read.charAt(idx_of_read));
                        }
                        
                        curr_pos++;
                        idx_of_read++;
                    }
                    break;
                case H:
                    System.out.println("found H line in Cigar");
                    return('!');
                case I:
                    idx_of_read += cigar.getLength();
                    if(debug){
                        System.out.println("Found I in cigar: ");
                        System.out.println("index_of_read + " + Integer.toString(cigar.getLength()));
                        System.out.println("final idx_of_read: " + Integer.toString(idx_of_read));

                    }
                    break;
                case M:    
                    final_pos = curr_pos + cigar.getLength();
                    if(debug){
                        System.out.println("Found M in cigar: ");
                        System.out.println();
                        System.out.println("index_of_read + " + Integer.toString(cigar.getLength()));
                        System.out.println("final idx_of_read: " + Integer.toString(idx_of_read));
                        System.out.println();
                        System.out.println("current_pos + " + Integer.toString(cigar.getLength()));
                        System.out.println("final curr_pos: " + Integer.toString(curr_pos));
                        System.out.println("\nStarting iteration:");

                    }
                    while(curr_pos < final_pos){
                        if(debug){
                            System.out.println("index is now:" + Integer.toString(idx_of_read));
                            System.out.println("ref pos is now:" + Integer.toString(curr_pos));
                            System.out.println("The base at this pos is:" + read.charAt(idx_of_read));
                        }
                        
                        if(curr_pos  == pos_int){            
//                            System.out.println("index is now:" + Integer.toString(idx_of_read));
//                            System.out.println("ref pos is now:" + Integer.toString(curr_pos));
//                            System.out.println("The base at this pos is:" + read.charAt(idx_of_read));
                            return(read.charAt(idx_of_read));
                        }
                        curr_pos++;
                        idx_of_read++;
                    }
                    break;
                case N:
                    curr_pos += cigar.getLength();
                    if(debug){
                        System.out.println("Found N in cigar: ");
                        System.out.println("current_pos + " + Integer.toString(cigar.getLength()));
                        System.out.println("final curr_pos: " + Integer.toString(curr_pos));

                    }
                    break;
                case P:
                    curr_pos += cigar.getLength();
                    if(debug){
                        System.out.println("Found P in cigar: ");
                        System.out.println("current_pos + " + Integer.toString(cigar.getLength()));
                        System.out.println("final curr_pos: " + Integer.toString(curr_pos));

                    }
                    break;
                case S:
                    idx_of_read += cigar.getLength();
                    if(debug){
                        System.out.println("Found S in cigar: ");
                        System.out.println("index_of_read + " + Integer.toString(cigar.getLength()));
                        System.out.println("final idx_of_read: " + Integer.toString(idx_of_read));

                    }
                    break;
                case X:
                    curr_pos += cigar.getLength();
                    idx_of_read += cigar.getLength();
                    if(debug){
                        System.out.println("Found X in cigar: ");
                        System.out.println();
                        System.out.println("index_of_read + " + Integer.toString(cigar.getLength()));
                        System.out.println("final idx_of_read: " + Integer.toString(idx_of_read));
                        System.out.println();
                        System.out.println("current_pos + " + Integer.toString(cigar.getLength()));
                        System.out.println("final curr_pos: " + Integer.toString(curr_pos));

                    }
                    break;
                    
            }
            idx_of_cigar++;
            
            if(curr_pos-1 > end_pos){
                System.out.println("Something goes wrong when bases based on cigar.");
                System.out.println(Integer.toString(curr_pos) +  ">" + Integer.toString(end_pos));
            }
            
        }
        //couldn't find anything, returning: '!'
        if(debug){
            //System.out.println("### coulnd't find overlapping base in this read. ###");
        }
        return('!');
    }
    
    

    
    /**
     *
     * @param individual_names
     * @return
     * @throws FileNotFoundException
     * @throws IOException
     */
    public static HashMap convert_individual_names(String[] individual_names, String coupling_location) throws FileNotFoundException, IOException{
        
        String coupling_loc = coupling_location;
        
        //This will be filled while reading the file
        ArrayList<String> sample_names_in_file     = new ArrayList<String>();
        ArrayList<String> individual_names_in_file = new ArrayList<String>();
        
        //This will be the one that is later returned
        HashMap ordered_sample_names = new HashMap();
        
        File coupling_file = new File(coupling_loc);
        
        FileInputStream fis;
        BufferedInputStream bis;
        DataInputStream dis;
        
        
        fis = new FileInputStream(coupling_file);

        // Here BufferedInputStream is added for fast reading.
        bis = new BufferedInputStream(fis);
        dis = new DataInputStream(bis);

        // dis.available() returns 0 if the file does not have more lines.
        int i = 0;
        while (dis.available() != 0) {
            
            String[] curr_line = dis.readLine().split("\t");
            individual_names_in_file.add(i, curr_line[0]); 
            sample_names_in_file.add(i, curr_line[1]);

        }

        // dispose all the resources after using them.
        fis.close();
        bis.close();
        dis.close();
        i=0;
        for(String i_name : individual_names){
            int index = individual_names_in_file.indexOf(i_name);
            if(index < 0){
                //System.out.println("When looking for \"" + i_name +"\" couln't find this in the coupling list."); // couldn't find this name in the following list:");
                //System.out.println(individual_names_in_file.toString());
                continue;
            }
            
            ordered_sample_names.put(sample_names_in_file.get(index), i);
            i++;
        }
        return(ordered_sample_names);
    
    }
    
// This was used when reading all bam_files at once, but the htsjdk did not allow for this.
//
//    private static ArrayList<SamReader> open_all_samfiles(File[] individual_files) throws Exception{
//        //don't know if i need this.
//        ArrayList<SamReader> all_bams = new ArrayList(); 
//        int i = 0;
//        
//        for(File i_file : individual_files){
//            
//            try{
//                all_bams.add();
//            } catch(Exception ex ){
//                System.err.println("Error in reading a file at: " + i_file.getAbsolutePath());
//                System.err.println("Iterator at: " + String.valueOf(i));                
//                System.err.println("Are you sure the bam files are correct?");
//                throw ex;
//            }
//            i++;
//        }
//    
//    
//        return(all_bams);
//    }
}




