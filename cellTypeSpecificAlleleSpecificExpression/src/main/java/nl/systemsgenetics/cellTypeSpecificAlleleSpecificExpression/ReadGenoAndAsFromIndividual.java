/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */

package nl.systemsgenetics.cellTypeSpecificAlleleSpecificExpression;


import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.SAMRecordIterator;
import htsjdk.samtools.SamReader;
import htsjdk.samtools.SamReaderFactory;
import htsjdk.samtools.CigarElement;

import java.io.BufferedInputStream;
import java.io.DataInputStream;
import java.io.File;
import java.io.FileInputStream;
import java.io.FileNotFoundException;
import java.io.IOException;
import java.io.InputStream;
import java.io.PrintWriter;
import java.nio.file.Files;
import java.nio.file.Path;
import java.nio.file.Paths;

import java.util.Arrays;
import java.util.HashMap;
import java.util.Set;
import java.util.ArrayList;

import org.molgenis.genotype.Alleles;
import org.molgenis.genotype.RandomAccessGenotypeData;
import org.molgenis.genotype.trityper.TriTyperGenotypeData;
import org.molgenis.genotype.variant.GeneticVariant;
import org.molgenis.genotype.vcf.VcfGenotypeData;

import org.apache.commons.io.FilenameUtils;
import org.apache.commons.io.IOUtils;

import org.jdom.IllegalDataException;



/**
 *
 * @author Adriaan van der Graaf
 */
public class ReadGenoAndAsFromIndividual {
    
    public static void readGenoAndAsFromIndividual(String loc_of_bam1,
                                                   String genotype_loc, 
                                                   String coupling_location, 
                                                   String outputLocation, 
                                                   String snpLocation) throws IOException, Exception{
        
        if(GlobalVariables.verbosity >= 10){
            //Print ASREADS header
            System.out.println(    "---- Starting ASREADS for the following settings: ----");
            System.out.println(    "\t input bam:         " +  loc_of_bam1);
            System.out.println(    "\t genotype location: " +  genotype_loc);
            System.out.println(    "\t coupling file:     " +  coupling_location);
            System.out.println(    "\t output location:   " +  outputLocation);
            if(!snpLocation.equals("")){
                System.out.println("\t snp Location:      " +  snpLocation);
            }else{
                System.out.println("\t snp Location:      " +  "NONE");
            }

            System.out.println("------------------------------------------------------");
        }
        
        
        

        //parse command line arguments
        String LocOfBam;
        LocOfBam = loc_of_bam1;
        System.out.println("Location of bam file: ");
        System.out.println(LocOfBam);
        
        if (!new File(LocOfBam).exists()){
           throw new IllegalArgumentException("ERROR! Location of bam file is not an existing file. Exitting.");
        }else{
            if(GlobalVariables.verbosity >= 10){
                System.out.println("Location of bam file is an existing file, will continue.");
            }
        }
        
        
        
        RandomAccessGenotypeData TTdataSet;
        VcfGenotypeData VCFdataSet;
        HashMap<String, GeneticVariant> variantIdMap;
        String[] individual_names;
        
        String tabixLoc = genotype_loc + ".tbi";
        
        //open vcf dataset
        //based on extention and existance of both files. 
        if(FilenameUtils.isExtension(genotype_loc, "gz") && 
           new File(tabixLoc).exists() &&
           new File(genotype_loc).exists() 
          ){
            try {
                VCFdataSet = new VcfGenotypeData(new File(genotype_loc), new File(tabixLoc), 0.99);
                variantIdMap = VCFdataSet.getVariantIdMap();
                individual_names = VCFdataSet.getSampleNames();
            } catch (IOException ex) {
                System.err.println("Error reading vcf dataset: " + genotype_loc);
                throw new IllegalArgumentException();
            }
        
        }else if(new File(genotype_loc + "/GenotypeMatrix.dat").exists() ) {
            //assuming trityper dataset based on the genotype matrix
            try {
                TTdataSet = new TriTyperGenotypeData(new File(genotype_loc));
                variantIdMap = TTdataSet.getVariantIdMap();
                individual_names = TTdataSet.getSampleNames();
            } catch (IOException ex) {
                System.err.println("Error reading trityper dataset: " + genotype_loc);
                throw new IllegalArgumentException();
            }
        
        }else{
            throw new IllegalDataException("could not find a Trityper or vcf file in the genotype location");
        }
        
        

        
        
        //get the variants in the variantIdMAP
       
        Set<String> snpNames = variantIdMap.keySet();
        
        ArrayList<String> SNPsToAnalyze;
        SNPsToAnalyze = new ArrayList<String>();
        
        //If available, read the file with rs numbers.
        if(!snpLocation.equals("")){
            ArrayList<String>  includeSNPs = UtilityMethods.readFileIntoStringArrayList(snpLocation);
            int snpsNotFound = 0;
            for(String snp_to_include : includeSNPs){
                if(snpNames.contains(snp_to_include)){
                    SNPsToAnalyze.add(snp_to_include);
                }else{
                    snpsNotFound++;
                }
            }
            
            if(GlobalVariables.verbosity >= 1){
                System.out.println("WARNING: Did not find " + Integer.toString(snpsNotFound) + " out of " + Integer.toString(includeSNPs.size()) + " SNPs in the include file.");
            }
        }else{
            for(String snp_to_include : snpNames){
                SNPsToAnalyze.add(snp_to_include);
            }
        }
        int totalSnps = SNPsToAnalyze.size();
       
        //String path = "/gcc/groups/lld/tmp01/projects/bamFiles/";
        //sample_map contains all the individuals that are in the sample file.
        HashMap sample_map = convert_individual_names(individual_names, coupling_location);
        
        if(GlobalVariables.verbosity >= 10){
        System.out.println("Sample names were loaded.");
        }
        if(GlobalVariables.verbosity >= 100){
            System.out.println(sample_map.toString());
        }
       
        //Twice because my files have the .MERGED.sorted.bam suffix attached to them.
        String sample_name = FilenameUtils.getBaseName(FilenameUtils.getBaseName(FilenameUtils.getBaseName(LocOfBam)));
        
        if(GlobalVariables.verbosity >= 10){
            System.out.println("sample_name: " + sample_name);
            System.out.println("sample_map:  " + sample_map.toString());
        }
        
        Object sample_idx = sample_map.get(sample_name);   
       
        if(sample_idx == null){
            throw new IllegalArgumentException("Couldn't find the filename in the sample names. Quitting.");
        }

        int sample_index = Integer.parseInt(sample_idx.toString());
        
        if(GlobalVariables.verbosity >= 10){
            System.out.println("sample_index: " + sample_index);
        }
        
        //bam file path and filename
        File sample_file = new File(LocOfBam);

        PrintWriter writer = new PrintWriter(outputLocation, "UTF-8");
        SamReader bam_file = SamReaderFactory.makeDefault().open(sample_file);

        if(GlobalVariables.verbosity >= 10){
            System.out.println("Initialized for reading bam file");
        }
        
        

        InputStream IndexInputStream = new FileInputStream(LocOfBam+ ".bai");
        byte[] indexFileByteArray = IOUtils.toByteArray(IndexInputStream);
        FindBamTranscribedRegions checkRegion;
        checkRegion = new FindBamTranscribedRegions(indexFileByteArray);
 
        
        //iterate over all snps
        int i = 0;
        for(String i_snp : SNPsToAnalyze){
            
            
            
            GeneticVariant thisVariant = variantIdMap.get(i_snp);
            
            
            String chromosome = thisVariant.getSequenceName();
            String position = String.valueOf(thisVariant.getStartPos());
          

            
            // We only do analyses if we find a SNP and it is biallelic
            // However this is trityper data, so if we use
            // the allele count is used for the check of something. 
            
            
            if(thisVariant.isSnp() & thisVariant.isBiallelic() ){
                
                  
                if(!checkRegion.bamHasOverlap(Integer.parseInt(chromosome), thisVariant.getStartPos())){
                    
                    writer.println(chromosome + "\t" + position + "\t" + i_snp + "\t" + 
                               thisVariant.getVariantAlleles().getAllelesAsChars()[0] + "\t" +
                               thisVariant.getVariantAlleles().getAllelesAsChars()[1] + "\t" +
                               "0\t0\t0" + "\t" + 
                               Arrays.toString(thisVariant.getSampleVariants().get(sample_index).getAllelesAsChars()) //+ "\t" +
                               //Boolean.toString(this_variant.getSamplePhasing().get(sample_index))
                            );

                }else{
                
                
                    String row_of_table = get_allele_specific_overlap_at_snp(thisVariant, 
                                                                            sample_index,
                                                                            chromosome,
                                                                            position,
                                                                            bam_file);


                    //commented out the phasing part.

                    writer.println(chromosome + "\t" + position + "\t" + i_snp + "\t" + 
                                   thisVariant.getVariantAlleles().getAllelesAsChars()[0] + "\t" +
                                   thisVariant.getVariantAlleles().getAllelesAsChars()[1] + "\t" +
                                   row_of_table + "\t" + 
                                   Arrays.toString(thisVariant.getSampleVariants().get(sample_index).getAllelesAsChars()) //+ "\t" +
                                   //Boolean.toString(this_variant.getSamplePhasing().get(sample_index))
                                );
                }
            }
         
            i++;
            
            if((i % 10000 == 0) && (GlobalVariables.verbosity >= 10)){

                System.out.printf("Finished %d (%3.1f %%) SNPs\r", i, (double)i / (double)totalSnps * 100.0 );
            
            }
        }
        
        System.out.println("Finished ASreads for: " + loc_of_bam1);
        System.out.println("Output is located in: " + outputLocation);
        
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
        
        int position_of_snp = Integer.parseInt(position);                       
        
        //Check to make sure the variant position is not 0.
        if(position_of_snp <= 0){
            System.out.println("A SNP was read with a position lower than 1. This is illegal");
            System.out.println("Please adapt your genotype files by removing SNPs with these illegal positions");
            System.out.println("\tchr: " + chromosome + " pos: " + position);
            throw new IllegalDataException("Variant Position was less than 1");
        }
        
        SAMRecordIterator all_reads_in_region;
        try{
            all_reads_in_region = bam_file.queryOverlapping(chromosome, position_of_snp, position_of_snp);
        } catch(IllegalArgumentException e){
            System.out.println("Found an error when trying the following input:");
            System.out.println("chr:\t"+chromosome);
            System.out.println("pos:\t"+ position);
            System.out.println("If these values look correct, please make sure your bam file is sorted AND indexed by samtools.");
            System.out.println("If the problem persists, perhaps the chromosome (or sequence) are not the same in the genotype or bam file");
            all_reads_in_region = null;
            System.exit(0);
            
        }
        
        String bases = "";
        
        while(all_reads_in_region.hasNext()){
            
            SAMRecord read_in_region = all_reads_in_region.next();
            
            Character base_in_read = get_base_at_position(read_in_region, pos_int);
            if(GlobalVariables.verbosity >= 100){
                System.out.println("base_in_read: " + base_in_read);
            }
            
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
        
        //This line below cost me a day to figure out the error.
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
        
        if(GlobalVariables.verbosity >= 100){
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
                    if(GlobalVariables.verbosity >= 100){
                        System.out.println("Found D in cigar: ");
                        System.out.println("current_pos + " + Integer.toString(cigar.getLength()));
                        System.out.println("final curr_pos: " + Integer.toString(curr_pos));

                    }
                    break;
                case EQ:
                    final_pos = curr_pos + cigar.getLength();
                    if(GlobalVariables.verbosity >= 100){
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
                        if(GlobalVariables.verbosity >= 100){
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
                    if(GlobalVariables.verbosity >= 10){
                        System.out.println("found H line in Cigar, skipping read");
                    }
                    
                    return('!');
                case I:
                    idx_of_read += cigar.getLength();
                    if(GlobalVariables.verbosity >= 100){
                        System.out.println("Found I in cigar: ");
                        System.out.println("index_of_read + " + Integer.toString(cigar.getLength()));
                        System.out.println("final idx_of_read: " + Integer.toString(idx_of_read));

                    }
                    break;
                case M:    
                    final_pos = curr_pos + cigar.getLength();
                    if(GlobalVariables.verbosity >= 100){
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
                        if(GlobalVariables.verbosity >= 100){
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
                    if(GlobalVariables.verbosity >= 100){
                        System.out.println("Found N in cigar: ");
                        System.out.println("current_pos + " + Integer.toString(cigar.getLength()));
                        System.out.println("final curr_pos: " + Integer.toString(curr_pos));

                    }
                    break;
                case P:
                    curr_pos += cigar.getLength();
                    if(GlobalVariables.verbosity >= 100){
                        System.out.println("Found P in cigar: ");
                        System.out.println("current_pos + " + Integer.toString(cigar.getLength()));
                        System.out.println("final curr_pos: " + Integer.toString(curr_pos));

                    }
                    break;
                case S:
                    idx_of_read += cigar.getLength();
                    if(GlobalVariables.verbosity >= 100){
                        System.out.println("Found S in cigar: ");
                        System.out.println("index_of_read + " + Integer.toString(cigar.getLength()));
                        System.out.println("final idx_of_read: " + Integer.toString(idx_of_read));

                    }
                    break;
                case X:
                    curr_pos += cigar.getLength();
                    idx_of_read += cigar.getLength();
                    if(GlobalVariables.verbosity >= 100){
                        
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
        if(GlobalVariables.verbosity >= 100){
            //System.out.println("### coulnd't find overlapping base in this read. ###");
        }
        return('!');
    }
    
    

    
    /**
     *
     * @param individual_names
     * @param coupling_location
     * @return
     * @throws FileNotFoundException
     * @throws IOException
     * 
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

            //print the individual al line for checking.
            //System.out.println("line: " + " " + curr_line[0] +" " + curr_line[1] );
            
        }

        // dispose all the resources after using them.
        fis.close();
        bis.close();
        dis.close();
        
        
        
        i=0;
        for(String i_name : individual_names){
            int index = individual_names_in_file.indexOf(i_name);
            
            if(index >= 0){
                //We find a name in the genotype folder in the individual names stuff.
                ordered_sample_names.put(sample_names_in_file.get(index), i);
            }
            i++;
            
        }
        return(ordered_sample_names);
    
    }

}




