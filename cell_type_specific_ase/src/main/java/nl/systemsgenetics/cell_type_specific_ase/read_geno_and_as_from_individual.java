/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package nl.systemsgenetics.cell_type_specific_ase;


import java.io.BufferedInputStream;
import java.io.DataInputStream;
import java.io.File;
import java.io.FileInputStream;
import java.io.FileNotFoundException;
import java.io.IOException;
import static java.lang.System.exit;
import java.lang.reflect.Array;
import java.util.Arrays;

import java.util.HashMap;
import java.util.HashSet;
import java.util.List;
import java.util.Set;
import org.molgenis.genotype.Allele;
import org.molgenis.genotype.Alleles;
import org.molgenis.genotype.RandomAccessGenotypeData;
import org.molgenis.genotype.Sample;
import org.molgenis.genotype.annotation.Annotation;
import org.molgenis.genotype.trityper.TriTyperGenotypeData;
import org.molgenis.genotype.util.FixedSizeIterable;
import org.molgenis.genotype.variant.GeneticVariant;
import org.molgenis.genotype.variant.GenotypeRecord;


import htsjdk.samtools.DefaultSAMRecordFactory;
import htsjdk.samtools.SAMFileWriter;
import htsjdk.samtools.SAMFileWriterFactory;
import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.SamInputResource;
import htsjdk.samtools.SamReader;
import htsjdk.samtools.SamReaderFactory;
import htsjdk.samtools.ValidationStringency;
import htsjdk.samtools.seekablestream.SeekableStream;




/**
 *
 * @author Adriaan van der Graaf
 */
public class read_geno_and_as_from_individual {
    
    public static void main(String[] args) throws IOException, Exception{
    
        String path_of_trityper = "/gcc/groups/lld/tmp01/projects/LL-imputed-20140306-TriTyper";
        
        RandomAccessGenotypeData dataset1;
        
        try {
            dataset1 = new TriTyperGenotypeData(new File(path_of_trityper));
        } catch (IOException ex) {
            System.err.println("Error reading dataset 1: " + path_of_trityper);
            throw ex;
        }
        
        
        HashMap<String, GeneticVariant> variantIdMap = dataset1.getVariantIdMap();
        Set<String> snpNames = variantIdMap.keySet();
        
        
        String[] individual_names;
        individual_names = dataset1.getSampleNames();
        
        String path = "/gcc/groups/lld/prm02/projects/rnaseq/alignment/";
        
        int i = 0;
        File[] individual_files = null;
        
        for(String i_individual : individual_names){
            
            String path_and_filename = path + i_individual + ".sorted.bam";
           
        
            individual_files[i] = new File(path_and_filename);

            
            i++;
        }
        
        SamReader[] all_bam_files;
        try{
            all_bam_files = open_all_samfiles(individual_files);
        }catch(Exception ex){
            System.err.println("I couldn't open all the bamfiles present in the trityper data");
            System.err.println("Please check that everything is correct, current individual names:");
            System.err.println(Arrays.toString(individual_names));
            throw(ex);
        }
        
        int number_of_individuals;         
        number_of_individuals = Arrays.asList(individual_names).size();
        System.out.println(Arrays.toString(individual_names));
        
        String[] sample_names = convert_individual_names(individual_names);
        
        System.out.println(Arrays.toString(sample_names));
        
//        i = 0;
//        for(String i_snp : snpNames){
//            
//            System.out.println(i_snp);
//            
//            GeneticVariant this_variant = variantIdMap.get(i_snp);
//            
//            String chromosome = this_variant.getSequenceName();
//            String position = String.valueOf(this_variant.getStartPos());
//            
//            // We only do analyses if we find a SNP and it is biallelic      
//            if(this_variant.isSnp() & this_variant.isBiallelic()){
//                
////                System.out.println(this_variant.getAllIds());
//                                
//                String row_of_table = get_allele_specific_overlap_at_snp(this_variant, 
//                                                                        individual_names,
//                                                                        chromosome,
//                                                                        position,
//                                                                        all_bam_files);
//                  
//                
//                i++;
//            }
//         
//        
//         if(i > 2){
//             break;
//         }
//        
//        }
    }
    
    
    public static String get_allele_specific_overlap_at_snp(GeneticVariant this_variant,
                                                            String[] individual_names, 
                                                            String chromosome,
                                                            String position, 
                                                            SamReader[] all_bam_files){
       
        Alleles all_variants = this_variant.getVariantAlleles();
        char ref_allele = all_variants.getAllelesAsChars()[0];
        char alt_allele = all_variants.getAllelesAsChars()[1];
        
        
        int num_of_individuals = Arrays.asList(individual_names).size();
        
        List variantsInPopulation = this_variant.getSampleVariants();
        if(num_of_individuals != variantsInPopulation.size() & variantsInPopulation.size() > 0){
            System.err.println("I could not find as many variants as there are individuals at SNP:");
            System.err.println(this_variant.getVariantId());
            System.err.println("This needs to be fixed.");
        }
        
       
        //now determine per individual the sample variants.
        //I'm assuming the ordering is correct
        
        String for_output = "";
        
        
        int i = 0;
        for(Object individual_allele : variantsInPopulation ){
            
            String allele = individual_allele.toString();
            
            String tempString = get_allele_specific_overlap_at_snp(all_bam_files[i],
                                                                   chromosome, 
                                                                   position, 
                                                                   allele);
            
//            for_output.concat();
        
        }
        
        
        return("hoihoi");
    }
    
    
    public static String[] convert_individual_names(String[] individual_names) throws FileNotFoundException, IOException{
        
        String coupling_loc = "/gcc/groups/lld/tmp01/projects/LL-imputed-20140306-TriTyper/gte_LL_final_genotyped.txt";
        
        //This will be filled while reading the file
        String[] sample_names_in_file = null;
        String[] individual_names_in_file = null;
        
        //This will be the one that is later returned
        String[] ordered_sample_names = null;
        
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
            individual_names_in_file[i] = curr_line[0];
            sample_names_in_file[i] = curr_line[1];

        }

        // dispose all the resources after using them.
        fis.close();
        bis.close();
        dis.close();
        i=0;
        for(String i_name : individual_names){
            int index = Arrays.asList(individual_names_in_file).indexOf(i_name);
            ordered_sample_names[i] = sample_names_in_file[index];
            i++;
        }
        return(ordered_sample_names);
    
    }
    
    
    
    

    public static String get_allele_specific_overlap_at_snp(SamReader bam_file,
                                                            String chromosome,
                                                            String position,
                                                            String allele){

        //TODO make sure all the stuff here works.
        
        return(null);
    }

    private static SamReader[] open_all_samfiles(File[] individual_files) throws Exception{
        //don't know if i need this.
        SamReader[] all_bams = null; 
        int i = 0;
        
        for(File i_file : individual_files){
            
            try{
                all_bams[i] = SamReaderFactory.makeDefault().open(i_file);
            } catch(Exception ex ){
                System.err.println("Error in reading a file at: " + i_file.getAbsolutePath());
                System.err.println("Are you sure the bam files are correct?");
                throw ex;
            }
            i++;
        }
    
    
        return(all_bams);
    }
}




