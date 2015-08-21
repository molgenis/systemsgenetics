/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package nl.systemsgenetics.cellTypeSpecificAlleleSpecificExpression;

import java.io.File;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.HashMap;
import java.util.List;
import java.util.Set;
import org.jdom.IllegalDataException;
import org.molgenis.genotype.Allele;
import org.molgenis.genotype.Alleles;
import org.molgenis.genotype.variant.GeneticVariant;
import org.molgenis.genotype.vcf.VcfGenotypeData;

/**
 *
 * @author adriaan
 */
class PhasedEntry {
    
    
    public PhasedEntry(String asLocations, String couplingLoc) throws IOException {
        /**
         * This method will perform a binomial test for some test region. 
         * later additional features will be add.
         * 
         * currently the flow of the program:
         * 1. read all SNPs from AS files
         * 2. read phasing and assign alleles for these snps
         * 3. load test regions and determine test snps.
         * 5. determine log likelihood for test-snps. (with some deduplication, hopefully)
         */

        // * 1. Read all SNPs from AS files

        ArrayList<String> allFiles = UtilityMethods.readFileIntoStringArrayList(asLocations);


        ReadAsLinesIntoIndividualSNPdata asReader = new ReadAsLinesIntoIndividualSNPdata(asLocations);
        
        HashMap<String, ArrayList<IndividualSnpData>> snpHashMap = new HashMap<String, ArrayList<IndividualSnpData>>();
        
        
        while (true) {
            //read some stuff from the files.
            ArrayList<IndividualSnpData> tempSNPdata;
            tempSNPdata = asReader.getIndividualsFromNextLine();
            if(tempSNPdata.isEmpty()) break;
            
            //I can safely assume all snps are the same per line, based on 
            //checks done in the getIndividualsFromNextLine.

            //take the SNP name and arraylist and put in the hashmap.
            snpHashMap.put(tempSNPdata.get(0).getSnpName() , tempSNPdata);
            
        }

        // 2. Read phasing info for these snps
        
        HashMap<String, ArrayList<IndividualSnpData>> phasedSnpHashMap;        
        phasedSnpHashMap = addPhasingToSNPHashMap(snpHashMap, couplingLoc);
        
        if(GlobalVariables.verbosity >= 10){
            System.out.println("Added phasing information to AS values of snps.");
        }
        
        // 3. Load test regions and determine test snps.
        
        ArrayList<GenomicRegion> allRegions;
        allRegions = ReadGenomicRegions("media/fast/GENETICA/implementInJava/CEUASREADS/genomicRegions.txt");
        
        
    }

    private HashMap<String, ArrayList<IndividualSnpData>> addPhasingToSNPHashMap(HashMap<String, ArrayList<IndividualSnpData>> snpHashMap, String couplingLoc) throws IOException {
            String vcfLoc;
            vcfLoc = "/media/fast/GENETICA/implementInJava/CEUGENOTYPES/special_vcf/" + 
                     "CEU.chr1.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.vcf.gz";
            
            String tabixLoc = vcfLoc + ".tbi";
            
            VcfGenotypeData genotypeData = new VcfGenotypeData(new File(vcfLoc), new File(tabixLoc), 0.99);
            

            //make a Hashmap with the coupling information correct.
            
            String[] sampleNames = genotypeData.getSampleNames();
            ArrayList<String> couplingList = UtilityMethods.readFileIntoStringArrayList(couplingLoc);
            
            HashMap<String, Integer> couplingMap = new HashMap<String, Integer>();
            
            for(String iSample : couplingList){
                
                String[] tempCouple = iSample.split("\t");
                boolean found = false;
                
                
                for(int i=0; i < sampleNames.length; i++){
                    
                    if(tempCouple[0].equals(sampleNames[i])){
                       couplingMap.put(tempCouple[1], i);
                       found =true;
                       break;
                    }
                }
                
                if(!found && GlobalVariables.verbosity >= 10){
                    System.out.println("couldn't find individual " + tempCouple[0] + " in sampleNames, continueing with the next.");
                }
            }
            if(GlobalVariables.verbosity >= 100){
                System.out.println("final coupling map:");
                System.out.println(couplingMap.toString());
            }
            
            HashMap<String, GeneticVariant> variantIdMap = genotypeData.getVariantIdMap();
            Set<String> phasedKeySet = variantIdMap.keySet();
            
            Set<String> genoKeySet = snpHashMap.keySet();
            
            int SNPsDone = 0; 
            
            for(String iKey : phasedKeySet ){
                
                GeneticVariant thisVariant = variantIdMap.get(iKey);
                
                List<Boolean> SamplePhasing;
                List<Alleles> Variants; 
                
                if(genoKeySet.contains(iKey)  && thisVariant.isSnp() && thisVariant.isBiallelic()){
                    
                    SamplePhasing = thisVariant.getSamplePhasing();
                    Variants = thisVariant.getSampleVariants();
                

                    //make sure we have phasing data.
                    if( SamplePhasing==null || Variants == null ){
                        if(GlobalVariables.verbosity >= 10){
                            System.out.println("Couldn't find snp: "+ iKey + " in genotyped data, continueing with the next.");
                        }
                        //no phasing so skip this snp
                        continue;
                    }

                    //add phasing infotmation to the arraylist of SNPdata:
                    ArrayList<IndividualSnpData> oldSnpData = snpHashMap.get(iKey);
                    ArrayList<IndividualSnpData> newSnpData = new ArrayList<IndividualSnpData>();
                    
                    
                    
                    for(IndividualSnpData iAS : oldSnpData){
                        
                        int i = couplingMap.get(iAS.getSampleName());

                        //make sure there is phasing data available:
                        if(SamplePhasing.get(i)){

                            char Alt = iAS.getAlternative();
                            
                            
                            char[] alleleChars;
                            alleleChars = Variants.get(i).getAllelesAsChars();


                            if(GlobalVariables.verbosity >= 100){
                                System.out.println(alleleChars);
                            }

                            //assuming the allele is reference
                            int first = 0;
                            int second = 0;

                            if(alleleChars[0] == Alt){
                                first = 1;
                            }
                            if(alleleChars[1] == Alt){
                                second = 1;
                            }

                            iAS.setPhasing(first, second);
                        }
                        newSnpData.add(iAS);
                    }
                    
                    //overwrite this into the snpHashMap
                    snpHashMap.put(iKey, newSnpData);
                    SNPsDone++;
                    if(SNPsDone % 10000 == 0 && GlobalVariables.verbosity >= 10){System.out.println("Finished " + SNPsDone + " Snps");}
                }
               
                
            }
            
        return snpHashMap;
    }

    private ArrayList<GenomicRegion> ReadGenomicRegions(String regionsFile) throws IOException {
        
        ArrayList<String> stringArray = UtilityMethods.readFileIntoStringArrayList(regionsFile);
        ArrayList<GenomicRegion> allRegions  = new ArrayList<GenomicRegion>();

        for(String iString : stringArray){
            String[] splitString = iString.split("\t");
            GenomicRegion tempRegion = new GenomicRegion();
            
            tempRegion.setAnnotation(splitString[0]); 
            tempRegion.setSequence(splitString[1]);
            tempRegion.setStartPosition(Integer.parseInt(splitString[2]));
            tempRegion.setStartPosition(Integer.parseInt(splitString[3])); 
            
            allRegions.add(tempRegion);
            
        }
        
        
        return allRegions;
    }
    
    

}

