/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package nl.systemsgenetics.cellTypeSpecificAlleleSpecificExpression;

import com.google.common.collect.Iterables;
import java.io.File;
import java.io.IOException;
import java.io.PrintWriter;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.HashMap;
import java.util.HashSet;
import java.util.List;
import java.util.Set;
import org.apache.commons.io.FilenameUtils;
import org.jdom.IllegalDataException;

import org.molgenis.genotype.Alleles;
import org.molgenis.genotype.variant.GeneticVariant;
import org.molgenis.genotype.vcf.VcfGenotypeData;
import umcg.genetica.containers.Pair;

/**
 *
 * @author adriaan
 */
class PhasedEntry {
    
    
    public PhasedEntry(String asLocations, 
                       String couplingLoc, 
                       String outputLocation, 
                       String cellPropLoc,
                       String dispersionLocation,
                       String phasingLocation, 
                       String regionLocation) throws IOException, Exception {
       
        /**
         * This method will perform a (beta) binomial test for some test region. 
         * later additional features will be add.
         * 
         * currently the flow of the program:
         * 1. read all SNPs from AS files and add overdispersion and cellprop to the sampless
         * 2. read phasing and assign alleles for these snps
         * 3. load test regions and determine test snps.
         * 4. determine log likelihood for test-snps. (with some deduplication to speed up the process.)
         */

        
        // 1. Read all SNPs from AS files

        ArrayList<String> allFiles = UtilityMethods.readFileIntoStringArrayList(asLocations);

        ReadAsLinesIntoIndividualSNPdata asReader = new ReadAsLinesIntoIndividualSNPdata(asLocations);
        
        HashMap<String, ArrayList<IndividualSnpData> > snpHashMap = new HashMap<String, ArrayList<IndividualSnpData>>();
       
        HashMap<String, String> posNameMap = new HashMap<String, String>();
        
        //first determine overdispersion values per SNP.
        
        ArrayList<BetaBinomOverdispInSample>  dispersionParameters = new ArrayList<BetaBinomOverdispInSample>();
        
        
        String dispersionOutput = FilenameUtils.getFullPath(outputLocation) + 
                                  FilenameUtils.getBaseName(outputLocation) +
                                  "_dispersionFile.txt";
        
        
        
        if(dispersionLocation == null){
            PrintWriter dispersionWriter = new PrintWriter(dispersionOutput, "UTF-8");       
            dispersionWriter.write("Filename\tdispersion\n");
            
            int i=0;
            for(String asLoc : allFiles){
                dispersionParameters.add(new BetaBinomOverdispInSample(asLoc));
                dispersionWriter.printf("%s\t%.6f\n", 
                                        dispersionParameters.get(i).getSampleName(), 
                                        dispersionParameters.get(i).getOverdispersion()[0]
                                    );
                                    
                i++;
            }

            dispersionWriter.close();
            if(GlobalVariables.verbosity >= 10){
                System.out.println("--------------------------------------------------");
                System.out.println("Finished dispersion estimates for all individuals.");
                System.out.println("--------------------------------------------------");        
            }
        } else {
            ArrayList<String> DispersionLines = UtilityMethods.readFileIntoStringArrayList(dispersionLocation);
            
            //remove header line.
            DispersionLines.remove(0);
            
            for(String Line : DispersionLines){
                
                String sampleName = Line.split("\t")[0];
                double[] dispersion;
                dispersion = new double[] { Double.parseDouble(Line.split("\t")[1]) };

                //the file is later checked for correctness in when writing the values into individual SNP data.
                BetaBinomOverdispInSample fillerForDispersion = new BetaBinomOverdispInSample(sampleName, dispersion);
                dispersionParameters.add(fillerForDispersion);
                
            }
            if(GlobalVariables.verbosity >= 10){
                System.out.println("-------------------------");
                System.out.println("Using dispersion data from a previously entered file");
                System.out.println("-------------------------");
            }
        
        }
        
        
        //See if a phenotype file is there.
        boolean hasPhenoValue = false;
        ArrayList<String> phenoString = new ArrayList<String>();
        if(cellPropLoc != null){
            hasPhenoValue = true;
            phenoString = UtilityMethods.readFileIntoStringArrayList(cellPropLoc);
                    
        }
        
        //second reading of the ASfiles.
        while (true) {
            
            //read some stuff from the files.
            ArrayList<IndividualSnpData> tempSNPdata;
            tempSNPdata = asReader.getIndividualsFromNextLine();
            
            //No more individuals in the 
            if(tempSNPdata.isEmpty()) break;
            
            //I'm, assuming, hopefully safely that all snps are the same 
            //per line, based on checks done in the getIndividualsFromNextLine.

            String snpName   = tempSNPdata.get(0).getSnpName();
            String chr       = tempSNPdata.get(0).getChromosome();
            String posString = tempSNPdata.get(0).getPosition();
            
            //add dispersionValues to the SNPs:
            for(int j= 0; j < tempSNPdata.size(); j++){
            
                if(!tempSNPdata.get(j).getSampleName().equals(dispersionParameters.get(j).getSampleName())){
                    System.out.println(tempSNPdata.get(j).getSampleName());
                    System.out.println(dispersionParameters.get(j).getSampleName());
                    throw new IllegalDataException("the name of the individual in the dispersion data is not the same as the individual name in the SNP");
                }
                tempSNPdata.get(j).setDispersion(dispersionParameters.get(j).getOverdispersion()[0]);
                
                if(hasPhenoValue){
                    tempSNPdata.get(j).setCellTypeProp(Double.parseDouble(phenoString.get(j)));
                }
            
            }
            
            posNameMap.put(chr + ":" + posString, snpName);
            
            //take the SNP name and arraylist and put in the hashmap.
            snpHashMap.put(chr + ":" + posString, tempSNPdata);
            
        }
        
        if(GlobalVariables.verbosity >= 10){
            System.out.println("all AS info Snps were read");
        }
      
        
        
        // 2. Load test regions and determine the snps in the region.
       
        if(GlobalVariables.verbosity >= 10){
            System.out.println("Starting the assignment of snps to regions.");
        }
        
        ArrayList<GenomicRegion> allRegions;
        allRegions = ReadGenomicRegions(regionLocation);
        
       
        
        // 3. Read phasing info for these snps

        Pair<HashMap<String, ArrayList<IndividualSnpData>>, ArrayList<GenomicRegion>> phasedPair;
        phasedPair = addPhasingToSNPHashMap(snpHashMap, couplingLoc, allRegions, phasingLocation);
        
        snpHashMap = phasedPair.getLeft();
        allRegions = phasedPair.getRight();
        
        
        if(GlobalVariables.verbosity >= 10){
            System.out.println("Added phasing information to AS values of snps.");
            if(allRegions.size() == 0){
                System.out.println("No regions to found. exiting");
            }
        }
        
        /**
         * 4.  Start testing, per region.:
         * 
         * 4.1 Detemine the test snp in the region, this will be the reference value
         * 4.2 Determine the heterozygotes for the test snp.
         * 4.3 switch alt and ref values of the heterozygotes in the test region 
         *      respective of the test snp. add the new list to a binomial test.
         * 4.4 do the test in the BinomialTest.java and others in the future.
         * 
         * 
        */
        
         //write output to these files.
           
        PrintWriter writerBinom = new PrintWriter(
                                                  FilenameUtils.getFullPath(outputLocation) + 
                                                  FilenameUtils.getBaseName(outputLocation) + 
                                                  "_Binomial_results.txt"
                                           , "UTF-8");
        
        writerBinom.println(writeASEheader());
        
        PrintWriter writerBetaBinom = new PrintWriter(
                                                    FilenameUtils.getFullPath(outputLocation) + 
                                                    FilenameUtils.getBaseName(outputLocation) + 
                                                    "_BetaBinomial_results.txt"
                                                , "UTF-8");
        
        writerBetaBinom.println(writeASEheader());
        
        //no header add, because not implemented yet.
        PrintWriter writerCTSBinom = new PrintWriter(
                                                  FilenameUtils.getFullPath(outputLocation) + 
                                                  FilenameUtils.getBaseName(outputLocation) + 
                                                  "_CellTypeSpecificBinomial_results.txt"
                                           , "UTF-8");
        
        PrintWriter writerCTSBetaBinom = new PrintWriter(
                                                  FilenameUtils.getFullPath(outputLocation) + 
                                                  FilenameUtils.getBaseName(outputLocation) + 
                                                  "_CellTypeSpecificBetaBinomial_results.txt"
                                           , "UTF-8");
        
        
        for(GenomicRegion iRegion : allRegions){
            
            if(GlobalVariables.verbosity > 10){
                System.out.println(iRegion.getAnnotation());
            }
            // I may want to change this into all test SNPS needs to be implemented still.
            // compared to all snps in the region.

            ArrayList<String> snpsInRegion = iRegion.getSnpInRegions();
            
            ArrayList<IndividualSnpData> allHetsInRegion = new ArrayList<IndividualSnpData>();
            
            //Don't want to do this in every iteration in the next loop.
            
            
            for( String regionSnp : snpsInRegion ){
                allHetsInRegion.addAll(UtilityMethods.isolateOnlyHeterozygotesFromIndividualSnpData(snpHashMap.get(regionSnp)));
            }
            
            
            
            HashSet<String> combinationsDone = new HashSet<String>();
            
            HashMap<String, BinomialTest> storedBinomTests = new HashMap<String, BinomialTest>();
            HashMap<String, BetaBinomialTest> storedBetaBinomTests = new HashMap<String, BetaBinomialTest>();
            
            ///PLEASE NOTE, COVARIATE SPECIFIC FUNCTIONALITY HAS NOT YET BEEN IMPLEMENTED.
            //Plan is to use this in the future but keeping them in
            HashMap<String, CTSbinomialTest> storedCTSBinomTests = new HashMap<String, CTSbinomialTest>();
            HashMap<String, CTSBetaBinomialTest> storedCTSBetaBinomTests = new HashMap<String, CTSBetaBinomialTest>();
          
            
            
            for( String testSnp : snpsInRegion ){

                
                ArrayList<IndividualSnpData> hetTestSnps =  UtilityMethods.isolateOnlyHeterozygotesFromIndividualSnpData(snpHashMap.get(testSnp));
                
                //Check if the snpp has phasing, but also see if there are heterozygous SNPs in the region.
                try{
                    if(!hetTestSnps.get(0).hasPhasing()){
                        System.out.println("\tno phasing");
                        continue;
                    }
                } catch(Exception e ){
                    continue;
                }
                
                StringBuilder inputIdA = new StringBuilder(); 
                StringBuilder inputIdB = new StringBuilder(); 

                ArrayList<String> hetTestNames = new ArrayList<String>();
                

                for(IndividualSnpData hetSample : hetTestSnps){
                    
                    if(hetSample.hasPhasing()){
                        inputIdA.append(hetSample.sampleName);
                        inputIdA.append(hetSample.getPhasingFirst());

                        inputIdB.append(hetSample.sampleName);
                        inputIdB.append(hetSample.getPhasingSecond());

                        hetTestNames.add(hetSample.sampleName);
                    }
                    
                }
                
                String refStringA = inputIdA.toString();
                String refStringB = inputIdB.toString();
                
                if(hetTestSnps.size() >= GlobalVariables.minHets){
                
                    //make sure I don't have to do two tests double.
                    if(combinationsDone.contains(refStringA)){
                       
                       BinomialTest binomForAddition         = storedBinomTests.get(refStringA);
                       BetaBinomialTest betaBinomForAddition = storedBetaBinomTests.get(refStringA);
                       
                       //there is duplication here to make sure it is stored under the correct name.
                       if(binomForAddition == null){
                       
                           binomForAddition = storedBinomTests.get(refStringB);
                           binomForAddition.addAdditionalSNP(hetTestSnps.get(0).snpName, hetTestSnps.get(0).position);
                           
                           betaBinomForAddition = storedBetaBinomTests.get(refStringB);
                           betaBinomForAddition.addAdditionalSNP(hetTestSnps.get(0).snpName, hetTestSnps.get(0).position);
                           
                           storedBinomTests.put(refStringB, binomForAddition);
                           storedBetaBinomTests.put(refStringB, betaBinomForAddition);
                           
                       } else {                       
                           binomForAddition.addAdditionalSNP(hetTestSnps.get(0).snpName, hetTestSnps.get(0).position);
                           storedBinomTests.put(refStringA, binomForAddition);
                           
                           betaBinomForAddition.addAdditionalSNP(hetTestSnps.get(0).snpName, hetTestSnps.get(0).position);
                           storedBetaBinomTests.put(refStringA, betaBinomForAddition);
                       }
                       
                       continue;
                    }
                
                    
                    ArrayList<IndividualSnpData> phasedSNPsForTest = new ArrayList<IndividualSnpData>();
                    
                    Set<String> uniqueGeneSNPnames = new HashSet<String>();
                    
                    // The following loop determines which SNPs will be used for
                    // test data.
                    
                    for(int j = 0 ; j < allHetsInRegion.size() ; j++){
                        
                        
                        IndividualSnpData thisHet = allHetsInRegion.get(j);
                        
                        
                        int snpPos  = Integer.parseInt(thisHet.position);
                        
                        //First check if the Heterozygote is in the test region and if it has enough reads.
                        if((snpPos < iRegion.getStartPosition() || 
                            snpPos > iRegion.getEndPosition())  ||
                            !(thisHet.hasPhasing())             ||
                            !UtilityMethods.valid_heterozygote(thisHet)){
                            continue;
                        }
                       
                        String sampleName = thisHet.sampleName;
                        uniqueGeneSNPnames.add(thisHet.snpName);
                        
                        if(!hetTestNames.contains(thisHet.sampleName) || !thisHet.hasPhasing()){
                            continue;
                        }
                        
                        //this is the heterozygote to compare to.
                        IndividualSnpData hetToCompareTo = hetTestSnps.get(hetTestNames.indexOf(sampleName));
                        if(!hetToCompareTo.hasPhasing()) continue;
                        
                        if(hetToCompareTo.getPhasingFirst() != thisHet.getPhasingFirst()){
                            // because it is a heterozygote, we can assume that 
                            // first is 0 and second is 1, or the other way around.
                            // if the first in  this snp doesn't match the 
                            // first in test, we will have to switch the ref, alt
                            // alleles
                            int temp = thisHet.refNum;
                            thisHet.refNum = thisHet.altNum;
                            thisHet.altNum = temp;
                        }
                        
                        phasedSNPsForTest.add(thisHet);
                    
                    }
                    
                    //Make sure we have heterozygous SNPs in the gene region.ww
                    if(phasedSNPsForTest.isEmpty()){
                        continue;
                    }
                    
                    if(GlobalVariables.verbosity >= 10){
                        System.out.println("\n----------------------------------------");
                        System.out.println("Testing Region:                " + iRegion.getAnnotation());
                        System.out.println("With the following test SNP:   " + hetTestSnps.get(0).snpName);
                        System.out.println("Using the following gene SNPs: ");
                        int whatSNP = 0;
                        System.out.print("\t[ ");
                        for(String snpName : uniqueGeneSNPnames){
                            System.out.print(snpName);
                            if( (whatSNP % 4 == 3) && (whatSNP != uniqueGeneSNPnames.size() - 1) ) {
                                System.out.print(",\n\t  ");
                            }else if(whatSNP != uniqueGeneSNPnames.size() - 1){
                                System.out.print(", ");
                            }
                            whatSNP += 1; 
                        }
                        System.out.println(" ]");
                        System.out.println("----------------------------------------\n");
                    }
                    
                    
                    BinomialTest thisBinomTest;
                    thisBinomTest = BinomialTest.phasedBinomialTest(phasedSNPsForTest, iRegion,hetTestSnps.size() );
                    thisBinomTest.addAdditionalSNP(hetTestSnps.get(0).snpName, hetTestSnps.get(0).position);                    
                    thisBinomTest.setGenotype(hetTestSnps.get(0).genotype);
                    storedBinomTests.put(refStringA, thisBinomTest);
                    
                    
                    BetaBinomialTest thisBetaBinomTest;
                    thisBetaBinomTest = BetaBinomialTest.phasedBetaBinomialTest(phasedSNPsForTest, iRegion, hetTestSnps.size());
                    thisBetaBinomTest.addAdditionalSNP(hetTestSnps.get(0).snpName, hetTestSnps.get(0).position);
                    thisBetaBinomTest.setGenotype(hetTestSnps.get(0).genotype);
                    storedBetaBinomTests.put(refStringA, thisBetaBinomTest);
                    

                    //make sure we don't have to do the computationally intensive tests again.
                    combinationsDone.add(refStringA);
                    combinationsDone.add(refStringB);
                    
                }
            }
            
            for(String thisTestName : storedBinomTests.keySet()){
                
                BinomialTest thisBinomTest = storedBinomTests.get(thisTestName);
                BetaBinomialTest thisBetaBinomTest = storedBetaBinomTests.get(thisTestName);
                
                if(thisBinomTest.isTestPerformed()){
                    writerBinom.println(writeBinomialTestOutput(thisBinomTest));
                    writerBetaBinom.println(writeBetaBinomialTestOutput(thisBetaBinomTest));
                }
                
            }
            
            
        }
        
        //close the files
        writerBinom.close();
        writerBetaBinom.close();
        writerCTSBinom.close();
        writerCTSBetaBinom.close();
        
    }

    
    
    
    
    
    
    
    private Pair<HashMap<String, ArrayList<IndividualSnpData>>, ArrayList<GenomicRegion>> 
            addPhasingToSNPHashMap(HashMap<String, ArrayList<IndividualSnpData>> snpHashMap, 
                                   String couplingLoc, 
                                   ArrayList<GenomicRegion> genomicRegions,
                                   String vcfLoc
                               ) throws IOException {

        String tabixLoc = vcfLoc + ".tbi";

        VcfGenotypeData genotypeData = new VcfGenotypeData(new File(vcfLoc), new File(tabixLoc), GlobalVariables.variantProb);


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
                   found = true;
                   break;
                }
            }

            if(!found && GlobalVariables.verbosity >= 10){
                System.out.println("Couldn't find individual " + tempCouple[0] + " in sampleNames of VCF, continueing with the next.");
                System.out.println(Arrays.toString(sampleNames));
            }
        }
        
        if(GlobalVariables.verbosity >= 100){
            System.out.println("final coupling map:");
            System.out.println(couplingMap.toString());
        }

        // iterate over all regions, and 
        int snpsDone = 0;
        for(int regionIndicator = 0; regionIndicator < genomicRegions.size(); regionIndicator++){
            
            GenomicRegion iRegion = genomicRegions.get(regionIndicator);
            
            Iterable<GeneticVariant> VariantsInTestRegion;

            VariantsInTestRegion = genotypeData.getVariantsByRange(
                                                iRegion.getSequence(),
                                                iRegion.getTestStart() - 1 , //minus one because bigger than.
                                                iRegion.getTestEnd());
            
            Iterable<GeneticVariant> VariantsInGeneRegion;
            Iterable<GeneticVariant> VariantsInRegion;
            
            if(iRegion.HasTestRegion()){
                VariantsInGeneRegion = genotypeData.getVariantsByRange(
                                                    iRegion.getSequence(),
                                                    iRegion.getStartPosition() - 1 , //minus one because bigger than.
                                                    iRegion.getEndPosition());
                
                VariantsInRegion = Iterables.concat(VariantsInTestRegion, VariantsInGeneRegion);
            }else{
                VariantsInRegion = VariantsInTestRegion;
            }
            
            
           
            
            ArrayList<String> snpsInThisTestRegion = new ArrayList<String>();
            
            HashSet<String> snpsAdded = new HashSet<String>();
            
            for(GeneticVariant currentVariant : VariantsInRegion){
                
                if(GlobalVariables.verbosity >= 100){
                    System.out.println("Starting " + currentVariant.getPrimaryVariantId());
                }
                
                if(!(currentVariant.isSnp() && currentVariant.isBiallelic())){
                    //not a valid variant
                    if(GlobalVariables.verbosity >= 100){
                        System.out.println("continueing because not a valid variant (not biallelic or not a SNP)");
                    }
                    continue;
                }
                
                //12 June 2016, this part will be rewritten to be dependent 
                //on the position of a snp, not on the variant Id.
                
                //Commenting this out at the moment, probably not necessary.
//                if(currentVariant.getAllIds().isEmpty()){
//                    continue;
//                }
                
                
                //make sure there is no overlap between gene and testing variant.
                //Now based on (start) position of the variant.
                if(snpsAdded.contains(currentVariant.getSequenceName() + ":"  + Integer.toString(currentVariant.getStartPos()))){
                   continue;
                } else{
                    snpsAdded.add(currentVariant.getSequenceName() + ":"  + Integer.toString(currentVariant.getStartPos()));
                }
                
                
                String snpName = currentVariant.getPrimaryVariantId();
                String chr = currentVariant.getSequenceName();
                String posString =  Integer.toString(currentVariant.getStartPos());
                
                ArrayList<IndividualSnpData> oldSnpData = snpHashMap.get(chr + ":" + posString);
                

                
 
                
                if(oldSnpData == null){
                    if(GlobalVariables.verbosity >= 100){
                        System.out.println("Couldn't find SNP: " + snpName + " in AS reads, continueing with the next.");
                    }
                    continue;
                }
                
                if(!(oldSnpData.get(0).chromosome + ":" + oldSnpData.get(0).position).equals(chr + ":" + posString)){
                    if(GlobalVariables.verbosity >= 10){
                        System.out.println(oldSnpData.get(0).chromosome + ":" + oldSnpData.get(0).position);
                        
                        System.out.println("Couldn't find SNP: " + snpName + " in the correct position: " + chr + ":" + posString);
                        System.out.println(oldSnpData.get(0).getSnpName());
                    }
                    continue;
                }
                
                snpsInThisTestRegion.add(chr + ":" + posString);

                //If phasing is already add to this variant we may as well stop right here.
                if(oldSnpData.get(0).hasPhasing()){
                    //this assumes all individuals in oldSNPdata are phased at the same time.
                    //currently this is the case, but if it changes it, it might break.
                    continue;
                }
                
                
                ArrayList<IndividualSnpData> newSnpData = new ArrayList<IndividualSnpData>();
                List<Boolean> SamplePhasing = currentVariant.getSamplePhasing();
                List<Alleles> Variants = currentVariant.getSampleVariants();
                
                boolean passSNP = false;
                
                
                for(IndividualSnpData iAS : oldSnpData){
                    
                    int i = -1;
                    
                    try{
                        i = couplingMap.get(iAS.getSampleName());
                    }
                    catch(NullPointerException e){
                        System.out.println("Found a null pointer exception.");
                        System.out.println("couplingMap");
                        System.out.println(couplingMap.toString());
                        System.out.println("iAS sample Name");
                        System.out.println(iAS.getSampleName());
                        System.out.println("The values in the couplingMap should be exactly the same as the file names.");
                        System.out.println("Make sure the sample name is the same as in the coupling file.");
                        System.exit(2);
                    }
                    
                    //make sure there is phasing data available:
                    if(SamplePhasing.get(i)){

                        char Alt = iAS.getAlternative();

                        //Genotype IO used the phasing boolean, and the 
                        //order of the allele characters for the phasing.
                        char[] alleleChars;
                        alleleChars = Variants.get(i).getAllelesAsChars();


                        if(GlobalVariables.verbosity >= 100){
                            System.out.println(Arrays.toString(alleleChars));
                        }

                        //first assuming the allele is reference
                        //if not the case, then we change it.
                        int first = 0;
                        int second = 0;

                        if(alleleChars[0] == Alt){
                            first = 1;
                        }
                        if(alleleChars[1] == Alt){
                            second = 1;
                        }
                        
                        try{
                            iAS.setPhasing(first, second);
                        } catch(IllegalDataException e){
                            if(GlobalVariables.verbosity >= 10 && !passSNP){
                                System.out.println("Did not set phasing for variant" + snpName + " phasing does not match genotype.");
                            }
                            passSNP = true;
                        }
                    }
                    newSnpData.add(iAS);
                }
                
                //something went wrong in the SNP with phasing information and genotype from other data.
                if(passSNP){
                    newSnpData = oldSnpData;
                }
                
                
                //overwrite this into the snpHashMap
                
                
                snpHashMap.put(chr + ":" + posString, newSnpData);
                snpsDone++;
            }
            
            iRegion.setSnpInRegions(snpsInThisTestRegion);
            
            genomicRegions.set(regionIndicator, iRegion);
            
            if(GlobalVariables.verbosity >= 10 && regionIndicator % 1000 == 0){
                System.out.println("Finished adding phasing to " + Integer.toString(regionIndicator) + " regions, " + Integer.toString(snpsDone) + " snps.");
            }
        }
        
        //so I can output it easily
        Pair<HashMap<String, ArrayList<IndividualSnpData>>, ArrayList<GenomicRegion>> returnPair;
        
        returnPair = new Pair<HashMap<String, ArrayList<IndividualSnpData>>, ArrayList<GenomicRegion>>(snpHashMap, genomicRegions);
        
        return returnPair;
    }

    private ArrayList<GenomicRegion> ReadGenomicRegions(String regionsFile) throws IOException {
        
        ArrayList<String> stringArray = UtilityMethods.readFileIntoStringArrayList(regionsFile);        
        ArrayList<GenomicRegion> allRegions  = new ArrayList<GenomicRegion>();

        for(String iString : stringArray){
            String[] splitString = iString.split("\t");
            GenomicRegion tempRegion = new GenomicRegion();
            
            //the line is incorrect, will proceed with a warning.
            try{
                tempRegion.setAnnotation(splitString[0]); 
                tempRegion.setSequence(splitString[1]);
            } catch(ArrayIndexOutOfBoundsException e){

                System.out.println("The following line was not parsable as a valid region. Continueing with the regions parsed up untill now.");
                System.out.println(iString);
                return allRegions;
            }
            tempRegion.setStartPosition(Integer.parseInt(splitString[2]));
            tempRegion.setEndPosition(Integer.parseInt(splitString[3])); 
            
            //Try to see if there are test regions specified.
            try{
                //set the test region to what is found in the file
                tempRegion.setTestStart(Integer.parseInt(splitString[4]));
                tempRegion.setTestEnd(Integer.parseInt(splitString[5]));
                tempRegion.setHasTestRegion(true);
                
            } catch(IndexOutOfBoundsException e){
                //set the test region to be the genomic region.
                tempRegion.setTestStart(Integer.parseInt(splitString[2]));
                tempRegion.setTestEnd(Integer.parseInt(splitString[3]));
                tempRegion.setHasTestRegion(false);    
            }
            
            
            allRegions.add(tempRegion);
            
        }
        
        return allRegions;
    }
    
    public static String writeASEheader(){
        return "chr\t"
                + "<testStart>-"
                + "<testEnd>\t"
                + "<regionStart>-"
                + "<regionEnd>\t"
                + "regionName\t"
                + "pVal\t"
                + "chiSq\t"
                + "hetSNPs\t"
                + "refSum\t"
                + "altSum\t"
                + "binomRatio\t"
                + "genotype\t"
                + "allSNPpos\t"
                + "allSNPnames";
    
    }
    
    public static String writeBinomialTestOutput(BinomialTest thisTest){
        
        StringBuilder outputString = new StringBuilder(); 
    
        outputString.append(thisTest.getChromosome());
        outputString.append("\t");
        
        // The start and end of the test region
        outputString.append(Integer.toString(thisTest.testRegionStart));
        outputString.append("-");
        outputString.append(Integer.toString(thisTest.testRegionEnd));
        outputString.append("\t");
        
        // The start and end will be in the same field seperated by a dash
        outputString.append(Integer.toString(thisTest.startOfRegion));
        outputString.append("-");
        outputString.append(Integer.toString(thisTest.endOfRegion));
        outputString.append("\t");
        
        
        outputString.append(thisTest.RegionName);
        outputString.append("\t");
        
        outputString.append(String.format("% .9f",thisTest.getTestStatistics().getpVal()));
        outputString.append("\t");
        
        outputString.append(String.format("% 9.3f",thisTest.getTestStatistics().getChiSq()));
        outputString.append("\t");
        
        //number of hets in region
        outputString.append(Integer.toString(thisTest.totalTestSNPs));
        outputString.append("\t");
        
        
        int refSum = 0;
        int altSum = 0;
        for(int i : thisTest.getAsRef()) refSum += i;
        for(int i : thisTest.getAsAlt()) altSum += i;

        
        outputString.append(Integer.toString(refSum));
        outputString.append("\t");
        
        outputString.append(Integer.toString(altSum));
        outputString.append("\t");
        
        
        outputString.append(String.format("% 4.4f",thisTest.binomRatio));
        outputString.append("\t");

        outputString.append(thisTest.getGenotype());
        outputString.append("\t");
        
        
        ArrayList<String> allTestPositions = thisTest.getAdditionalPositions();
        
        for(String i : allTestPositions){
            outputString.append(i);
            outputString.append(";");
        }
        outputString.append("\t");
        
        ArrayList<String> allTestNames = thisTest.getAdditionalNames();
        
        for(String i : allTestNames){
            outputString.append(i);
            outputString.append(";");
        }
        outputString.append("\t");
        
        
        return outputString.toString();
    }
    
    public static String writeBetaBinomialTestOutput(BetaBinomialTest thisTest){
        
        StringBuilder outputString = new StringBuilder(); 
    
        outputString.append(thisTest.getChromosome());
        outputString.append("\t");
        
        // The start and end of the test region
        outputString.append(Integer.toString(thisTest.testRegionStart));
        outputString.append("-");
        outputString.append(Integer.toString(thisTest.testRegionEnd));
        outputString.append("\t");
        
        // The start and end will be in the same field seperated by a dash
        outputString.append(Integer.toString(thisTest.startOfRegion));
        outputString.append("-");
        outputString.append(Integer.toString(thisTest.endOfRegion));
        outputString.append("\t");
        
        outputString.append(thisTest.RegionName);
        outputString.append("\t");
        
        outputString.append(String.format("% .9f",thisTest.pVal));
        outputString.append("\t");
        
        outputString.append(String.format("% 9.3f", thisTest.chiSq));
        outputString.append("\t");
        
        //number of hets in region
        outputString.append(Integer.toString(thisTest.totalTestSNPs));
        outputString.append("\t");
        
        
        int refSum = 0;
        int altSum = 0;
        
        for(int i : thisTest.getAsRef()) refSum += i;
        for(int i : thisTest.getAsAlt()) altSum += i;

        
        outputString.append(Integer.toString(refSum));
        outputString.append("\t");
        
        outputString.append(Integer.toString(altSum));
        outputString.append("\t");
        
        outputString.append(String.format("% 4.4f",thisTest.binomRatio));
        outputString.append("\t");

        outputString.append(thisTest.getGenotype());
        outputString.append("\t");

        
        
        ArrayList<String> allTestPositions = thisTest.getAdditionalPositions();
        
        for(String i : allTestPositions){
            outputString.append(i);
            outputString.append(";");
        }
        outputString.append("\t");
        
        ArrayList<String> allTestNames = thisTest.getAdditionalNames();
        
        for(String i : allTestNames){
            outputString.append(i);
            outputString.append(";");
        }
        outputString.append("\t");
        
        
        return outputString.toString();
    }
    
}

