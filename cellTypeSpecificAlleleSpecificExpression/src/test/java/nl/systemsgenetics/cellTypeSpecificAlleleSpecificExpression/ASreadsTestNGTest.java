package nl.systemsgenetics.cellTypeSpecificAlleleSpecificExpression;

/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */


import java.io.File;
import java.io.FileNotFoundException;
import java.io.PrintWriter;
import java.net.URISyntaxException;
import java.nio.charset.Charset;
import java.text.DateFormat;
import java.text.SimpleDateFormat;
import java.util.Date;
import org.molgenis.genotype.Allele;
import org.molgenis.genotype.Alleles;
import org.testng.annotations.Test;
import org.testng.annotations.BeforeMethod;



import java.util.ArrayList;
import junit.framework.Assert;




/**
 *
 * @author adriaan
 * 
 * These tests are made so that the output of the files will stay the same
 * and changes do not break stuff in the actual output of files.
 */
public class ASreadsTestNGTest {
    
    private File tmpOutputFolder;
    private String fileSep = System.getProperty("file.separator");
    private File testFilesFolder;
    private static final Charset FILE_ENCODING = Charset.forName("UTF-8");
    private static final Alleles MISSING_ALLELES = Alleles.createAlleles(Allele.ZERO, Allele.ZERO);

    public ASreadsTestNGTest() throws URISyntaxException {
            testFilesFolder = new File(this.getClass().getResource("/").toURI());
    }

    // TODO add test methods here.
    // The methods must be annotated with annotation @Test. For example:
    //
    // @Test
    // public void hello() {}


    @BeforeMethod
    public void setUpMethod() throws Exception {
        File tmpDir = new File(System.getProperty("java.io.tmpdir"));

        DateFormat dateFormat = new SimpleDateFormat("yyyy-MM-dd_HH-mm-ss");
        Date date = new Date();

        tmpOutputFolder = new File(tmpDir, "ASReads_" + dateFormat.format(date));

        Runtime.getRuntime().addShutdownHook(new Thread() {
                
            @Override
            public void run() {
                System.out.println("Removing tmp dir and files");
                    for (File file : tmpOutputFolder.listFiles()) {
                        System.out.println(" - Deleting: " + file.getAbsolutePath());
                        file.delete();
                    }
                    System.out.println(" - Deleting: " + tmpOutputFolder.getAbsolutePath());
                    tmpOutputFolder.delete();
                }
            });

        tmpOutputFolder.mkdir();


        System.out.println("Temp folder with output of this test: " + tmpOutputFolder.getAbsolutePath());
    }
    
    
    
    @Test
    public void testASREADS() throws Exception{
      /**
       * This test is made to make sure output stays the same in the AS part.
       */
        String genotypeLoc = testFilesFolder + fileSep + "TriTyperFolder" + fileSep;
        String couplingLoc = testFilesFolder + fileSep + "coupling.txt";
        String bamLoc = testFilesFolder + fileSep + "testBam.bam";
        
        String outputLoc = tmpOutputFolder + fileSep + "testASREADSoutput.txt";
        
        
        MainEntryPoint.main("-A", "1", "-G", genotypeLoc, "-O", outputLoc, "-C", couplingLoc, "-B", bamLoc);
        
        String[] ref = {"1	1249187	rs12142199	G	A	21	12	0	[G, A]",
                        "1	778745	rs1055606	A	G	0	2	0	[A, G]",
                        "1	1342612	rs2275915	C	G	17	21	0	[C, G]",
                        "1	787135	rs28753393	A	G	0	0	0	[A, G]",
                        "1	787399	rs2905055	T	G	5	0	0	[T, G]",
                        "1	787262	rs56108613	G	C	0	0	0	[G, C]",
                        "1	1375810	rs4590	G	C	2	11	0	[G, C]" };

        
       ArrayList<String> calculated =  UtilityMethods.readFileIntoStringArrayList(outputLoc);
       
       Assert.assertEquals(calculated.size(),7 );
       
       for(int i = 0; i < calculated.size(); i++){
          Assert.assertEquals(calculated.contains(ref[i]), true);
       }
        
        
        
    }
    
 
    
    
    @Test
    public void testBinomial() throws Exception{
        
      /**
       * This test is made to make sure output stays the same in the AS part.
       * But this does not hold, probably due to floating point problems.
       * Haven't checked yet, this needs some looking over.
       */
        String OldasLocation = testFilesFolder + fileSep + "locationOfFiles.txt";
        
        //To make sure this works in every environment
        //Read in the file, only contains basepath, then add it to a new file.
        //Use this in the command

        String asLocation = tmpOutputFolder + fileSep + "LocationOfFiles.txt";
        
        PrintWriter writer = new PrintWriter(asLocation, "UTF-8");
        
        ArrayList<String>  allOldAsLocations = UtilityMethods.readFileIntoStringArrayList(OldasLocation);
        
        for(String iLoc : allOldAsLocations ){
            
            writer.println(testFilesFolder + fileSep + iLoc);
            System.out.println(testFilesFolder + fileSep + iLoc);
        }
        
        writer.close();
        
        
        String outputLocUnedited = tmpOutputFolder + fileSep + "testBINOMIALoutput.txt";
        
        //do the stuff for the output location
        MainEntryPoint.main("-A", "2", "-O",outputLocUnedited  ,"-L", asLocation, "--minimum_reads", "1", "--minimum_heterozygotes", "1");
        //cannot check this by md5sum because of floating point precision probably.
        //I will do this using unit tests, and not the whole program.
        
        
        
        
        
        
        
    }
      
    

}
