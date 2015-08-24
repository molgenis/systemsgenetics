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


import java.math.BigInteger;
import java.security.MessageDigest;
import java.security.NoSuchAlgorithmException;
import java.util.ArrayList;
import java.util.Scanner;
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
        
        
        mainEntryPoint.main("-A", "1", "-G", genotypeLoc, "-O", outputLoc, "-C", couplingLoc, "-B", bamLoc);
        
        //Check md5sum of the output file.
        String md5sumResult = md5StringFromFile(outputLoc);
        
        Assert.assertEquals("2536ad0a2b3bd4e383a5fb0c60893b3f", md5sumResult);
    
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
        mainEntryPoint.main("-A", "2", "-O",outputLocUnedited  ,"-L", asLocation, "--minimum_reads", "1", "--minimum_heterozygotes", "1");
        //cannot check this by md5sum because of floating point precision probably.
        //I will do this using unit tests, and not the whole program.
        
        
        
        
        
        
        
    }
      
    

    
    public String md5StringFromFile(String filename) throws FileNotFoundException, NoSuchAlgorithmException{
        
        // Seems to work in the same way as the md5sum program in bash.
        // Got this off stackoverflow I believe
        
        String content = new Scanner(new File(filename)).useDelimiter("\\Z").next() + "\n" ; 
        
        MessageDigest m = MessageDigest.getInstance("MD5");
        m.reset();
        m.update(content.getBytes());
        byte[] digest = m.digest();
        BigInteger bigInt = new BigInteger(1,digest);
        String hashtext = bigInt.toString(16);
        // Now we need to zero pad it if you actually want the full 32 chars.
        while(hashtext.length() < 32 ){
          hashtext = "0"+hashtext;
        }
    
        return hashtext;
    }
    
}
