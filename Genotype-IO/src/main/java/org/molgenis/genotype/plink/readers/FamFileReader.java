/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package org.molgenis.genotype.plink.readers;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileInputStream;
import java.io.FileNotFoundException;
import java.io.IOException;
import java.io.InputStreamReader;
import java.nio.charset.Charset;
import java.util.ArrayList;
import java.util.LinkedHashMap;
import java.util.Map;
import java.util.regex.Pattern;
import static org.molgenis.genotype.GenotypeData.DOUBLE_PHENOTYPE_SAMPLE_ANNOTATION_NAME;
import static org.molgenis.genotype.GenotypeData.FATHER_SAMPLE_ANNOTATION_NAME;
import static org.molgenis.genotype.GenotypeData.MOTHER_SAMPLE_ANNOTATION_NAME;
import static org.molgenis.genotype.GenotypeData.SEX_SAMPLE_ANNOTATION_NAME;
import org.molgenis.genotype.Sample;
import org.molgenis.genotype.annotation.SexAnnotation;

/**
 *
 * @author patri
 */
public class FamFileReader {
    
    private static final Charset FILE_ENCODING = Charset.forName("UTF-8");
    private static final Pattern SEPARATOR_PATTERN = Pattern.compile("[ \\t]+");
    private static final org.apache.log4j.Logger LOGGER = org.apache.log4j.Logger.getLogger(FamFileReader.class);

    public static ArrayList<Sample> readFamFile(String famFilePath) throws FileNotFoundException, IOException {
        return readFamFile(new File(famFilePath));
    }
    
    public static ArrayList<Sample> readFamFile(File famFile) throws FileNotFoundException, IOException {

        ArrayList<Sample> samples = new ArrayList<Sample>();

        BufferedReader famFileReader = new BufferedReader(new InputStreamReader(new FileInputStream(famFile), FILE_ENCODING));

        String line;
        while ((line = famFileReader.readLine()) != null) {

            String[] elements = SEPARATOR_PATTERN.split(line);

            Map<String, Object> annotationValues = new LinkedHashMap<String, Object>();
            annotationValues.put(FATHER_SAMPLE_ANNOTATION_NAME, elements[2]);
            annotationValues.put(MOTHER_SAMPLE_ANNOTATION_NAME, elements[3]);
            annotationValues.put(SEX_SAMPLE_ANNOTATION_NAME, SexAnnotation.getSexAnnotationForPlink(Byte.parseByte(elements[4])));
            annotationValues.put(DOUBLE_PHENOTYPE_SAMPLE_ANNOTATION_NAME, Double.parseDouble(elements[5]));

            samples.add(new Sample(elements[1], elements[0], annotationValues));

        }

        LOGGER.info("Read " + samples.size() + " samples from " + famFile.getAbsolutePath());

        famFileReader.close();
        
        return samples;

    }

}