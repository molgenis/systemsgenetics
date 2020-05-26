/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package umcg.genetica.text.parsing;

import java.io.IOException;
import java.util.regex.Matcher;
import java.util.regex.Pattern;
import umcg.genetica.io.text.TextFile;

/**
 *
 * @author juha
 */
public class Demographics {

    private String annotationFile;

    public Demographics(String annotationFile) {
        this.annotationFile = annotationFile;
    }

    public static void main(String[] args) throws IOException {

//    new Demographics("/Data/GeneExpression/SampleAnnotation/GPL96GPL570/GPL96GPL570SampleAnnotation.txt").checkColumns("age");
//    new Demographics("/Data/GeneExpression/SampleAnnotation/GPL96GPL570/GPL96GPL570SampleAnnotation.txt").writeMatchingSamples("age", "/Data/Sasha/GPL96GPL570AgeSamples.txt");
//    new Demographics("/Data/GeneExpression/SampleAnnotation/GPL96GPL570/GPL96GPL570SampleAnnotation.txt").writeAgesMethylation("/Data/Sasha/GPL96GPL570AgeSamplesFinal.txt", 27, true);
//    new Demographics("/Data/GeneExpression/SampleAnnotation/GPL96GPL570/GPL96GPL570SampleAnnotation.txt").combineAgesAndRanges("/Data/Sasha/GPL96GPL570AgeSamplesFinal.txt", "/Data/Sasha/GPL96GPL570AgeSamplesFinalRanges.txt",
//        "/Data/Sasha/GPL96GPL570AgeSamplesWithRangesAveraged.txt");
//    new Demographics("/Data/GeneExpression/SampleAnnotation/GPL96GPL570/GPL96GPL570SampleAnnotation.txt").writeGenders("/Data/Sasha/GPL96GPL570FemaleSamples.txt", "/Data/Sasha/GPL96GPL570MaleSamples.txt", 27);
//    new Demographics("/Data/GeneExpressionFinal/SampleAnnotation/GPL8490/GPL8490_family_annotation_mesh2012.txt").writeGenders("/Data/MJ/GPL8490FemaleSamples.txt", "/Data/MJ/GPL8490MaleSamples.txt", 27);
//    new Demographics("/Data/GeneExpressionFinal/SampleAnnotation/GPL8490/GPL8490_family_annotation_mesh2012.txt").writeAges("/Data/MJ/GPL8490Ages.txt", 27);
//    new Demographics("D:\\UMCG\\Methylation_GPL8490\\GPL8490_raw_3112012\\GPL8490_family_annotation.txt").writeAgesMethylation("D:\\UMCG\\Methylation_GPL8490\\GPL8490_raw_3112012\\GPL8490Ages_New.txt", 27, true);
    new Demographics("D:\\UMCG\\Methylation_GPL8490\\GPL8490_raw_3112012\\GPL8490_family_annotation.txt").writeGenders("D:\\UMCG\\Methylation_GPL8490\\GPL8490_raw_3112012\\GPL8490FemaleSamples.txt", "D:\\UMCG\\Methylation_GPL8490\\GPL8490_raw_3112012\\GPL8490MaleSamples.txt", 27);
//    new Demographics("/Data/GeneExpression/SampleAnnotation/GPL96GPL570/GPL96GPL570SampleAnnotation.txt").writeMatchingSamples("female", "male", new int[]{7, 11, 30, 31, 33, 34}, "/Data/Sasha/GPL96GPL570FemaleSamples.txt");
//    new Demographics("/Data/GeneExpression/SampleAnnotation/GPL96GPL570/GPL96GPL570SampleAnnotation.txt").writeMatchingSamples("male", "female", new int[]{7, 11, 30, 31, 33, 34}, "/Data/Sasha/GPL96GPL570MaleSamples.txt");
    }

    private void writeGenders(String femaleFileName, String maleFileName, int gseCol) throws IOException {

        Pattern femaleP = Pattern.compile("\\bfemale\\b", Pattern.CASE_INSENSITIVE);
        Pattern femaleP2 = Pattern.compile("\\b(gender|sex)[ ]*[:=]?[ ]*f", Pattern.CASE_INSENSITIVE);
        Pattern maleP = Pattern.compile("\\bmale\\b", Pattern.CASE_INSENSITIVE);
        Pattern maleP2 = Pattern.compile("\\b(gender|sex)[ ]*[:=]?[ ]*m(?!atched)", Pattern.CASE_INSENSITIVE);

        TextFile fOut = new TextFile(femaleFileName, true);
        TextFile mOut = new TextFile(maleFileName, true);

        TextFile tf = new TextFile(annotationFile, false);
        String line = tf.readLine();
        while ((line = tf.readLine()) != null) {
            String[] split = line.split("\t");
            for (int i = 0; i < split.length; i++) {
                Matcher femaleM = femaleP.matcher(split[i]);
                Matcher femaleM2 = femaleP2.matcher(split[i]);
                Matcher maleM = maleP.matcher(split[i]);
                Matcher maleM2 = maleP2.matcher(split[i]);
                String gse = split[gseCol].trim().replace("\t", " ");
                if (femaleM.find(0) && !maleM.find(0)) {
                    fOut.writeln(split[12] + "\t" + gse + "\t" + split[0] + "\t" + split[3] + "\t" + split[13] + "\t" + split[i]);
                    break;
                } else if (maleM.find(0) && !femaleM.find(0)) {
                    mOut.writeln(split[12] + "\t" + gse + "\t" + split[0] + "\t" + split[3] + "\t" + split[13] + "\t" + split[i]);
                    break;
                } else if (femaleM2.find()) {
                    if (!split[i].toLowerCase().contains("gender: f/m")) {
                        fOut.writeln(split[12] + "\t" + gse + "\t" + split[0] + "\t" + split[3] + "\t" + split[13] + "\t" + split[i]);
                        break;
                    }
                } else if (maleM2.find()) {
                    if (!split[i].toLowerCase().contains("gender: m/f")) {
                        mOut.writeln(split[12] + "\t" + gse + "\t" + split[0] + "\t" + split[3] + "\t" + split[13] + "\t" + split[i]);
                        break;
                    }
                }
            }
        }
        tf.close();
        fOut.close();
        mOut.close();
    }

    private void writeAges(String fileName, int gseCol) throws IOException {

        Pattern ageP = Pattern.compile("\\bage[ ]*[:=]?[ ]*([0-9]+[\\.]?[0-9]*)", Pattern.CASE_INSENSITIVE);
        Pattern rangeP = Pattern.compile("\\bage[ ]*[:=]?[ ]*([0-9]+[ ]*(to|-)[ ]*[0-9]+[ ]*(?!week))", Pattern.CASE_INSENSITIVE);
        Pattern monthP = Pattern.compile("\\bage[ ]*[:=]?[ ]*([0-9]+[\\.]?[0-9]*[ ]*m(?!enopaus|iller))", Pattern.CASE_INSENSITIVE);
        Pattern weekP = Pattern.compile("\\bage[ ]*[:=]?[ ]*([0-9]+[\\.]?[0-9]*[ ]*(week|gestational week))", Pattern.CASE_INSENSITIVE);
        Pattern dayP = Pattern.compile("\\bage[ ]*[:=]?[ ]*([0-9]+[\\.]?[0-9]*[ ]*days(?! of symptoms))", Pattern.CASE_INSENSITIVE);

        TextFile out = new TextFile(fileName, true);
        TextFile rangeOut = new TextFile(fileName.replace(".txt", "") + "Ranges.txt", true);

        TextFile tf = new TextFile(annotationFile, false);
        String line = tf.readLine();
        
        while ((line = tf.readLine()) != null) {
            String[] split = line.split("\t");
            for (int i = 0; i < split.length; i++) {
                Matcher m1 = ageP.matcher(split[i]);
                String gse = split[gseCol].trim().replace("\t", " ");
                if (m1.find()) {
                    Matcher rangeM = rangeP.matcher(split[i]);
                    Matcher monthM = monthP.matcher(split[i]);
                    Matcher weekM = weekP.matcher(split[i]);
                    Matcher dayM = dayP.matcher(split[i]);
                    if (rangeM.find()) {
                        rangeOut.writeln(split[12] + "\t" + gse + "\t" + split[0] + "\t" + rangeM.group(1) + "\t" + split[3] + "\t" + split[13] + "\t" + split[i]);
                        break;
                    } else if (!monthM.find() && !weekM.find() && !dayM.find()) {
                        out.writeln(split[12] + "\t" + gse + "\t" + split[0] + "\t" + m1.group(1) + "\t" + split[3] + "\t" + split[13] + "\t" + split[i]);
                        break;
                    }
                }
            }
        }
        tf.close();
        out.close();
        rangeOut.close();
    }

    private void writeAgesMethylation(String fileName, int gseCol, boolean specialSelection) throws IOException {

        Pattern agePatern = Pattern.compile("\\bage[ ]*\\(*(y|yrs|years|months|[0-9]{4})*\\)*[:=]?[ ]*([0-9]+[\\.]?[0-9]*)", Pattern.CASE_INSENSITIVE);
        Pattern monthPatern1 = Pattern.compile("\\bage[ ]+\\(*(months|m)+\\)*[:=]?[ ]*([0-9]+[\\.]?[0-9]*)", Pattern.CASE_INSENSITIVE);
        Pattern monthPatern2 = Pattern.compile("\\bage[ ]*[:=]?[ ]*([0-9]+[\\.]?[0-9]*)[ ]*\\(*(month|m(?!enopaus|iller))+\\)*", Pattern.CASE_INSENSITIVE);
        Pattern rangeP = Pattern.compile("\\bage[ ]*[:=]?[ ]*([0-9]+[ ]*(to|-)[ ]*[0-9]+[ ]*(?!week))", Pattern.CASE_INSENSITIVE);

        Pattern weekPatern = Pattern.compile("\\bage[ ]*[:=]?[ ]*([0-9]+[\\.]?[0-9]*)[ ]*\\(*week[s]*\\)*", Pattern.CASE_INSENSITIVE);
        Pattern dayPatern = Pattern.compile("\\bage[ ]*[:=]?[ ]*([0-9]+[\\.]?[0-9]*)[ ]*\\(*day[s]*\\)*", Pattern.CASE_INSENSITIVE);

        Pattern ageSpecial1Patern = Pattern.compile("\\bageatdraw[:=]?[ ]*([0-9]+[\\.]?[0-9]*)", Pattern.CASE_INSENSITIVE);
        Pattern ageSpecial2Patern = Pattern.compile("\\bageatdiagnosis[:=]?[ ]*([0-9]+[\\.]?[0-9]*)", Pattern.CASE_INSENSITIVE);
        Pattern ageSpecial3Patern = Pattern.compile("\\bdurationt1d[:=]?[ ]*([0-9]+[\\.]?[0-9]*)", Pattern.CASE_INSENSITIVE);
        Pattern ageSpecial4Patern = Pattern.compile(".*age at collection \\(months\\): ([0-9]+[\\.]?[0-9]*).*", Pattern.CASE_INSENSITIVE);
        Pattern ageSpecial5Patern = Pattern.compile("\\bageatrecruitment[:=]?[ ]*([0-9]+[\\.]?[0-9]*)", Pattern.CASE_INSENSITIVE);
        Pattern maternalAgePatern = Pattern.compile("\\bmaternal age", Pattern.CASE_INSENSITIVE);


        TextFile out = new TextFile(fileName, true);
        TextFile rangeOut = new TextFile(fileName.replace(".txt", "Ranges.txt"), true);

        int numberMatches = 0;

        TextFile tf = new TextFile(annotationFile, false);
        String line = tf.readLine();
        
        while ((line = tf.readLine()) != null) {
            String[] split = line.split("\t");
            for (int i = 0; i < split.length; i++) {
                Matcher m1 = agePatern.matcher(split[i]);
                String gse = split[gseCol].trim().replace("\t", " ");
                if (m1.find()) {
                    Matcher monthM = monthPatern1.matcher(split[i]);
                    Matcher monthM2 = monthPatern2.matcher(split[i]);
                    Matcher maternalM = maternalAgePatern.matcher(split[i]);
                    Matcher weekM = weekPatern.matcher(split[i]);
                    Matcher dayM = dayPatern.matcher(split[i]);
                    Matcher rangeM = rangeP.matcher(split[i]);

                    if (rangeM.find(0)) {
                        rangeOut.writeln(split[12] + "\t" + gse + "\t" + split[0] + "\t" + rangeM.group(1) + "\t" + split[3] + "\t" + split[13] + "\t" + split[i]);
                        continue;
                    }

                    if (dayM.find(0) || weekM.find(0) || maternalM.find(0)) {
                        continue;
                    }

                    if (!monthM.find(0) && !monthM2.find(0) && Double.parseDouble(m1.group(2)) > 125) {
                        continue;
                    }

                    if (!monthM.find(0) && !monthM2.find(0)) {
                        out.writeln(split[12] + "\t" + gse + "\t" + split[0] + "\t" + m1.group(2) + "\t" + split[3] + "\t" + split[13] + "\t" + split[i]);
                        numberMatches++;
                        break;
                    } else if (monthM.find(0) || monthM2.find(0)) {
                        out.writeln(split[12] + "\t" + gse + "\t" + split[0] + "\t" + Double.parseDouble(m1.group(2)) / 12 + "\t" + split[3] + "\t" + split[13] + "\t" + split[i]);
                        numberMatches++;
                        break;
                    }
                }

                if (specialSelection) {
                    Matcher special = ageSpecial1Patern.matcher(split[i]);
                    Matcher special0 = ageSpecial4Patern.matcher(split[i]);
                    Matcher special1 = ageSpecial2Patern.matcher(split[i]);
                    Matcher special2 = ageSpecial3Patern.matcher(split[i]);
                    Matcher special3 = ageSpecial5Patern.matcher(split[i]);

                    if (special.find(0)) {
                        out.writeln(split[12] + "\t" + gse + "\t" + split[0] + "\t" + special.group(1) + "\t" + split[3] + "\t" + split[13] + "\t" + split[i]);
                        numberMatches++;
                        break;
                    } else if(special3.find(0)){
                        out.writeln(split[12] + "\t" + gse + "\t" + split[0] + "\t" + special3.group(1) + "\t" + split[3] + "\t" + split[13] + "\t" + split[i]);
                        numberMatches++;
                        break;
                    } else if(special0.find(0)){
                        out.writeln(split[12] + "\t" + gse + "\t" + split[0] + "\t" + Double.parseDouble(special0.group(1)) / 12 + "\t" + split[3] + "\t" + split[13] + "\t" + split[i]);
                        numberMatches++;
                        break;
                    } else if (special1.find(0) && special2.find(0)) {
                        out.writeln(split[12] + "\t" + gse + "\t" + split[0] + "\t" + (Double.parseDouble(special1.group(1)) + Double.parseDouble(special2.group(1))) + "\t" + split[3] + "\t" + split[13] + "\t" + split[i]);
                        numberMatches++;
                        break;
                    }
                }

            }
        }

        System.out.println(numberMatches);
        tf.close();
        out.close();
        rangeOut.close();
    }

    private void writeMatchingSamples(String matchThisWord, String fileName) throws IOException {

        Pattern p1 = Pattern.compile("\\b" + matchThisWord + "\\b", Pattern.CASE_INSENSITIVE);

        TextFile out = new TextFile(fileName, true);

        TextFile tf = new TextFile(annotationFile, false);
        String line = tf.readLine();
        while ((line = tf.readLine()) != null) {
            String[] split = line.split("\t");
            for (int i = 0; i < split.length; i++) {
                Matcher m1 = p1.matcher(split[i]);
                if (m1.find()) {
                    out.writeln(split[0] + "\t" + split[i]);
                }
            }
        }
        tf.close();
        out.close();
    }

    private void writeMatchingSamples(String matchThisWord, String dontMatchThisWord, int[] colsToSearch, String fileName) throws IOException {

        Pattern p1 = Pattern.compile("\\b" + matchThisWord + "\\b", Pattern.CASE_INSENSITIVE);
        Pattern p2 = Pattern.compile("\\b" + dontMatchThisWord + "\\b", Pattern.CASE_INSENSITIVE);

        TextFile out = new TextFile(fileName, true);

        TextFile tf = new TextFile(annotationFile, false);
        String line = tf.readLine();
        while ((line = tf.readLine()) != null) {
            String[] split = line.split("\t");
            boolean matches = false;
            for (int col : colsToSearch) {
                if (split.length > col) {
                    Matcher m1 = p1.matcher(split[col]);
                    Matcher m2 = p2.matcher(split[col]);
                    if (m1.find() && !m2.find()) {
                        matches = true;
                        break;
                    }
                }
            }
            if (matches) {
                out.writeln(split[0]);
            }
        }
        tf.close();
        out.close();
    }

    private void checkColumns(String query) throws IOException {
        TextFile tf = new TextFile(annotationFile, false);
        String line = tf.readLine();
        int[] counts = new int[line.split("\t").length];
        while ((line = tf.readLine()) != null) {
            String[] split = line.split("\t");
            for (int col = 0; col < split.length; col++) {
                if (col < counts.length && split[col].contains(query)) {
                    counts[col]++;
                }
            }
        }
        for (int col = 0; col < counts.length; col++) {
            if (counts[col] > 0) {
                System.out.println(col + "\t" + counts[col]);
            }
        }
    }

    private void combineAgesAndRanges(String agefile, String rangefile, String outfile) throws IOException {

        TextFile in = new TextFile(agefile, false);
        TextFile out = new TextFile(outfile, true);
        String line = null;
        while ((line = in.readLine()) != null) {
            out.writeln(line);
        }
        in.close();
        in = new TextFile(rangefile, false);
        while ((line = in.readLine()) != null) {
            String[] split = line.split("\t");
            String range = split[3].trim();
            int indexOf = range.indexOf("-");
            double age = -1;
            if (indexOf > 0) {
                int age1 = Integer.parseInt(range.substring(0, indexOf));
                int age2 = Integer.parseInt(range.substring(indexOf + 1));
                if (age2 - age1 < 19) {
                    age = (age2 + age1) / 2d;
                }
            } else {
                indexOf = range.indexOf(" to ");
                if (indexOf > 0) {
                    int age1 = Integer.parseInt(range.substring(0, indexOf));
                    int age2 = Integer.parseInt(range.substring(indexOf + 4));
                    if (age2 - age1 < 19) {
                        age = (age2 + age1) / 2d;
                    }
                }
            }
            if (age > 0) {
                String delim = "";
                for (int i = 0; i < 3; i++) {
                    out.write(delim + split[i]);
                    delim = "\t";
                }
                out.write(delim + age);
                for (int i = 4; i < split.length; i++) {
                    out.write(delim + split[i]);
                }
                out.writeln();
            }
        }
        in.close();
        out.close();
    }
}
