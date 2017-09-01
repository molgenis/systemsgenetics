package umcg.genetica.io.fasta;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileReader;
import java.io.IOException;
import java.util.Collections;
import java.util.HashSet;
import java.util.LinkedHashMap;
import java.util.Set;
import java.util.regex.Matcher;
import java.util.regex.Pattern;

/**
 *
 * @author Patrick Deelen
 */
public class ReferenceGenomeFasta {

    private final LinkedHashMap<String, LargeByteArray> chromosomes;
    private static final Pattern FASTA_HEADER_PATTERN = Pattern.compile("^\\>(\\S+).*:(\\d+):\\d+$");
    private static final Pattern FASTA_HEADER_PATTERN_ENSEMBL = Pattern.compile("^\\>(\\S+).*:(\\d+):\\d+ REF$");
    private static final Pattern FASTA_HEADER_CHR_PATTERN = Pattern.compile("^\\>(\\S+)");

    /**
     * CHR 1-22, X, Y and MT
     */
    public static final Set<String> HUMAN_NORMAL_CHR;

    static {
        HashSet<String> humanChr = new HashSet<String>(24);
        humanChr.add("1");
        humanChr.add("2");
        humanChr.add("3");
        humanChr.add("4");
        humanChr.add("5");
        humanChr.add("6");
        humanChr.add("7");
        humanChr.add("8");
        humanChr.add("9");
        humanChr.add("10");
        humanChr.add("11");
        humanChr.add("12");
        humanChr.add("13");
        humanChr.add("14");
        humanChr.add("15");
        humanChr.add("16");
        humanChr.add("17");
        humanChr.add("18");
        humanChr.add("19");
        humanChr.add("20");
        humanChr.add("21");
        humanChr.add("22");
        humanChr.add("X");
        humanChr.add("Y");
        humanChr.add("MT");

        HUMAN_NORMAL_CHR = Collections.unmodifiableSet(humanChr);
    }

    public ReferenceGenomeFasta(File fastaFile, Set<String> chrIncludeFilter) throws IOException, Exception {

        chromosomes = new LinkedHashMap<String, LargeByteArray>(32);

        BufferedReader fasteReader = new BufferedReader(new FileReader(fastaFile));

        String line;
        LargeByteArray currentSeq = null;
        long pos = 0;
        boolean skip = false;

        while ((line = fasteReader.readLine()) != null) {

            if (line.length() == 0) {
                continue;
            }

            if (line.charAt(0) == '>') {

//				System.out.println("Previous chr bases read: " + pos);
                Matcher headerChrMatcher = FASTA_HEADER_CHR_PATTERN.matcher(line);

                if (!headerChrMatcher.find()) {
                    throw new Exception("Error parsing reference genome fasta header: " + line);
                }

                String chr = headerChrMatcher.group(1);

                if (!chrIncludeFilter.contains(chr)) {
                    skip = true;
                    continue;
                }

                skip = false;

                Matcher headerMatcher = FASTA_HEADER_PATTERN.matcher(line);

                if (!headerMatcher.matches()) {
                    throw new Exception("Error parsing reference genome fasta header: " + line);
                }
//				System.out.println("----");
//				System.out.println(line);

                long lenght = Long.parseLong(headerMatcher.group(2));

                if (chr.equals("Y") && lenght == 59034049) {
                    lenght = 59373566;
                }

//				System.out.println(chr + " " + lenght);
                pos = 0;
                currentSeq = new LargeByteArray(lenght);
                chromosomes.put(chr, currentSeq);

//				System.out.println(currentSeq.getSize());
            } else if (!skip) {

                for (int n = 0; n < line.length(); ++n) {
                    currentSeq.setQuick(pos++, (byte) line.charAt(n));
                }

            }

        }

    }

    public ReferenceGenomeFasta(File fastaFile) throws IOException, Exception {

        chromosomes = new LinkedHashMap<String, LargeByteArray>(32);

        BufferedReader fasteReader = new BufferedReader(new FileReader(fastaFile));

        String line;
        LargeByteArray currentSeq = null;
        long pos = 0;

        while ((line = fasteReader.readLine()) != null) {

            if (line.length() == 0) {
                continue;
            }

            if (line.charAt(0) == '>') {

//				System.out.println("Previous chr bases read: " + pos);
//                                System.out.println(line);
                Matcher headerChrMatcher = FASTA_HEADER_CHR_PATTERN.matcher(line);

                if (!headerChrMatcher.find()) {
                    throw new Exception("Error parsing reference genome fasta header: " + line);
                }

                String chr = headerChrMatcher.group(1);

                Matcher headerMatcher = FASTA_HEADER_PATTERN.matcher(line);

                if (!headerMatcher.matches()) {
                    throw new Exception("Error parsing reference genome fasta header: " + line);
                }
//				System.out.println("----");
//				System.out.println(line);

                long lenght = Long.parseLong(headerMatcher.group(2));

                if (chr.equals("Y") && lenght == 59034049) {
                    lenght = 59373566;
                }

//				System.out.println(chr + " " + lenght);
                pos = 0;
                currentSeq = new LargeByteArray(lenght);
                chromosomes.put(chr, currentSeq);

//				System.out.println(currentSeq.getSize());
            } else {
                for (int n = 0; n < line.length(); ++n) {
                    currentSeq.setQuick(pos++, (byte) line.charAt(n));
                }
            }
        }

    }

    public ReferenceGenomeFasta(File[] fastaFiles) throws IOException, Exception {

        chromosomes = new LinkedHashMap<String, LargeByteArray>(32);
        for (File f : fastaFiles) {
            System.out.println(f.getAbsolutePath());
            BufferedReader fasteReader = new BufferedReader(new FileReader(f));

            String line;
            LargeByteArray currentSeq = null;
            long pos = 0;
            while ((line = fasteReader.readLine()) != null) {

                if (line.length() == 0) {
                    continue;
                }

                if (line.charAt(0) == '>') {

                    //				System.out.println("Previous chr bases read: " + pos);
                    //                                System.out.println(line);
                    Matcher headerChrMatcher = FASTA_HEADER_CHR_PATTERN.matcher(line);

                    if (!headerChrMatcher.find()) {
                        throw new Exception("Error parsing reference genome fasta header: " + line);
                    }

                    String chr = headerChrMatcher.group(1);

                    Matcher headerMatcher = FASTA_HEADER_PATTERN.matcher(line);

                    if (!headerMatcher.matches()) {
                        headerMatcher = FASTA_HEADER_PATTERN_ENSEMBL.matcher(line);
                        if (!headerMatcher.matches()) {
                            throw new Exception("Error parsing reference genome fasta header: " + line);
                        }
                    }
                    //				System.out.println("----");
                    //				System.out.println(line);
//                    System.out.println(headerMatcher.group(2));
                    long lenght = Long.parseLong(headerMatcher.group(2));

                    if (chr.equals("Y") && lenght == 59034049) {
                        lenght = 59373566;
                    }
                    
                    if (chr.equals("Y") && lenght == 56887902) {
                        lenght = 57227415;
                    }
                    //				System.out.println(chr + " " + lenght);
                    pos = 0;
                    currentSeq = new LargeByteArray(lenght);
                    chromosomes.put(chr, currentSeq);

                    //				System.out.println(currentSeq.getSize());
                } else {
                    for (int n = 0; n < line.length(); ++n) {
                        currentSeq.setQuick(pos++, (byte) line.charAt(n));
                    }
                }
            }
        }
    }

    public char getNucleotide(String chr, long pos) throws Exception {
        LargeByteArray chrNucleotides = chromosomes.get(chr);
        if (chrNucleotides == null) {
            throw new Exception("Chr " + chr + " not found in reference fasta");
        }

        if (pos > chrNucleotides.getSize()) {
            throw new Exception("Chr " + chr + " is shorter than < " + pos);
        }

        return (char) chrNucleotides.getQuick(pos - 1);
    }

    public StringBuilder getNucleotides(String chr, long posStart, long posStop) throws Exception {
        LargeByteArray chrNucleotides = chromosomes.get(chr);

        if (chrNucleotides == null) {
            throw new Exception("Chr " + chr + " not found in reference fasta");
        }

        if (posStart > chrNucleotides.getSize() || posStop > chrNucleotides.getSize()) {
            throw new Exception("Chr " + chr + " is shorter than < " + posStart + " or " + posStop);
        }

        StringBuilder nucleotides = new StringBuilder();

        for (long i = posStart - 1; i < posStop; i++) {
            nucleotides.append((char) chrNucleotides.getQuick(i));
        }

        return nucleotides;
    }

    public Set<String> getChromosomes() {
        return Collections.unmodifiableSet(chromosomes.keySet());
    }

    public boolean loadedChr(String chr) {
        return chromosomes.containsKey(chr);
    }
}
