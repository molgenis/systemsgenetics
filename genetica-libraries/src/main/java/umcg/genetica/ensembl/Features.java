/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package umcg.genetica.ensembl;

import umcg.genetica.containers.Transcript;
import umcg.genetica.containers.Exon;
import umcg.genetica.containers.Gene;
import umcg.genetica.containers.Chromosome;
import java.io.IOException;
import java.util.HashMap;
import umcg.genetica.io.text.TextFile;

/**
 *
 * @author harmjan
 */
public class Features {

    private HashMap<String, Chromosome> chromosomeHash = new HashMap<String, Chromosome>();
    private HashMap<String, Gene> geneHash = new HashMap<String, Gene>();
    private HashMap<String, Transcript> transcriptHash = new HashMap<String, Transcript>();
    private HashMap<String, Exon> exonHash = new HashMap<String, Exon>();

    /**
     * @return the chromosomeHash
     */
    public HashMap<String, Chromosome> getChromosomeHash() {
        return chromosomeHash;
    }

    /**
     * @param chromosomeHash the chromosomeHash to set
     */
    public void setChromosomeHash(HashMap<String, Chromosome> chromosomeHash) {
        this.chromosomeHash = chromosomeHash;
    }

    /**
     * @return the geneHash
     */
    public HashMap<String, Gene> getGeneHash() {
        return geneHash;
    }

    /**
     * @param geneHash the geneHash to set
     */
    public void setGeneHash(HashMap<String, Gene> geneHash) {
        this.geneHash = geneHash;
    }

    /**
     * @return the transcriptHash
     */
    public HashMap<String, Transcript> getTranscriptHash() {
        return transcriptHash;
    }

    /**
     * @param transcriptHash the transcriptHash to set
     */
    public void setTranscriptHash(HashMap<String, Transcript> transcriptHash) {
        this.transcriptHash = transcriptHash;
    }

    /**
     * @return the exonHash
     */
    public HashMap<String, Exon> getExonHash() {
        return exonHash;
    }

    /**
     * @param exonHash the exonHash to set
     */
    public void setExonHash(HashMap<String, Exon> exonHash) {
        this.exonHash = exonHash;
    }


    /*
     * 0	Ensembl Gene ID 1	Ensembl Transcript ID 2	Ensembl Protein ID
     *
     * 3	Chromosome Name 4	Gene Start (bp) 5	Gene End (bp) 6	Transcript Start (bp) 7	Transcript End (bp) 8	Strand 9	Associated Gene Name 10	Associated Gene DB
     * 11	5' UTR Start 12	5' UTR End 13	3' UTR Start 14	3' UTR End 15	CDS Length 16	Transcript count 17	Description 18	Gene Biotype 19	Ensembl Exon ID 20	Exon
     * Chr Start (bp) 21	Exon Chr End (bp) 22	Constitutive Exon 23	Exon Rank in Transcript 24	phase 25	cDNA coding start 26	Genomic coding start 27	cDNA coding
     * end 28	Genomic coding end 29	CDS End 30	CDS Start
     */
    public void loadAnnotation(String fileName) throws IOException {
        System.out.println("Loading sequence feature annotation from: " + fileName);

        TextFile in = new TextFile(fileName, TextFile.R);

        String[] linesplit = in.readLineElems(TextFile.tab);
        HashMap<String, Integer> elementToColId = new HashMap<String, Integer>();
        for (int i = 0; i < linesplit.length; i++) {
            System.out.println(linesplit[i] + "\t" + i);
            elementToColId.put(linesplit[i], i);
        }


        int lnCounter = 0;
        String[] elems = in.readLineElems(TextFile.tab);
        while (elems != null) {

            /*
             * Ensembl Gene ID Ensembl Transcript ID Ensembl Protein ID Chromosome Name Gene Start (bp) Gene End (bp) Transcript Start (bp) Transcript End (bp)
             * Strand Associated Gene Name Gene Biotype Description Transcript count CDS Length 3' UTR End 3' UTR Start 5' UTR End 5' UTR Start Ensembl Exon ID
             * Exon Chr Start (bp) Exon Chr End (bp) Exon Rank in Transcript phase Associated Gene DB
             *
             *
             *
             *
             *
             */
            if (elems.length > 1) {
                String gene = elems[elementToColId.get("Ensembl Gene ID")];
                String transcript = elems[elementToColId.get("Ensembl Transcript ID")];
                String protein = elems[elementToColId.get("Ensembl Protein ID")];
                String chromosome = elems[elementToColId.get("Chromosome Name")];
                Integer genestart = Integer.parseInt(elems[elementToColId.get("Gene Start (bp)")]);
                Integer genestop = Integer.parseInt(elems[elementToColId.get("Gene End (bp)")]);
                Integer transcriptstart = Integer.parseInt(elems[elementToColId.get("Transcript Start (bp)")]);
                Integer transcriptstop = Integer.parseInt(elems[elementToColId.get("Transcript End (bp)")]);
                Integer strand = Integer.parseInt(elems[elementToColId.get("Strand")]);
                String hgnc = elems[elementToColId.get("Associated Gene Name")];
                //                String db               = elems[10];
                //                Integer utr5start       = Integer.parseInt(elems[11]);
                //                Integer utr5end         = Integer.parseInt(elems[12]);
                //                Integer utr3start       = Integer.parseInt(elems[13]);
                //                Integer utr3end         = Integer.parseInt(elems[14]);
                //                Integer cdslen          = Integer.parseInt(elems[15]);
                //                Integer numTranscripts  = Integer.parseInt(elems[16]);
                //                String Desc  = Integer.parseInt(elems[17]);
                //                String biotype          = elems[18];
                String exon = elems[elementToColId.get("Ensembl Exon ID")];
                Integer exonstart = Integer.parseInt(elems[elementToColId.get("Exon Chr Start (bp)")]);
                Integer exonend = Integer.parseInt(elems[elementToColId.get("Exon Chr End (bp)")]);
//                    String exonconstitutive = elems[23];
                Integer exonrankintranscript = Integer.parseInt(elems[elementToColId.get("Exon Rank in Transcript")]);
                // String phase                   = elems[24];
//                    Integer cDNACodingStart         = Integer.parseInt(elems[25]);
//                    Integer genomicCodingStart      = Integer.parseInt(elems[26]);
//                    Integer cDNACodingEnd           = Integer.parseInt(elems[27]);
//                    Integer genomicCodingEnd        = Integer.parseInt(elems[28]);
//                    Integer CDSStart                = Integer.parseInt(elems[29]);
//                    Integer CDSEnd                  = Integer.parseInt(elems[30]);

                Chromosome currChr = null;
                Gene currGen = null;
                Transcript currTra = null;
                Exon currExo = null;

                if (chromosome.trim().length() > 0) {
                    if (chromosomeHash.get(chromosome) == null) {
                        Chromosome tmpChr = new Chromosome(chromosome);
                        tmpChr.setName(chromosome);
                        chromosomeHash.put(chromosome, tmpChr);
                        System.out.println("Adding chromosome: " + chromosome);
                    }
                    currChr = chromosomeHash.get(chromosome);
                }

                if (gene.trim().length() > 0) {
                    if (geneHash.get(gene) == null) {
                        Gene tmpGen = new Gene();
                        tmpGen.setName(gene);
                        tmpGen.setStart(genestart);
                        tmpGen.setEnd(genestop);
                        tmpGen.setStrand(strand);
                        tmpGen.setAnnotation(hgnc);
                        geneHash.put(gene, tmpGen);
                    }
                    currGen = geneHash.get(gene);
                }

                if (transcript.trim().length() > 0) {
                    if (transcriptHash.get(transcript) == null) {
                        Transcript tmpTra = new Transcript();
                        tmpTra.setName(transcript);
                        tmpTra.setStart(transcriptstart);
                        tmpTra.setEnd(transcriptstop);
                        tmpTra.setStrand(strand);
                        tmpTra.setProtein(protein);
                        transcriptHash.put(transcript, tmpTra);
                    }
                    currTra = transcriptHash.get(transcript);
                }

                if (exon.trim().length() > 0) {
                    if (exonHash.get(exon) == null) {
                        Exon tmpExo = new Exon();
                        tmpExo.setName(exon);
                        tmpExo.setStart(exonstart);
                        tmpExo.setEnd(exonend);

                        tmpExo.setStrand(strand);
                        exonHash.put(exon, tmpExo);
                    }
                    currExo = exonHash.get(exon);
                }



                if (currGen != null) {
                    currChr.addGene(currGen);
                    currGen.setParentChromosome(currChr);
                }

                if (currGen != null && currTra != null) {
                    currGen.addTranscript(currTra);
                    currTra.setParentGene(currGen);
                    currTra.setParentChromosome(currChr);
                }

                if (currTra != null && currExo != null) {
                    currTra.addExon(currExo);
                    currExo.setParentTranscript(currTra);
                    currExo.setParentChromosome(currChr);
                    currTra.setExonRank(currExo, exonrankintranscript);
                }

            }

            if (lnCounter % 100000 == 0) {
                System.out.print(".");
            }
            elems = in.readLineElems(TextFile.tab);
            lnCounter++;
        }
        System.out.println("\tDone.");

        in.close();
        System.out.println("Loaded " + chromosomeHash.size() + " chromosomes, " + geneHash.size() + " genes, " + transcriptHash.size() + " transcripts, " + exonHash.size() + " exons.");


    }

    public Exon getExon(String match) {
        return exonHash.get(match);
    }

    public Transcript getTranscript(String match) {
        return transcriptHash.get(match);
    }
}
