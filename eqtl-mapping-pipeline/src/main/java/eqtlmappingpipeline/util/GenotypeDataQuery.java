/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package eqtlmappingpipeline.util;

import umcg.genetica.io.text.TextFile;
import umcg.genetica.io.trityper.SNP;
import umcg.genetica.io.trityper.SNPLoader;
import umcg.genetica.io.trityper.TriTyperGenotypeData;
import umcg.genetica.io.trityper.util.BaseAnnot;

/**
 *
 * @author harmjan
 */
public class GenotypeDataQuery {

    public void getSNPsInRegion(String dataset, int chr, int positionA, int positionB) {
        if (positionB < positionA) {
            int tmp = positionA;
            positionA = positionB;
            positionB = tmp;
        }

        try {
            TriTyperGenotypeData tgd = new TriTyperGenotypeData();
            tgd.load(dataset);

            String[] snps = tgd.getSNPs();
            int pos = 0;
            for (String s : snps) {
                int snpchr = tgd.getChr(pos);
                int snppos = tgd.getChrPos(pos);
                if (snpchr == chr && snppos >= positionA && snppos <= positionB) {
                    System.out.println(s + "\t" + chr + "\t" + snppos);
                }
                pos++;
            }
        } catch (Exception e) {
            e.printStackTrace();
        }


    }

    public void getSNPMAF(String dataset, String snpfile) {
        try {
            TriTyperGenotypeData tgd = new TriTyperGenotypeData();
            tgd.load(dataset);

            TextFile tf = new TextFile(snpfile, TextFile.R);
            String[] snpstoquery = tf.readAsArray();
            tf.close();

            SNPLoader loader = tgd.createSNPLoader();

            for (String snp : snpstoquery) {
                Integer snpid = tgd.getSnpToSNPId().get(snp);
                if (snpid == -9) {
                    System.out.println("SNP " + snp + " does not exist in dataset");
                } else {

                    SNP snpObj = tgd.getSNPObject(snpid);

                    loader.loadGenotypes(snpObj);

                    System.out.println(snp + "\t" + BaseAnnot.toString(snpObj.getAlleles()[0])
                            + "/" + BaseAnnot.toString(snpObj.getAlleles()[1])
                            + "\t" + BaseAnnot.toString(snpObj.getMinorAllele())
                            + "\t" + snpObj.getMAF());
                }

            }
        } catch (Exception e) {
            e.printStackTrace();
        }


    }

    public void getSNPStatsForAllSNPs(String dataset) {
        try {
            TriTyperGenotypeData tgd = new TriTyperGenotypeData();
            tgd.load(dataset);

            SNPLoader loader = tgd.createSNPLoader();


            String[] snpstoquery = tgd.getSNPs();

            if (loader.hasDosageInformation()) {
                System.out.println("SNP\tAlleles\tMinorAllele\tMAF\tDosageMAF\tCR\tHWEP\tNumAA\tNumAB\tNumBB\tTotalCalled\tNotCalled");
            } else {
                System.out.println("SNP\tAlleles\tMinorAllele\tMAF\tCR\tHWEP\tNumAA\tNumAB\tNumBB\tTotalCalled\tNotCalled");
            }


            for (String snp : snpstoquery) {
                Integer snpid = tgd.getSnpToSNPId().get(snp);
                if (snpid == -9) {
                    System.out.println("SNP " + snp + " does not exist in dataset");
                } else {

                    SNP snpObj = tgd.getSNPObject(snpid);

                    loader.loadGenotypes(snpObj);

                    if (loader.hasDosageInformation()) {
                        loader.loadDosage(snpObj);
                    }

                    int nrGenotypes = 0;
                    int numNULL = 0;
                    int numAA = 0;
                    int numAB = 0;
                    int numBB = 0;


                    byte[] genotypes = snpObj.getGenotypes();
                    for (int g = 0; g < genotypes.length; g++) {
                        switch (genotypes[g]) {
                            case -1:
                                numNULL++;
                                break;
                            case 0:
                                numAA++;
                                nrGenotypes++;
                                break;
                            case 1:
                                numAB++;
                                nrGenotypes++;
                                break;
                            case 2:
                                numBB++;
                                nrGenotypes++;
                                break;
                        }
                    }

                    if (loader.hasDosageInformation()) {
                        System.out.println(snp + "\t" + BaseAnnot.toString(snpObj.getAlleles()[0])
                                + "/" + BaseAnnot.toString(snpObj.getAlleles()[1])
                                + "\t" + BaseAnnot.toString(snpObj.getMinorAllele())
                                + "\t" + snpObj.getMAF()
                                + "\t" + snpObj.getDosageMAF()
                                + "\t" + snpObj.getCR()
                                + "\t" + snpObj.getHWEP()
                                + "\t" + numAA
                                + "\t" + numAB
                                + "\t" + numBB
                                + "\t" + nrGenotypes
                                + "\t" + numNULL);
                    } else {
                        System.out.println(snp + "\t" + BaseAnnot.toString(snpObj.getAlleles()[0])
                                + "/" + BaseAnnot.toString(snpObj.getAlleles()[1])
                                + "\t" + BaseAnnot.toString(snpObj.getMinorAllele())
                                + "\t" + snpObj.getMAF()
                                + "\t" + snpObj.getCR()
                                + "\t" + snpObj.getHWEP()
                                + "\t" + numAA
                                + "\t" + numAB
                                + "\t" + numBB
                                + "\t" + nrGenotypes
                                + "\t" + numNULL);
                    }

                }

            }
        } catch (Exception e) {
            e.printStackTrace();
        }


    }

    void getSNPStatsForAllSNPs(String dataset, String in2) {
        try {
            TriTyperGenotypeData tgd = new TriTyperGenotypeData();
            tgd.load(dataset);

            SNPLoader loader = tgd.createSNPLoader();


            String[] snpstoquery = new String[]{in2};

            if (loader.hasDosageInformation()) {
                System.out.println("SNP\tAlleles\tMinorAllele\tMAF\tDosageMAF\tCR\tHWEP\tNumAA\tNumAB\tNumBB\tTotalCalled\tNotCalled");
            } else {
                System.out.println("SNP\tAlleles\tMinorAllele\tMAF\tCR\tHWEP\tNumAA\tNumAB\tNumBB\tTotalCalled\tNotCalled");
            }


            for (String snp : snpstoquery) {
                Integer snpid = tgd.getSnpToSNPId().get(snp);
                if (snpid == -9) {
                    System.out.println("SNP " + snp + " does not exist in dataset");
                } else {

                    SNP snpObj = tgd.getSNPObject(snpid);

                    loader.loadGenotypes(snpObj);

                    if (loader.hasDosageInformation()) {
                        loader.loadDosage(snpObj);
                    }

                    int nrGenotypes = 0;
                    int numNULL = 0;
                    int numAA = 0;
                    int numAB = 0;
                    int numBB = 0;


                    byte[] genotypes = snpObj.getGenotypes();
                    for (int g = 0; g < genotypes.length; g++) {
                        switch (genotypes[g]) {
                            case -1:
                                numNULL++;
                                break;
                            case 0:
                                numAA++;
                                nrGenotypes++;
                                break;
                            case 1:
                                numAB++;
                                nrGenotypes++;
                                break;
                            case 2:
                                numBB++;
                                nrGenotypes++;
                                break;
                        }
                    }

                    if (loader.hasDosageInformation()) {
                        System.out.println(snp + "\t" + BaseAnnot.toString(snpObj.getAlleles()[0])
                                + "/" + BaseAnnot.toString(snpObj.getAlleles()[1])
                                + "\t" + BaseAnnot.toString(snpObj.getMinorAllele())
                                + "\t" + snpObj.getMAF()
                                + "\t" + snpObj.getDosageMAF()
                                + "\t" + snpObj.getCR()
                                + "\t" + snpObj.getHWEP()
                                + "\t" + numAA
                                + "\t" + numAB
                                + "\t" + numBB
                                + "\t" + nrGenotypes
                                + "\t" + numNULL);
                    } else {
                        System.out.println(snp + "\t" + BaseAnnot.toString(snpObj.getAlleles()[0])
                                + "/" + BaseAnnot.toString(snpObj.getAlleles()[1])
                                + "\t" + BaseAnnot.toString(snpObj.getMinorAllele())
                                + "\t" + snpObj.getMAF()
                                + "\t" + snpObj.getCR()
                                + "\t" + snpObj.getHWEP()
                                + "\t" + numAA
                                + "\t" + numAB
                                + "\t" + numBB
                                + "\t" + nrGenotypes
                                + "\t" + numNULL);
                    }

                }

            }
        } catch (Exception e) {
            e.printStackTrace();
        }


    }
}
