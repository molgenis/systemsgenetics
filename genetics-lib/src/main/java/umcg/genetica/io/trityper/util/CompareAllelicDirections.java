/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package umcg.genetica.io.trityper.util;

import umcg.genetica.io.trityper.SNP;

/**
 *
 * @author harmjan
 */
public class CompareAllelicDirections {

    public static Boolean[] compare(SNP[] snps) {
        int firstDatasetToPassQC = -1;
        byte[] firstDatasetPassinQCAlleles = null;
        Boolean[] flipAlleleOutput = new Boolean[snps.length];
        byte firstminor = -1;

        for (int d = 0; d < snps.length; d++) {
            SNP dSNP = snps[d];

            if (dSNP != null) {
                // check if the alleles are identical with previuously loaded SNP...
                if (firstDatasetToPassQC == -1) {
                    firstDatasetToPassQC = d;
                    firstDatasetPassinQCAlleles = dSNP.getAlleles();
                    byte minor = dSNP.getMinorAllele();
                    if (firstDatasetPassinQCAlleles[1] == minor) {
                        flipAlleleOutput[d] = false;
                    } else {
                        flipAlleleOutput[d] = true;
                    }

                    firstminor = minor;

                } else {

                    byte[] allelesToCompare = dSNP.getAlleles();
                    int nrAllelesIdentical = 0;
                    byte minor = dSNP.getMinorAllele();

                    boolean flipalleles = false;
                    int minorAlleleNum = 0;
                    if (allelesToCompare[0] != minor) {
                        minorAlleleNum = 1;
                    }

                    for (int a = 0; a < 2; a++) {
                        for (int b = 0; b < 2; b++) {
                            if (firstDatasetPassinQCAlleles[a] == allelesToCompare[b]) {
                                nrAllelesIdentical++;
                            }
                        }
                    }

                    if (nrAllelesIdentical != 2) {
                        //Alleles are different, take complimentary:
                        allelesToCompare = BaseAnnot.convertToComplementaryAlleles(allelesToCompare);
                        minor = BaseAnnot.getComplement(minor);
                    }

                    nrAllelesIdentical = 0;

                    for (int a = 0; a < 2; a++) {
                        for (int b = 0; b < 2; b++) {
                            if (firstDatasetPassinQCAlleles[a] == allelesToCompare[b]) {
                                nrAllelesIdentical++;
                            }
                        }
                    }

                    if (nrAllelesIdentical != 2) {
                        return null;
                    } else {
                        if (minor != firstminor) {
                            // error or warning or whatever
                            //                        System.out.println("WARNING: minor allele is different for identical SNP: "+dSNP.getName() + ", probably due to high MAF.\nWill conform to allelic direction of dataset: "+m_gg[firstDatasetToPassQC].getSettings().name);
                            //                        double[] allelefreq = dSNP.getAlleleFreq();
                            //                        byte[] origAlleles = snps[firstDatasetToPassQC].getAlleles();
                            //                        double[] origAlleleFreq = snps[firstDatasetToPassQC].getAlleleFreq();
                            //                        System.out.println("Reference MAF:"+snps[firstDatasetToPassQC].getMAF()+"\tAssessed MAF:"+dSNP.getMAF());
                            //                        for(int i=0; i<2; i++){
                            //                            System.out.println("ref ds: "+m_gg[firstDatasetToPassQC].getSettings().name+"\t"+BaseAnnot.toString(origAlleles[i])+"\t("+origAlleleFreq[i]+")\tAssessed: "+m_gg[d].getSettings().name+"\t"+BaseAnnot.toString(allelesToCompare[i])+"\t("+allelefreq[i]+")");
                            //                        }
                            //                        System.out.println("");
                            // take the orientation of the first dataset..., which is dataset
                            if (minorAlleleNum == 1) {
                                flipalleles = true;
                            } else {
                                flipalleles = false;
                            }

                        } else {
                            if (allelesToCompare[0] == minor) {
                                flipalleles = true;
                            } else {
                                flipalleles = false;
                            }
                        }

                        flipAlleleOutput[d] = flipalleles;
                    }



                }
            }
        }

        return flipAlleleOutput;
    }
}
