package mbqtl;

import cern.colt.matrix.tdouble.DoubleMatrix1D;
import mbqtl.data.Dataset;
import mbqtl.stat.RankArray;
import mbqtl.vcf.VCFTabix;
import mbqtl.vcf.VCFVariant;
import org.apache.commons.math3.linear.SingularMatrixException;
import org.apache.commons.math3.stat.regression.OLSMultipleLinearRegression;
import umcg.genetica.console.ProgressBar;
import umcg.genetica.containers.Pair;
import umcg.genetica.enums.Chromosome;
import umcg.genetica.features.Feature;
import umcg.genetica.math.matrix2.DoubleMatrixDataset;
import umcg.genetica.math.stats.Correlation;
import umcg.genetica.math.stats.Descriptives;
import umcg.genetica.math.stats.VIF;

import java.io.IOException;
import java.util.*;

public class QTLRegression extends QTLAnalysis {

    public QTLRegression(String vcfFile, int chromosome, String linkfile, String snpLimitFile, String geneLimitFile, String snpGeneLimitFile, String geneExpressionDataFile, String geneAnnotationFile, String outputPrefix) throws IOException {
        super(vcfFile, chromosome, linkfile, snpLimitFile, geneLimitFile, snpGeneLimitFile, geneExpressionDataFile, geneAnnotationFile, outputPrefix);
    }

    private boolean debug = false;
    private boolean rankData = true;
    private boolean replaceMissingGenotypes = true;
    private int cisWindow = 1000000;

    public void setCisWindow(int cisWindow) {
        this.cisWindow = cisWindow;
    }

    public void setRankData(boolean rankData) {
        this.rankData = rankData;
    }

    public void setReplaceMissingGenotypes(boolean replaceMissingGenotypes) {
        this.replaceMissingGenotypes = replaceMissingGenotypes;
    }


    class GeneObj implements Comparable<GeneObj> {
        int chr;
        int pos;
        int id;
        String name;

        @Override
        public boolean equals(Object o) {
            if (this == o) return true;
            if (o == null || getClass() != o.getClass()) return false;
            GeneObj geneObj = (GeneObj) o;
            return chr == geneObj.chr && pos == geneObj.pos && Objects.equals(name, geneObj.name);
        }

        @Override
        public int hashCode() {
            return Objects.hash(chr, pos);
        }

        @Override
        public int compareTo(GeneObj geneObj) {
            if (this.equals(geneObj)) {
                return 0;
            } else {
                if (this.chr > geneObj.chr) {
                    return 1;
                } else if (this.chr < geneObj.chr) {
                    return -1;
                } else {
                    if (this.pos > geneObj.pos) {
                        return 1;
                    } else {
                        return -1;
                    }
                }
            }
        }
    }


    public void run() throws IOException {

        System.out.println();
        System.out.println("-----------------------------------------------------");
        System.out.println("Cis-QTL regression");
        System.out.println("-----------------------------------------------------");
        System.out.println("MAF:\t" + mafthreshold);
        System.out.println("Call-rate:\t" + callratethreshold);
        System.out.println("HWE-P:\t" + hwepthreshold);
        System.out.println("Min nr datasets:\t" + minNumberOfDatasets);
        System.out.println("Ranking data:\t" + rankData);
//        System.out.println("Re-ranking data:\t");
        System.out.println("Replacing missing genotypes: " + replaceMissingGenotypes);
        System.out.println("Min observations: " + minObservations);
        System.out.println();

        Chromosome chromosomeObj = Chromosome.parseChr("" + chromosome);
        System.out.println("Processing: " + vcfFile);

        // sort the genes by position for easier VCF access
        ArrayList<GeneObj> genes = new ArrayList<>();
        for (int g = 0; g < expressionData.genes.length; g++) {
            String gene = expressionData.genes[g];
            Integer geneAnnotationId = geneAnnotation.getGeneId(gene);
            if (geneAnnotationId != null) {
                int pos = geneAnnotation.getStartPos(geneAnnotationId);
                int chr = geneAnnotation.getChr(geneAnnotationId);
                GeneObj obj = new GeneObj();
                obj.chr = chr;
                obj.pos = pos;
                obj.id = g;
                genes.add(obj);
            }
        }
        Collections.sort(genes);

        // iterate genes
        ProgressBar pb = new ProgressBar(genes.size(), "Regressing QTLs...");
//		for (int g = 0; g < expressionData.genes.length; g++) {
        int prevchr = -1;
        VCFTabix tabix = null;
        for (int gid = 0; gid < genes.size(); gid++) {
            GeneObj geneObj = genes.get(gid);
            int geneID = geneObj.id;
            String gene = expressionData.genes[geneID];
            HashSet<String> snpToRemove = snpGeneLimitSet.get(gene);
            if (snpToRemove != null) {
                Integer geneAnnotationId = geneAnnotation.getGeneId(gene);
                int geneChr = geneAnnotation.getChr(geneAnnotationId);
                chromosomeObj = Chromosome.parseChr("" + geneChr);
                double[] expData = expressionData.data[geneID];

                // define CIS window
                int pos = geneAnnotation.getStartPos(geneAnnotationId);
                String geneSymbol = geneAnnotation.getSymbol(geneAnnotationId);
                int start = pos - cisWindow;
                if (start < 0) {
                    start = 0;
                }
                int stop = pos + cisWindow;
                Feature cisRegion = new Feature(chromosomeObj, start, stop);

                if (geneChr != prevchr && origvcfFile.contains("CHR")) {
                    if (tabix != null) {
                        tabix.close();
                    }
                    System.out.println();
                    System.out.println("Changing chromosome to " + geneChr);
                    String actualVCFFile = origvcfFile.replace("CHR", "" + geneChr);
                    updateDatasets(actualVCFFile);
                    tabix = new VCFTabix(actualVCFFile);
                    prevchr = geneChr;
                }
                // attempt to open VCF file for chromosome

                // sample order in another VCF can be different, so re-initialize dataset definitions

                // load variants for gene
                ArrayList<VCFVariant> variants = new ArrayList<>();
                HashSet<String> loadedVariants = new HashSet<>();
                if (snpAnnotation == null) {
                    // TABIX doesn't allow for lookup of variants by ID, so we have to iterate the cisregion of the gene
                    Iterator<VCFVariant> snpIterator = tabix.getVariants(cisRegion, genotypeSamplesToInclude, snpToRemove);
                    while (snpIterator.hasNext()) {
                        VCFVariant variant = snpIterator.next();
                        if (variant != null) {
                            String variantId = variant.getId();
                            if (snpToRemove.contains(variantId) && !loadedVariants.contains(variantId)) {
                                variants.add(variant);
                                loadedVariants.add(variantId);
                            }
                        }
                    }
                } else {
                    for (String variantName : snpToRemove) {
                        Integer id = snpAnnotation.getId(variantName);
                        if (id == null) {
                            System.out.println("Error: SNP annotation loaded, but no annotation for " + variantName);
                            System.exit(0);
                        }
                        int snppos = snpAnnotation.getPos(id);
                        cisRegion.setStart(snppos - 1);
                        cisRegion.setStop(snppos + 1);
                        Iterator<VCFVariant> snpIterator = tabix.getVariants(cisRegion, genotypeSamplesToInclude, snpToRemove);
                        while (snpIterator.hasNext()) {
                            VCFVariant variant = snpIterator.next();
                            if (variant != null) {
                                String variantId = variant.getId();
                                if (snpToRemove.contains(variantId) && !loadedVariants.contains(variantId)) {
                                    variants.add(variant);
                                    loadedVariants.add(variantId);
                                }
                            }
                        }
                    }
                }

//                System.out.println(gene + "\t" + variants.size());
//                int datasetsRemoved = 0;
                for (int d = 0; d < datasets.length; d++) {
                    Dataset thisDataset = datasets[d];

                    double[] datasetExp = thisDataset.select(expData, thisDataset.expressionIds);
                    double[] datasetExpRaw = thisDataset.select(expData, thisDataset.expressionIds);
                    double[] datasetExpRanked = datasetExp;
                    if (rankData) {
                        RankArray ranker = new RankArray();
                        datasetExpRanked = ranker.rank(datasetExp, true); // does this work with NaNs? answer: no
                    }

                    double meanY = JSci.maths.ArrayMath.mean(datasetExpRaw);
                    double varianceY = JSci.maths.ArrayMath.variance(datasetExpRaw);

                    // get genotypes
                    ArrayList<double[]> xs = new ArrayList<>();
                    ArrayList<String> includedVariants = new ArrayList<>();
                    for (int v = 0; v < variants.size(); v++) {
                        VCFVariant variant = variants.get(v);
                        final double[] genotypes = getGenotype(variant.getGenotypesAsByteVector());
                        final double[] dosages = getDosage(variant.getDosage());
                        double[] dosagesForDataset = thisDataset.select(dosages, thisDataset.genotypeIds); // select required dosages
                        double[] genotypesForDataset = thisDataset.select(genotypes, thisDataset.genotypeIds); // select required genotype IDs

                        VariantQCObj qcobj = checkVariant(genotypesForDataset);
                        if (qcobj.passqc) {
                            // replace missing genotypes
                            // only replace missing genotypes on variants that pass the qc thresholds
//                            for(int z=0;z<dosagesForDataset.length;z++){
//                                System.out.println(dosagesForDataset[z]+"\t"+genotypesForDataset[z]);
//                            }
                            double meanDosage = Util.meanGenotype(dosagesForDataset);
                            double meanGenotype = Util.meanGenotype(genotypesForDataset);
//                            System.out.println(meanDosage);
//                            System.out.println(meanGenotype);
//                            System.exit(-1);
                            for (int i = 0; i < dosagesForDataset.length; i++) {
                                if (genotypesForDataset[i] == -1) {
                                    genotypesForDataset[i] = meanGenotype;
                                    dosagesForDataset[i] = meanDosage;
                                }
                            }
//                            double[] dosagesForDatasetCenterScale = Util.centerScale(dosagesForDataset);
//                            xs.add(dosagesForDatasetCenterScale);
                            xs.add(dosagesForDataset);
                            includedVariants.add(variant.getId());
                        }
                    }

                    if (!includedVariants.isEmpty()) {
                        // TODO: check for missingness in expression data
//                        double[] expdataForDataset = Util.centerScale(datasetExpRanked);
                        double[] expdataForDataset = datasetExpRanked; // Util.centerScale(datasetExpRanked);

                        int nrIndividuals = xs.get(0).length;
                        DoubleMatrixDataset<String, String> xcovars = new DoubleMatrixDataset<>(xs.size(), xs.get(0).length);
                        for (int i = 0; i < nrIndividuals; i++) {
                            xcovars.getHashCols().put("Ind" + i, i);
                        }

                        for (int i = 0; i < xs.size(); i++) {
                            double[] vals = xs.get(i);
                            if (debug) {
                                double r = Correlation.correlate(expdataForDataset, vals);
                                System.out.println(d + "\t" + r);
                            }
                            xcovars.getRow(i).assign(vals);
                            String snpid = includedVariants.get(i);
                            xcovars.getHashRows().put(snpid, i);
                        }

                        // transpose (samples should be on rows)
                        try {
                            xcovars = xcovars.viewDice();
                        } catch (IllegalArgumentException e) {
                            System.out.println("Issue with genotype matrix: ");
                            System.out.println(xs.size() + " variants for gene " + gene);
                            System.out.println(includedVariants.size() + " SNP IDS");
                            for (int v = 0; v < includedVariants.size(); v++) {
                                System.out.println(includedVariants.get(v));
                            }
                            System.out.println(xcovars.rows() + "\t" + xcovars.columns());
                            for (int r = 0; r < xcovars.rows(); r++) {
                                System.out.println(r + "\t" + xcovars.getRowObjects().get(r));
                            }
                            for (int c = 0; c < xcovars.columns(); c++) {
                                System.out.println(c + "\t" + xcovars.getColObjects().get(c));
                            }
                            System.exit(-1);
                        }
                        // check whether there are more predictors than data rows
                        if (xcovars.rows() < xcovars.columns()) {
                            // remove the rows with lowest variance
                            int toRemove = (xcovars.columns() - xcovars.rows()) + 1;
                            System.out.println("\nWarning: " + thisDataset.name + " has more predictors than datapoints for gene " + gene + ": " + xcovars.rows() + "x" + xcovars.columns() + " removing " + toRemove + " lowest variance covars");
                            xcovars = removeCovarWithLowestVariance(xcovars, toRemove);
                            System.out.println("\nWarning: " + thisDataset.name + " gene had few covars for " + gene + ". Remaining covars: " + xcovars.rows() + "x" + xcovars.columns());
                        }

                        try {
                            // prevent aliasing; correct for variance inflation.
                            if (xcovars.columns() > 1) {
                                VIF vif = new VIF();
                                int prevCovars = xcovars.columns();
                                xcovars = vif.vifCorrect(xcovars, (1 - 1E-4));
                                int currentCovars = xcovars.columns();
//                            if (logout != null) {
//                                logout.writeln(gene + "\t had " + prevCovars + " before VIF, and " + currentCovars + " after.");
//                            }
                            }

                            // initialize OLS
                            OLSMultipleLinearRegression ols = new OLSMultipleLinearRegression();
                            boolean singular = true;
                            double[] residuals = null;
                            double[][] covars = xcovars.getMatrixAs2dDoubleArray();
                            double rsq = 0;
                            while (singular) {
                                if (covars[0].length > 0) {
                                    ols.newSampleData(expdataForDataset, covars);
                                    try {
                                        // use OLS to determine regression coefficients
                                        residuals = ols.estimateResiduals();
//                                        datasetsRemoved++;
                                        singular = false;
                                        rsq = ols.calculateRSquared(); // I'm assuming this is an appropriate approximation of the explained variance.
                                        // debug: check whether residuals are correlated to genotype?
                                    } catch (SingularMatrixException e) {
                                        // remove lowest variance covariate
                                        // covars has samples on rows, covars on cols

                                        System.err.println("WARNING: singular matrix exception when regressing eQTLs for: " + gene + " with " + covars[0].length + " covariates (variants). Removing lowest variance covariate.");

                                        if (covars[0].length > 1) {
                                            covars = removeCovarWithLowestVariance(covars, 1);

                                        } else {
                                            System.err.println("WARNING: could not resolve covariate issue for: " + gene + " keeping original data.");
//                                        if (logout != null) {
//                                            logout.writeln("WARNING: could not resolve covariate issue for: " + gene + " keeping original data.");
//                                        }
                                            singular = false;
                                            rsq = 0;
                                        }
                                    }
                                } else {
                                    // nothing more to do, all covariates have some issue or another
                                    singular = false;
                                }
                            }

                            if (covars[0].length > 0) {
                                if (rsq < 0) {
                                    if (rsq < -1E-9) {
                                        System.out.println("Warning: large negative r-squared: " + rsq + ". MeanY: " + meanY + ", varY: " + varianceY + ", SumSqTotal: " + ols.calculateTotalSumOfSquares() + ", SumSqResid: " + ols.calculateResidualSumOfSquares());
//                                    if (logout != null) {
//                                        logout.writeln("Warning: large negative r-squared: " + rsq + ". MeanY: " + meanY + ", varY: " + varianceY + ", SumSqTotal: " + ols.calculateTotalSumOfSquares() + ", SumSqResid: " + ols.calculateResidualSumOfSquares());
//                                    }
                                    }
                                    rsq = 0d;
                                } else if (rsq > 1) {
                                    System.out.println("Warning: r-squared > 1.0: " + rsq + ". MeanY: " + meanY + ", varY" + varianceY + ", SumSqTotal: " + ols.calculateTotalSumOfSquares() + ", SumSqResid: " + ols.calculateResidualSumOfSquares());
//                                if (logout != null) {
//                                    logout.writeln("Warning: r-squared > 1.0: " + rsq + ". MeanY: " + meanY + ", varY" + varianceY + ", SumSqTotal: " + ols.calculateTotalSumOfSquares() + ", SumSqResid: " + ols.calculateResidualSumOfSquares());
//                                }
                                    rsq = 1d;
                                }
//                            explainedVariancePerEQTLProbe[d][(int) Math.round(rsq * 100d)]++;

//                            if (logout != null) {
//                                SpearmansCorrelation sp = new SpearmansCorrelation();
//
//                                RankDoubleArray rda = new RankDoubleArray();
//                                double[] ry = rda.rank(y);
//                                double[] correlcoeff = new double[xcovars.columns()];
//
////								String[] indsY = currentDataset.getExpressionData().getIndividuals();
////								String[] indsX = currentDataset.getGenotypeData().getIndividuals();
//
////								for (int q = 0; q < y.length; q++) {
////									String outln = indsY[q] + "\t" + indsX[q] + "\t" + indsX[currentDataset.getExpressionToGenotypeIdArray()[q]] + "\t" + y[q] + "\t" + ry[q];
////									for (int xc = 0; xc < xcovars.columns(); xc++) {
////										outln += "\t" + xcovars.getElementQuick(q, xc);
////									}
////									System.out.println(q + "\t" + outln);
////								}
//                                for (int c = 0; c < xcovars.columns(); c++) {
//                                    correlcoeff[c] = sp.correlation(xcovars.getCol(c).toArray(), ry);
//                                }
//                                String logln = gene + "\tNr SNPs: " + xcovars.columns() + "\tMeanY: " + meanY + "\tVarY: " + varianceY + "\trsq: " + rsq + "\tcorrel: " + Strings.concat(correlcoeff, Strings.tab);
//                                logout.writeln(logln);
////								System.out.println(logln);
////								System.exit(-1);
//                            }

                                double meanUpdated = JSci.maths.ArrayMath.mean(residuals);
                                double stdDevRatio = JSci.maths.ArrayMath.standardDeviation(residuals) / Math.sqrt(varianceY);

                                if (!Double.isNaN(meanUpdated) && !Double.isNaN(stdDevRatio) && stdDevRatio > 0) {
                                    for (int s = 0; s < residuals.length; s++) {
                                        residuals[s] -= meanUpdated;
                                        residuals[s] /= stdDevRatio;
                                        residuals[s] += meanY;
                                    }

                                    if (debug) {
                                        double[] vals = xs.get(0);
                                        RankArray ranker = new RankArray();
//                                    double[] resranked = ranker.rank(residuals, true);
                                        double r = Correlation.correlate(residuals, vals);
                                        System.out.println(d + "\tResidual: " + r);
                                    }
                                    // copy the data back to the matrix somehow.
                                    for (int q = 0; q < thisDataset.expressionIds.length; q++) {
                                        int id = thisDataset.expressionIds[q];
                                        expressionData.data[geneID][id] = residuals[q];
                                    }

                                    if (debug) {
//                                    // DEBUG STUFF:
                                        double[] expData2 = expressionData.data[geneID];
                                        double[] datasetExp2 = thisDataset.select(expData2, thisDataset.expressionIds);

                                        VCFVariant variant = variants.get(0);
                                        final double[] genotypes = getGenotype(variant.getGenotypesAsByteVector());
                                        final double[] dosages = getDosage(variant.getDosage());
                                        double[] dosagesForDataset = thisDataset.select(dosages, thisDataset.genotypeIds); // select required dosages
                                        double[] genotypesForDataset = thisDataset.select(genotypes, thisDataset.genotypeIds);
                                        double meanDosage = Util.meanGenotype(dosagesForDataset);
                                        double meanGenotype = Util.meanGenotype(genotypesForDataset);
//                            System.out.println(meanDosage);
//                            System.out.println(meanGenotype);
//                            System.exit(-1);
                                        RankArray ranker = new RankArray();
                                        for (int i = 0; i < dosagesForDataset.length; i++) {
                                            if (genotypesForDataset[i] == -1) {
                                                genotypesForDataset[i] = meanGenotype;
                                                dosagesForDataset[i] = meanDosage;
                                            }
                                        }
                                        dosagesForDataset = Util.centerScale(dosagesForDataset);
                                        datasetExp2 = ranker.rank(datasetExp2, true); // does this work with NaNs? answer: no
                                        datasetExp2 = Util.centerScale(datasetExp2);

                                        double r2 = Correlation.correlate(datasetExp2, dosagesForDataset);
                                        System.out.println(d + "\tResidual replaced/ranked: " + r2);
                                    }
//                                nrEQTLGenesRegressedOut[d]++;
//                                nrEQTLsRegressedOut[d] += xcovars.columns();
                                } else {
                                    String logln = "Error: " + gene + "\tNr SNPs: " + xcovars.columns() + "\tMeanY: " + meanY + "\tVarY: " + varianceY + "\trsq: " + rsq + "\tmeanUpdated: " + meanUpdated + "\tstdevRatio: " + stdDevRatio;
//                                if (logout != null) {
//                                    logout.writeln(logln);
//                                }
                                }
                            }
                        } catch (Exception e) {
                            e.printStackTrace();
                        }
                    }


                }
                //System.out.println(gene + "\tDS removed: " + datasetsRemoved);
            }
            pb.set(gid);

        }

        // save the residualized expressiondata somewhere
        System.out.println();
        System.out.println("Saving residuals: " + outputPrefix + "-QTLsRemoved.txt.gz");
        expressionData.save(outputPrefix + "-QTLsRemoved.txt.gz");
        System.out.println("Done");
        if (tabix != null) {
            tabix.close();
        }
    }


    private double[][] removeCovarWithLowestVariance(double[][] covars, int nrToRemove) {
        int minCovar = -1;
        double minVar = Double.MAX_VALUE;

        int startlen = covars[0].length;
        int endlen = startlen - nrToRemove;
        if (endlen <= 1) {
            return covars;
        }
        while (covars[0].length > endlen) {
            for (int c = 0; c < covars[0].length; c++) {
                double[] col = new double[covars.length];
                for (int r = 0; r < covars[0].length; r++) {
                    col[r] = covars[r][c];
                }
                double var = Descriptives.variance(col);
                if (var < minVar) {
                    minVar = var;
                    minCovar = c;
                }
            }

            double[][] tmpcovars = new double[covars.length][covars[0].length - 1];
            int cctr = 0;
            for (int c = 0; c < covars[0].length; c++) {
                if (c != minCovar) {
                    for (int r = 0; r < covars[0].length; r++) {
                        tmpcovars[r][cctr] = covars[r][c];
                    }
                    cctr++;
                }
            }
            covars = tmpcovars;
        }
        return covars;
    }

    private DoubleMatrixDataset<String, String> removeCovarWithLowestVariance(DoubleMatrixDataset<String, String> covars, int nrToRemove) {
        ArrayList<Pair<Integer, Double>> pairs = new ArrayList<>();
        for (int col = 0; col < covars.columns(); col++) {
            DoubleMatrix1D colobj = covars.getCol(col);
            double var = Descriptives.variance(colobj.toArray());
            pairs.add(new Pair<>(col, var, Pair.SORTBY.RIGHT));
        }

        // sort ascending
        Collections.sort(pairs);
        double[][] output = new double[covars.rows()][covars.columns() - nrToRemove];
        for (int r = 0; r < covars.rows(); r++) {
            int n = 0;
            // start iterating from nrToRemove
            for (int c = nrToRemove; c < pairs.size(); c++) {
                output[r][n] = covars.getElement(r, pairs.get(c).getLeft());
                n++;
            }
        }

        ArrayList<String> newCols = new ArrayList<>();
        for (int c = nrToRemove; c < pairs.size(); c++) {
            newCols.add(covars.getColObjects().get(pairs.get(c).getLeft()));
        }

        DoubleMatrixDataset<String, String> ds = new DoubleMatrixDataset<>();
        ds.setMatrix(output);
        try {
            ds.setRowObjects(covars.getRowObjects());
            ds.setColObjects(newCols);
        } catch (Exception e) {
            e.printStackTrace();
        }
        return ds;
    }

}
