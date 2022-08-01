package mbqtl;

import com.itextpdf.text.DocumentException;
import mbqtl.data.Dataset;
import mbqtl.gfx.QTLPanel;
import mbqtl.stat.PVal;
import mbqtl.stat.RankArray;
import mbqtl.vcf.VCFTabix;
import mbqtl.vcf.VCFVariant;
import umcg.genetica.containers.Triple;
import umcg.genetica.enums.Chromosome;
import umcg.genetica.features.Feature;
import umcg.genetica.graphics.Grid;
import umcg.genetica.graphics.panels.ScatterplotPanel;
import umcg.genetica.math.stats.Correlation;
import umcg.genetica.math.stats.ZScores;
import umcg.genetica.util.Primitives;

import java.io.IOException;
import java.text.DecimalFormat;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Iterator;
import java.util.stream.IntStream;

public class MbQTLPlot extends QTLAnalysis {

    private int nrPermutations = 1000;
    private int cisWindow = 1000000;
    private long randomSeed = 123456789;
    private boolean rankData = true;
    private boolean outputAll = false;
    private boolean replaceMissingGenotypes = false;

    public MbQTLPlot(String vcfFile, int chromosome, String linkfile, String snpLimitFile, String geneLimitFile, String snpGeneLimitFile, String geneExpressionDataFile, String geneAnnotationFile, String outputPrefix) throws IOException {
        super(vcfFile, chromosome, linkfile, snpLimitFile, geneLimitFile, snpGeneLimitFile, geneExpressionDataFile, geneAnnotationFile, outputPrefix);
    }

    public void setNrPermutations(int nrPermutations) {
        this.nrPermutations = nrPermutations;
    }

    public void setCisWindow(int cisWindow) {
        this.cisWindow = cisWindow;
    }

    public void setRandomSeed(long randomSeed) {
        this.randomSeed = randomSeed;
    }

    public void setRankData(boolean rankData) {
        this.rankData = rankData;
    }

    public void setOutputAll(boolean outputAll) {
        this.outputAll = outputAll;
    }

    public void setReplaceMissingGenotypes(boolean replaceMissingGenotypes) {
        this.replaceMissingGenotypes = replaceMissingGenotypes;
    }

    public void plot() throws IOException {
        System.out.println();
        System.out.println("-----------");
        System.out.println("QTL plotter");
        System.out.println("-----------");


        for (int g = 0; g < expressionData.genes.length; g++) {
            String gene = expressionData.genes[g];
            Integer geneAnnotationId = geneAnnotation.getGeneId(gene);


            if (geneAnnotationId != null) {
                // get chromosome
                int chr = geneAnnotation.getChr(geneAnnotationId);
                Chromosome chromosomeObj = Chromosome.parseChr("" + chr);
                System.out.println(g + "/" + expressionData.genes.length + " - Plotting " + gene + " on chr " + chr);
                String tmpVcf = origvcfFile;
                if (origvcfFile.contains("CHR")) {
                    tmpVcf = origvcfFile.replaceAll("CHR", "" + chr);
                }
                System.out.println("Opening " + tmpVcf);
                VCFTabix tabix = new VCFTabix(tmpVcf);


                double[] expData = expressionData.data[g];

                // define CIS window
                int pos = geneAnnotation.getPos(geneAnnotationId);
                String geneSymbol = geneAnnotation.getSymbol(geneAnnotationId);
                int start = pos - cisWindow;
                if (start < 0) {
                    start = 0;
                }
                int stop = pos + cisWindow;
                Feature cisRegion = new Feature(chromosomeObj, start, stop);
                System.out.println("Region: " + cisRegion.toString());

                // split expression data per dataset
                double[][] expressionPerDataset = new double[datasets.length][];
                IntStream.range(0, datasets.length).forEach(d -> {
                    Dataset thisDataset = datasets[d];
                    double[] datasetExp = thisDataset.select(expData, thisDataset.expressionIds);
                    double[] datasetExpRanked = datasetExp;
                    if (rankData) {
                        RankArray ranker = new RankArray();
                        datasetExpRanked = ranker.rank(datasetExp, true); // does this work with NaNs? answer: no
                    }
                    expressionPerDataset[d] = datasetExpRanked;

                });

                System.out.println("Querying tabix..");
                Iterator<VCFVariant> snpIterator = tabix.getVariants(cisRegion, genotypeSamplesToInclude, snpLimitSet);

                while (snpIterator.hasNext()) {
                    VCFVariant variant = snpIterator.next();

                    if (variant != null) {
                        String variantId = variant.getId();
                        if ((snpLimitSet == null || snpLimitSet.contains(variantId)) || (snpGeneLimitSet == null || (snpGeneLimitSet.containsKey(gene) && snpGeneLimitSet.get(gene).contains(variantId)))) {
                            final double[] genotypes = getGenotype(variant.getGenotypesAsByteVector());
                            final double[] dosages = getDosage(variant.getDosage());

                            // split genotype data per dataset, perform QC
                            double[][] genotypesPerDataset = new double[datasets.length][];
                            double[][] dosagesPerDataset = new double[datasets.length][];
                            VariantQCObj[] qcobjs = new VariantQCObj[datasets.length];
                            IntStream.range(0, datasets.length).forEach(d -> {
                                Dataset thisDataset = datasets[d];
                                dosagesPerDataset[d] = thisDataset.select(dosages, thisDataset.genotypeIds); // select required dosages
                                genotypesPerDataset[d] = thisDataset.select(genotypes, thisDataset.genotypeIds); // select required genotype IDs

                                VariantQCObj qcobj = checkVariant(genotypesPerDataset[d]);
                                if (qcobj.passqc) {
                                    if (replaceMissingGenotypes) {
                                        // only replace missing genotypes on variants that pass the qc thresholds
                                        double meanDosage = Util.meanGenotype(dosagesPerDataset[d]);
                                        double meanGenotype = Util.meanGenotype(genotypesPerDataset[d]);
                                        for (int i = 0; i < dosagesPerDataset[d].length; i++) {
                                            if (genotypesPerDataset[d][i] == -1) {
                                                genotypesPerDataset[d][i] = meanGenotype;
                                                dosagesPerDataset[d][i] = meanDosage;
                                            }
                                        }
                                    }

                                    // prune the data here once, to check if there are enough values to go ahead with this snp/gene combo
                                    // but only if the variant is passing the QC in the first place for the samples selected in this dataset
                                    Triple<double[], double[], double[]> prunedDatasetData = pruneMissingValues(genotypesPerDataset[d], dosagesPerDataset[d], expressionPerDataset[d]);

                                    // check the variant again, taking into account missingness in the expression data
                                    qcobj = checkVariant(prunedDatasetData.getLeft());

                                    // require minimum number of observations, otherwise kick out dataset from analysis
                                    if (prunedDatasetData.getLeft().length < minObservations) {
                                        qcobj.passqc = false;
                                    }
                                }
                                qcobjs[d] = qcobj;
                            });

                            int nrPassingQC = 0;
                            for (int d = 0; d < qcobjs.length; d++) {
                                if (qcobjs[d].passqc) {
                                    nrPassingQC++;
                                }
                            }
                            if (nrPassingQC >= minNumberOfDatasets) {
//                                Grid grid = new Grid(200, 200, 1, nrPassingQC, 100, 100);
                                Grid qtlPlotGrid = new Grid(300, 300, 1, datasets.length + 1, 100, 50);
                                QTLPanel metaAnalysisPanel = new QTLPanel(1, 1);
                                ArrayList<QTLPanel> panels = new ArrayList<>();

                                // iterate datasets
                                String geneStr = gene + " (" + geneSymbol + ")";
                                ArrayList<Double> metaValsX = new ArrayList<>();
                                ArrayList<Double> metaValsY = new ArrayList<>();

                                double[] zscores = new double[nrPassingQC];
                                int[] samplesizes = new int[nrPassingQC];

                                int passqcctr = 0;
                                for (int d = 0; d < datasets.length; d++) {
                                    Dataset thisDataset = datasets[d];
                                    double[] datasetGt = genotypesPerDataset[d]; // thisDataset.select(genotypes, thisDataset.genotypeIds); // select required genotype IDs
                                    VariantQCObj qcobj = qcobjs[d]; // check maf, hwep, call-rate, number of genotypes per genotype group
                                    if (!qcobj.passqc) {
                                        QTLPanel dsPanel = new QTLPanel(1, 1);
                                        dsPanel.setNotTested();
                                        dsPanel.setDatasetDetails(thisDataset.name, geneStr, variantId, 0, 1, 0);
                                        panels.add(dsPanel);
                                    } else {
                                        double[] datasetExp = expressionPerDataset[d];
                                        double[] datasetExpCopy = new double[datasetExp.length];
                                        System.arraycopy(datasetExp, 0, datasetExpCopy, 0, datasetExpCopy.length);

                                        double[] datasetDs = dosagesPerDataset[d];
//                                        System.out.println(gene);
//                                        for (int q = 0; q < datasetGt.length; q++) {
//                                            System.out.println(thisDataset.name + "\t" + datasetGt[q] + "\t" + datasetDs[q] + "\t" + datasetExp[q]);
//                                        }

                                        // prune the data (remove missing values)
                                        // can't prune the data earlier (would save a lot of compute time) because shuffling is performed over all available samples for this dataset
                                        // this is because the order of permuted samples should be equal across all SNPs
                                        Triple<double[], double[], double[]> prunedDatasetData = pruneMissingValues(datasetGt, datasetDs, datasetExpCopy);

                                        // re-rank data here? original EMP does not, but it is the right thing to do...
                                        double[] datasetExpPruned = prunedDatasetData.getRight();
//                                    if (rankData) {
//                                        RankArray ranker = new RankArray();
//                                        datasetExpPruned = ranker.rank(datasetExpPruned, true); // does this work with NaNs? answer: no
//                                    }
                                        double[] origExp = datasetExpPruned;
                                        datasetExpPruned = Util.centerScale(datasetExpPruned);
                                        double[] datasetDsPruned = Util.centerScale(prunedDatasetData.getMiddle());
                                        double[] datasetGtPruned = Util.centerScale(prunedDatasetData.getLeft());

                                        if (datasetExpPruned.length == 0) {
                                            System.out.println(datasets[d].name + ": exp: " + datasetExpCopy.length + "\tgt: " + datasetDs.length + "\t has no genotypes after pruning, but passqc: " + qcobj.passqc);
                                            System.out.println(datasets[d].name + ": cr: " + qcobj.cr + "\thwep: " + qcobj.hwep + "\tmaf: " + qcobj.maf);
                                            QTLPanel dsPanel = new QTLPanel(1, 1);
                                            dsPanel.setNotTested();
                                            dsPanel.setDatasetDetails(thisDataset.name, geneStr, variantId, 0, 1, 0);
                                            panels.add(dsPanel);
                                        } else {
                                            // perform correlation
                                            double r = Correlation.correlate(datasetDsPruned, datasetExpPruned);
                                            double p = PVal.getPvalue(r, datasetExpPruned.length - 2);
                                            double z = ZScores.pToZTwoTailed(p); // p value is already two-tailed, so need to use this other p-value conversion method... :/; returns negative z-scores by default
                                            if (r > 0) {
                                                z *= -1; // flip z-score if correlation is positive because p-value conversion returns only negative z-scores
                                            }

                                            zscores[passqcctr] = z;
                                            samplesizes[passqcctr] = datasetExpPruned.length;
                                            // add panel to grid

//                                        System.out.println(gene);
//                                        for (int q = 0; q < datasetDsPruned.length; q++) {
//                                            System.out.println(thisDataset.name + "\t" + datasetGtPruned[q] + "\t" + datasetDsPruned[q] + "\t" + datasetExpPruned[q]);
//                                        }
                                            QTLPanel dsPanel = new QTLPanel(1, 1);
                                            double[] xvals = Arrays.stream(prunedDatasetData.getMiddle()).toArray();
                                            double[] yvals = origExp;
                                            for (int v = 0; v < xvals.length; v++) {
                                                metaValsX.add(xvals[v]);
                                                metaValsY.add(yvals[v]);
                                            }
                                            dsPanel.setData(xvals, yvals);
                                            dsPanel.setAlleles(variant.getAlleles());
                                            dsPanel.setDatasetDetails(thisDataset.name, geneStr, variantId, z, p, r);
                                            panels.add(dsPanel);
                                            passqcctr++;
                                        }


                                    }
                                }

                                // define meta-analysis panel
                                metaAnalysisPanel.setData(Primitives.toPrimitiveArr(metaValsX), Primitives.toPrimitiveArr(metaValsY));
                                metaAnalysisPanel.setAlleles(variant.getAlleles());

                                double metaZ = ZScores.getWeightedZ(zscores, samplesizes);
                                double metaP = ZScores.zToP(metaZ);

                                metaAnalysisPanel.setDatasetDetails("Meta-analysis", geneStr, variantId, metaZ, metaP, Double.NaN);
                                qtlPlotGrid.addPanel(metaAnalysisPanel);
                                for (QTLPanel p : panels) {
                                    qtlPlotGrid.addPanel(p);
                                }

                                gene = gene.replaceAll("\\|", "_");
                                gene = gene.replaceAll("\\.", "-");
                                geneSymbol = geneSymbol.replaceAll("\\|", "_");
                                geneSymbol = geneSymbol.replaceAll("\\.", "-");
                                if (gene.length() > 20) {
                                    gene = gene.substring(0, 20);
                                }
                                String snp = variant.getId().replaceAll(":", "_");
                                try {
//                                    String fileout = outputPrefix + "-" + gene + "-" + snp + ".pdf";
//                                    System.out.println("Plotting: " + fileout);
//                                    grid.draw(fileout);
                                    String fileout = outputPrefix + "qtlplot-" + chr + "_" + geneAnnotation.getPos(geneAnnotationId) + "-" + gene + "-" + geneSymbol + "_" + snp + ".pdf";
                                    System.out.println("Plotting: " + fileout);
                                    qtlPlotGrid.draw(fileout);

                                } catch (DocumentException e) {
                                    e.printStackTrace();
                                }
                            }
                        }
                    }
                }
            }
        }

    }
}
