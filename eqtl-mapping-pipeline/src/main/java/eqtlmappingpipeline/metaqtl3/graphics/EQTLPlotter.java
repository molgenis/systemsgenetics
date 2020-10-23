/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package eqtlmappingpipeline.metaqtl3.graphics;

import cern.colt.matrix.tint.IntMatrix2D;
import umcg.genetica.math.stats.Descriptives;
import eqtlmappingpipeline.metaqtl3.containers.Settings;
import eqtlmappingpipeline.metaqtl3.containers.Result;
import eqtlmappingpipeline.metaqtl3.containers.WorkPackage;

import java.awt.Color;
import java.awt.Graphics2D;
import java.awt.RenderingHints;
import java.awt.image.BufferedImage;
import java.io.File;
import java.text.DecimalFormat;
import java.text.DecimalFormatSymbols;
import java.util.ArrayList;

import umcg.genetica.util.Primitives;
import umcg.genetica.graphics.ViolinBoxPlot;
import umcg.genetica.io.trityper.SNP;
import umcg.genetica.io.trityper.TriTyperExpressionData;
import umcg.genetica.io.trityper.TriTyperGeneticalGenomicsDataset;
import umcg.genetica.io.trityper.util.BaseAnnot;
import umcg.genetica.math.stats.Regression;
import umcg.genetica.util.RankArray;

/**
 * @author harm-jan
 */
public class EQTLPlotter {

    private static final int FILE_TYPE_PNG = 1;
    private static final int FILE_TYPE_PDF = 2;
    private int outputPlotsFileType = FILE_TYPE_PDF;
    private static final DecimalFormatSymbols dfs = new DecimalFormatSymbols(java.util.Locale.US);
    private static final DecimalFormat df1 = new DecimalFormat("###.#;-###.#", dfs);
    private static final DecimalFormat df2 = new DecimalFormat("###.######;-###.######", dfs);
    private static final DecimalFormat df3 = new DecimalFormat("##.00;-##.00", dfs);
    private static final DecimalFormat df5 = new DecimalFormat("0.##E0", dfs); //0.##E0
    private static final DecimalFormat df6 = new DecimalFormat("##.###;-##.###", dfs);
    private final TriTyperGeneticalGenomicsDataset[] m_gg;
    private final String m_outputDir;
    private final String[] m_probeList;
    private final IntMatrix2D m_probeTranslation;
    private final boolean m_cisOnly;
    private final boolean m_parametricAnalysis;
    private final RankArray m_rda;
    private static final Color red = new Color(255, 0, 0);
    private static final Color black = new Color(0, 0, 0);
    private static final Color darkgray = new Color(100, 100, 100);
    private static final Color white = new Color(255, 255, 255);
    private static final Color grey = new Color(125, 125, 125);
    private static final Color blue = new Color(0, 0, 255);
    private static final Color lightgray = new Color(200, 200, 200);
    private static final Color medgray = new Color(150, 150, 150);
    private final boolean m_giveTiesSameRank;

    public EQTLPlotter(TriTyperGeneticalGenomicsDataset[] gg, Settings settings, String[] probeList, IntMatrix2D probeTranslation) {
        this.m_gg = gg;
        this.m_outputDir = settings.plotOutputDirectory;
        m_probeList = probeList;
        if (settings.cisAnalysis && !settings.transAnalysis) {
            m_cisOnly = true;
        } else {
            m_cisOnly = false;
        }
        m_probeTranslation = probeTranslation;

        m_parametricAnalysis = settings.performParametricAnalysis;
        m_rda = new RankArray();
        m_giveTiesSameRank = settings.equalRankForTies;
    }

    /**
     * Draw individual eQTL plot for a given SNP-probe combination for all
     * loaded datasets and provide additional summary statistics. This method
     * uses data contained within GeneticalGenomicsDataset. If you want to plot
     * a particular SNP-Probe eQTL, first load the actual SNP, using
     * GeneticalGenomicsDataset.loadSNP, then invoke this method
     *
     * @param file                                       Output file (depending on 'outputPlotsFileType' this should
     *                                                   be a PNG or PDF file)
     * @param snpInformation                             SNP information (SNPName (chr. SNPChr, SNPChrPos))
     * @param probeInformation                           Probe information (ProbeName (chr. ProbeChr,
     *                                                   ProbeChrStartPos-ProbeChrEndPos), GeneName)
     * @param pValueOverall                              Overall P-Value
     * @param pValueOverallAbsolute                      Overall P-Value when using absolute Z-Scores
     *                                                   in the meta-analysis
     * @param gg                                         Genetical Genomics Datasets
     * @param datasetIncluded                            Which Genetical Genomics Datasets should we
     *                                                   include
     * @param takeNegativeCorrelationAsAllelesAreFlipped For which Genetical
     *                                                   Genomics Datasets should we take negative correlations, as the actual
     *                                                   alleles have been flipped
     * @param probeName                                  Unique probe identifier, used to get the expression data
     *                                                   for this eQTL
     */
    //public void drawPlot(File file, String snpInformation, String probeInformation, double pValueOverall, double pValueOverallAbsolute, GGDataset[] gg, SNP[] loadedSNP, boolean[] datasetIncluded, boolean[] takeNegativeCorrelationAsAllelesAreFlipped, String probeName, boolean performParametricAnalysis, boolean onlyPerformCiseQTLAnalysis) {
    public void draw(WorkPackage wp, int pid) {
        System.setProperty("java.awt.headless", "true");


        boolean jitter = true;
        //Init image:
        int width = 40 + m_gg.length * 200;
        if (width < 440) {
            width = 440;
        }
        int height = 340;
        int margin = 30;
        int x0 = margin;
        int x1 = width - margin;
        int innerWidth = x1 - x0;
        int y0 = margin + 40;
        int y1 = height - margin - 12 - 10;
        int innerHeight = y1 - y0;

        double metaPvalue = wp.results.pvalues[pid];

        SNP[] snps = wp.getSnps();

        String snpName = "";
        byte snpChr = 0;
        int snpChrpos = 0;

        int[] probes = wp.getProbes();
        int probeId = pid;

        if (probes != null && probes.length == wp.results.pvalues.length) {
            probeId = probes[pid];
        }

        String probeName = m_probeList[probeId];
        for (SNP snp : snps) {
            if (snp != null) {
                snpName = snp.getName();
                snpChr = snp.getChr();
                snpChrpos = snp.getChrPos();
                break;
            }
        }
        String probeAnnot = "";
        for (int d = 0; d < snps.length; d++) {
            if (m_probeTranslation.get(d, probeId) != -9) {
                Integer realProbeId = m_probeTranslation.get(d, probeId);
                TriTyperExpressionData expressionData = m_gg[d].getExpressionData();
                probeAnnot = expressionData.getAnnotation()[realProbeId] + ", Chr: " + expressionData.getChr()[realProbeId] + " (" + expressionData.getChrStart()[realProbeId] + " - " + expressionData.getChrStop()[realProbeId] + ")";
                break;
            }
        }

        double logPValue = -Math.log10(metaPvalue);
        String logPValueString = df3.format(logPValue);
        if (logPValueString.startsWith(".")) {
            logPValueString = "0" + logPValueString;
        }

        String snpNameFix = snpName.replaceAll("/", "_").replaceAll(":", "_");
        String probeNameFix = probeName.replace("/", "_").replace(";", "_").replace("|", "").replace("__", "_");
        String fileName = m_outputDir + "" + logPValueString + "-" + snpNameFix + "-" + probeNameFix + ".pdf";
//        fileName = m_outputDir + "" + logPValueString + "-" + snpName + "-" + probeName + ".pdf";
        File file = new File(fileName);

        Graphics2D g2d = null;
        BufferedImage bi = null;
        com.itextpdf.text.Document document = null;
        com.itextpdf.text.pdf.PdfContentByte cb = null;
        com.itextpdf.text.pdf.PdfWriter writer = null;
        if (outputPlotsFileType == FILE_TYPE_PNG) {
            bi = new java.awt.image.BufferedImage(width, height, java.awt.image.BufferedImage.TYPE_INT_RGB);
            g2d = bi.createGraphics();
        } else {
            com.itextpdf.text.Rectangle rectangle = new com.itextpdf.text.Rectangle(width, height);
            document = new com.itextpdf.text.Document(rectangle);

            try {
                writer = com.itextpdf.text.pdf.PdfWriter.getInstance(document, new java.io.FileOutputStream(file));
                document.open();
                cb = writer.getDirectContent();
                cb.saveState();
                g2d = cb.createGraphics(width, height);
            } catch (Exception e) {
                System.out.println("Cannot write to PDF file!:\t" + file.getAbsolutePath());
                System.exit(-1);
            }

        }

        g2d.setRenderingHint(RenderingHints.KEY_ANTIALIASING, RenderingHints.VALUE_ANTIALIAS_ON);
        g2d.setColor(white);
        g2d.fillRect(0, 0, width, height);

        //Draw basic eQTL information:
        g2d.setColor(black);
        g2d.setComposite(java.awt.AlphaComposite.getInstance(java.awt.AlphaComposite.SRC_OVER, 1.0f));
        g2d.setFont(new java.awt.Font("Arial", java.awt.Font.BOLD, 9));
        g2d.drawString("SNP", margin, 15);
        g2d.drawString("Probe", margin, 25);
        g2d.drawString("P-Value", margin, 35);

        g2d.setFont(new java.awt.Font("Arial", java.awt.Font.PLAIN, 9));
        g2d.drawString(snpName + " - Chr: " + snpChr + " (" + snpChrpos + ")", margin + 40, 15);
        g2d.drawString(probeName + " - " + probeAnnot, margin + 40, 25);


        String pValueOverallString = df5.format(metaPvalue);
        if (metaPvalue > 0.001) {
            pValueOverallString = df6.format(metaPvalue);
        }

        g2d.drawString(pValueOverallString, margin + 40, 35);

        int numDatasets = m_gg.length;

        Boolean[] flipalleles = wp.getFlipSNPAlleles();
        Result results = wp.results;
        //Draw individual plots:
        for (int d = 0; d < numDatasets; d++) {
            //Draw dataset:
            g2d.setColor(black);
            g2d.setComposite(java.awt.AlphaComposite.getInstance(java.awt.AlphaComposite.SRC_OVER, 1.0f));
            g2d.setFont(new java.awt.Font("Arial", java.awt.Font.BOLD, 10));
            g2d.drawString(m_gg[d].getSettings().name, margin + d * 200, 53);

            TriTyperGeneticalGenomicsDataset currentDataset = m_gg[d];
            SNP currentSNP = snps[d];
            Integer probe = currentDataset.getExpressionData().getProbeToId().get(probeName);
            if (currentSNP != null && probe != -9) {

                //Define x-axis data:

                int itr = 0;
                int[] nrSamplesPerX = new int[3];
                int nrSamplesWithData = results.numSamples[d];
                double[] x = new double[nrSamplesWithData];
                int[] xBinary = new int[nrSamplesWithData];
                Boolean[] indIsFemale = new Boolean[x.length];
                Boolean[] isFemale = currentDataset.getGenotypeData().getIsFemale();
                int numSamples = currentDataset.getTotalGGSamples();
                int[] indWGA = currentDataset.getExpressionToGenotypeIdArray();
                byte[] genotypes = currentSNP.getGenotypes();
                double[] rawData = currentDataset.getExpressionData().getMatrix()[probe];
                double[] y = new double[nrSamplesWithData];


                ArrayList<Double> valsAA = new ArrayList<Double>();
                ArrayList<Double> valsAB = new ArrayList<Double>();
                ArrayList<Double> valsBB = new ArrayList<Double>();

                int minYAA = Integer.MAX_VALUE;
                int maxYAA = -Integer.MAX_VALUE;
                int minYAB = Integer.MAX_VALUE;
                int maxYAB = -Integer.MAX_VALUE;
                int minYBB = Integer.MAX_VALUE;
                int maxYBB = -Integer.MAX_VALUE;


                for (int s = 0; s < numSamples; s++) {
                    int ind = indWGA[s];
                    if (ind != -1) {
                        double valX = genotypes[ind];
                        if (valX != -1) {
                            if (flipalleles[d]) {
                                nrSamplesPerX[(int) (2 - valX)]++;
                            } else {
                                nrSamplesPerX[(int) valX]++;
                            }
                            xBinary[itr] = (int) valX;
                            if (wp.getFlipSNPAlleles()[d]) {
                                xBinary[itr] = 2 - (int) valX;
                            }

                            if (currentSNP.hasDosageInformation()) {
                                valX = currentSNP.getDosageValues()[ind];
                            }
                            if (flipalleles[d]) {
                                valX = 2 - valX;
                            }
                            x[itr] = valX;
                            indIsFemale[itr] = isFemale[ind];
                            itr++;
                        }
                    }
                }

                //Define y-axis data:
                if (nrSamplesWithData == numSamples) {
                    //All genotypes have been succesfully called, use quick approach:
                    y = rawData;
                } else {
                    //Not all genotypes have been succesfully called, use slow approach:
                    itr = 0;
                    for (int s = 0; s < numSamples; s++) {
                        int ind = indWGA[s];
                        if (ind != -1) {
                            int valX = currentSNP.getGenotypes()[ind];
                            if (valX != -1) {
                                y[itr] = rawData[s];
                                itr++;
                            }
                        }
                    }
                }


                //Draw graph:
                g2d.setColor(black);

                //Get minimal and maximal expression for this probe:
                double minExpression = Double.MAX_VALUE;
                double maxExpression = Double.MIN_VALUE;
                for (int i = 0; i < y.length; i++) {
                    if (y[i] < minExpression) {
                        minExpression = y[i];
                    }
                    if (y[i] > maxExpression) {
                        maxExpression = y[i];
                    }
                }

                //Draw regression line:
                double[] correlationValues = Regression.getLinearRegressionCoefficients(x, y);
                g2d.setComposite(java.awt.AlphaComposite.getInstance(java.awt.AlphaComposite.SRC_OVER, 0.25f));
                g2d.setStroke(new java.awt.BasicStroke(2));
                int pixelY1 = y1 - (int) ((correlationValues[1] - minExpression) / (maxExpression - minExpression) * (double) innerHeight);
                int pixelY2 = y1 - (int) (((2 * correlationValues[0] + correlationValues[1]) - minExpression) / (maxExpression - minExpression) * (double) innerHeight);
                int pixelXGroup0 = x0 + d * 200 + (int) (0 * 50.0d);
                int pixelXGroup2 = x0 + d * 200 + (int) (2 * 50.0d);
                g2d.drawLine(pixelXGroup0 + 6, pixelY1, pixelXGroup2 + 6, pixelY2);
                g2d.setStroke(new java.awt.BasicStroke(1));

                //Draw individual measurements:
                g2d.setComposite(java.awt.AlphaComposite.getInstance(java.awt.AlphaComposite.SRC_OVER, 0.50f));
                boolean[][][] positionToPlot = new boolean[3][41][innerHeight + 5];
                g2d.setFont(new java.awt.Font("Arial", java.awt.Font.PLAIN, 6));
                for (int i = 0; i < x.length; i++) {

                    int pixelX = x0 + d * 200 + (int) (x[i] * 50.0d);
                    if (jitter) {
                        pixelX = x0 + d * 200 + (int) (xBinary[i] * 50.0d);
                    }
                    int posY = (int) ((y[i] - minExpression) / (maxExpression - minExpression) * (double) innerHeight);
                    int pixelY = y1 - posY;

                    switch (xBinary[i]) {
                        case 0:
                            valsAA.add(y[i]);
                            if (pixelY > maxYAA) {
                                maxYAA = pixelY;
                            }
                            if (pixelY < minYAA) {
                                minYAA = pixelY;
                            }
                            break;
                        case 1:
                            valsAB.add(y[i]);
                            if (pixelY > maxYAB) {
                                maxYAB = pixelY;
                            }
                            if (pixelY < minYAB) {
                                minYAB = pixelY;
                            }
                            break;
                        case 2:
                            valsBB.add(y[i]);
                            if (pixelY > maxYBB) {
                                maxYBB = pixelY;
                            }
                            if (pixelY < minYBB) {
                                minYBB = pixelY;
                            }
                            break;
                    }


                    if (!currentSNP.hasDosageInformation() || jitter) {
                        //Draw expression dots in non-overlapping manner:
                        int offsetToPlot = 40;
                        for (int offset = 0; offset < 41; offset++) {
                            boolean allClear = true;
                            for (int q = 0; q < 4; q++) {
                                if (positionToPlot[xBinary[i]][offset][posY + q]) {
                                    allClear = false;
                                    break;
                                }
                            }
                            if (allClear) {
                                for (int q = 0; q < 4; q++) {
                                    positionToPlot[xBinary[i]][offset][posY + q] = true;
                                }
                                offsetToPlot = offset;
                                break;
                            }
                        }
                        double diffOffset = 0;
                        if (offsetToPlot > 0) {
                            double sign = offsetToPlot % 2;
                            if (sign == 0) {
                                sign = -1;
                            }
                            diffOffset = 1 * sign * (offsetToPlot - (offsetToPlot + 1) % 2 + 1);
                        }
                        pixelX += diffOffset;
                    }

                    if (indIsFemale[i] == null) {
                        g2d.setColor(grey);
                    } else if (indIsFemale[i]) {
                        g2d.setColor(red);
                    } else {
                        g2d.setColor(blue);
                    }

                    g2d.fillOval(pixelX + 6, pixelY, 2, 2);
                }

                // draw Spearman's rank correlation coefficient:
                double correlation = results.correlations[d][pid];
                g2d.setColor(black);
                g2d.setFont(new java.awt.Font("Arial", java.awt.Font.PLAIN, 10));
                g2d.setComposite(java.awt.AlphaComposite.getInstance(java.awt.AlphaComposite.SRC_OVER, 1.0f));
                String correlationString = df6.format(correlation);
                String r2String = df6.format(correlation * correlation);
                String allele0 = "";
                String allele1 = "";

                byte[] loadedSNPAlleles = currentSNP.getAlleles();
                allele0 = BaseAnnot.toString(loadedSNPAlleles[0]);
                allele1 = BaseAnnot.toString(loadedSNPAlleles[1]);

                if (wp.getFlipSNPAlleles()[d]) {
                    String alleleTemp = allele0;
                    allele0 = allele1;
                    allele1 = alleleTemp;
                }

                // overlay the boxplots
                ViolinBoxPlot boxplotter = new ViolinBoxPlot();
                java.awt.AlphaComposite alphaComposite100 = java.awt.AlphaComposite.getInstance(java.awt.AlphaComposite.SRC, 1.00f);
                g2d.setComposite(alphaComposite100);
                g2d.setColor(new Color(0, 0, 0));

                // g2d, int x, int y, int width, int height, double[] vals, double minValue, double maxValue, boolean drawOutliers

                double[] aaArr = toArray(valsAA);
                if (aaArr.length > 0) {
                    boxplotter.drawBoxPlot(g2d, margin + d * 200 + 0 - 5, minYAA, 25, (maxYAA - minYAA), aaArr, Primitives.min(aaArr), Primitives.max(aaArr), false);
                }

                double[] abArr = toArray(valsAB);

                if (abArr.length > 0) {
                    boxplotter.drawBoxPlot(g2d, margin + d * 200 + 50 - 5, minYAB, 25, (maxYAB - minYAB), abArr, Primitives.min(abArr), Primitives.max(abArr), false);
                }

                double[] bbArr = toArray(valsBB);

                if (bbArr.length > 0) {

                    boxplotter.drawBoxPlot(g2d, margin + d * 200 + 100 - 5, minYBB, 25, (maxYBB - minYBB), bbArr, Primitives.min(bbArr), Primitives.max(bbArr), false);
                }

                //Draw alleles:
                g2d.setColor(black);
                g2d.drawString(allele0 + allele0, margin + d * 200 + 0, height - 15 - 10);
                g2d.drawString(allele0 + allele1, margin + d * 200 + 50, height - 15 - 10);
                g2d.drawString(allele1 + allele1, margin + d * 200 + 100, height - 15 - 10);

                //Draw nr samples per genotype group:
                g2d.setColor(black);
                g2d.drawString("(" + nrSamplesPerX[0] + ")", margin + d * 200 + 0 + 17, height - 15 - 10);
                g2d.drawString("(" + nrSamplesPerX[1] + ")", margin + d * 200 + 50 + 17, height - 15 - 10);
                g2d.drawString("(" + nrSamplesPerX[2] + ")", margin + d * 200 + 100 + 17, height - 15 - 10);

                //Draw correlation:
                g2d.setColor(black);
                g2d.setFont(new java.awt.Font("Arial", java.awt.Font.BOLD, 9));
                g2d.drawString("Corr:", margin + d * 200, height - 2 - 10);
                g2d.drawString("R2:", margin + d * 200 + 75, height - 2 - 10);
                g2d.setFont(new java.awt.Font("Arial", java.awt.Font.PLAIN, 9));
                g2d.drawString(correlationString, margin + d * 200 + 37, height - 2 - 10);
                g2d.drawString(r2String, margin + d * 200 + 75 + 37, height - 2 - 10);

                //Draw Z-Score and P-Value:
                double zScore = results.zscores[d][pid];
                if (wp.getFlipSNPAlleles()[d]) {
                    zScore *= -1;
                }

                String zScoreString = df5.format(zScore);
                if (zScore < 100) {
                    zScoreString = df6.format(zScore);
                }

                double pValue = Descriptives.convertZscoreToPvalue(zScore);

                String pValueString = df5.format(pValue);
                if (pValue > 0.001) {
                    pValueString = df6.format(pValue);
                }

                g2d.setColor(black);
                g2d.setFont(new java.awt.Font("Arial", java.awt.Font.BOLD, 9));
                g2d.drawString("Z-Score:", margin + d * 200, height - 2);
                g2d.drawString("P-Value:", margin + d * 200 + 75, height - 2);
                g2d.setFont(new java.awt.Font("Arial", java.awt.Font.PLAIN, 9));
                g2d.drawString(zScoreString, margin + d * 200 + 37, height - 2);
                g2d.drawString(pValueString, margin + d * 200 + 75 + 37, height - 2);

                // System.out.println(currentDataset.getName() + "\t" + snpInformation + "\t" + zScore + "\t" + pValue + "\t" + correlation + "\t" + nrSamplesWithData);

                //Draw average expression rank for this probe, to permit comparison with other datasets, where eQTL effect is different:
                int probeCount = currentDataset.getExpressionData().getProbes().length;
                double[] means = new double[probeCount];
                for (int p = 0; p < probeCount; p++) {
                    means[p] = currentDataset.getExpressionData().getOriginalProbeMean()[probe];
                }

                double probeOriginalMeanRank = m_rda.rank(means, m_giveTiesSameRank)[probe];
                for (int p = 0; p < probeCount; p++) {
                    means[p] = currentDataset.getExpressionData().getProbeMean()[probe];
                }
                double probeMeanRank = m_rda.rank(means, m_giveTiesSameRank)[probe];
                g2d.setColor(lightgray); // 
                g2d.drawLine(margin + d * 200 + 140, y0, margin + d * 200 + 140, y1);
                g2d.drawLine(margin + d * 200 + 138, y0, margin + d * 200 + 142, y0);
                g2d.drawLine(margin + d * 200 + 138, y1, margin + d * 200 + 142, y1);
                g2d.setColor(medgray);
                int pixelYOriginal = y1 - (int) (probeOriginalMeanRank / (double) (probeCount) * (double) innerHeight);
                g2d.drawLine(margin + d * 200 + 135, pixelYOriginal, margin + d * 200 + 145, pixelYOriginal);
                g2d.setColor(darkgray);
                int pixelY = y1 - (int) (probeMeanRank / (double) (probeCount) * (double) innerHeight);
                g2d.drawLine(margin + d * 200 + 135, pixelY, margin + d * 200 + 145, pixelY);
                //System.out.println(d + "\t" + pixelYOriginal + "\t" + pixelY + "\t" + probeOriginalMeanRank + "\t" + probeMeanRank + "\t" + gg[d].probeOriginalMean[probe] + "\t" + gg[d].probeMean[probe] + "\t" + gg[d].probeOriginalVariance[probe] + "\t" + gg[d].probeVariance[probe]);


            } else { // if !datasetincluded[d]

                //eQTL data is not available for this datasets, provide information why not:
                g2d.setColor(black);
                g2d.setComposite(java.awt.AlphaComposite.getInstance(java.awt.AlphaComposite.SRC_OVER, 1.0f));

                g2d.setFont(new java.awt.Font("Arial", java.awt.Font.BOLD, 9));
                g2d.drawString("eQTL not available", margin + d * 200, 70);

                g2d.setFont(new java.awt.Font("Arial", java.awt.Font.PLAIN, 9));
                currentSNP = snps[d];
                if (currentSNP == null) {
                    //SNP has not been genotyped
                    g2d.drawString("SNP has not been genotyped.", margin + d * 200, 90);
                } else {
                    //SNP does not pass QC:
                    if (!currentSNP.passesQC()) {
                        g2d.drawString("SNP does not pass QC:", margin + d * 200, 90);
                        String callRateString = df1.format(currentSNP.getCR() * 100);
                        g2d.drawString("- Call rate: " + callRateString + "%", margin + d * 200, 100);
                        String mafString = df1.format(currentSNP.getMAF() * 100);
                        g2d.drawString("- MAF: " + mafString + "%", margin + d * 200, 110);
                        String hweString = df2.format(currentSNP.getHWEP());
                        g2d.drawString("- HWE P-Value: " + hweString, margin + d * 200, 120);
                    } else {
                        if (currentDataset.getExpressionData().getProbeToId().get(probeName) != -9) {
                            if (m_cisOnly) {
                                g2d.setColor(red);
                                g2d.drawString("Cis-probe maps too far away.", margin + d * 200, 90);
                                g2d.drawString("You probably have different probe", margin + d * 200, 100);
                                g2d.drawString("mappings in the various datasets!", margin + d * 200, 110);
                                g2d.setColor(red);
                            } else {
                                if (currentDataset.getExpressionData().getProbeToId().get(probeName) != -9) {
                                    g2d.setColor(red);
                                    g2d.drawString("Unknown why not included!", margin + d * 200, 90);
                                    g2d.drawString("This suggests a bug.", margin + d * 200, 100);
                                    g2d.setColor(black);
                                }
                            }
                        }
                    }
                    //Expression probe is not present in this dataset:
                    if (currentDataset.getExpressionData().getProbeToId().get(probeName) == -9) {
                        if (currentSNP.passesQC()) {
                            g2d.drawString("Probe not present.", margin + d * 200, 90);
                        } else {
                            g2d.drawString("Probe not present.", margin + d * 200, 140);
                        }
                    }

                }


            }
        } // end for each dataset

        //Save image:
        if (outputPlotsFileType == FILE_TYPE_PNG) {
            try {
                javax.imageio.ImageIO.write(bi, "png", file);
            } catch (Exception e) {
                System.out.println(e.getMessage());
                System.out.println(e.getStackTrace());
            }
        } else {

            g2d.dispose();
            cb.restoreState();
            document.close();
            writer.close();
        }


    }

    private double[] toArray(ArrayList<Double> vals) {
        double[] rrrr = new double[vals.size()];
        for (int i = 0; i < rrrr.length; i++) {
            rrrr[i] = vals.get(i);
        }
        return rrrr;
    }
}
