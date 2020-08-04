/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package nl.systemsgenetics.depict2.io;

import cern.colt.matrix.tdouble.DoubleMatrix1D;

import java.io.File;
import java.io.FileNotFoundException;
import java.io.FileOutputStream;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Collection;
import java.util.List;
import java.util.Map;

import hep.aida.tdouble.ref.DoubleVariableAxis;
import nl.systemsgenetics.depict2.Depict2;
import nl.systemsgenetics.depict2.Depict2Options;
import nl.systemsgenetics.depict2.Depict2Step2Results;
import nl.systemsgenetics.depict2.Depict2Step3Results;
import nl.systemsgenetics.depict2.gene.Gene;
import nl.systemsgenetics.depict2.pathway.PathwayAnnotations;
import nl.systemsgenetics.depict2.pathway.PathwayDatabase;
import nl.systemsgenetics.depict2.pathway.PathwayEnrichments;
import org.apache.poi.common.usermodel.HyperlinkType;
import org.apache.poi.ss.SpreadsheetVersion;
import org.apache.poi.ss.usermodel.CellStyle;
import org.apache.poi.ss.usermodel.CellType;
import org.apache.poi.ss.usermodel.CreationHelper;
import org.apache.poi.ss.usermodel.DataFormat;
import org.apache.poi.ss.usermodel.Font;
import org.apache.poi.ss.usermodel.Hyperlink;
import org.apache.poi.ss.usermodel.Workbook;
import org.apache.poi.ss.util.AreaReference;
import org.apache.poi.ss.util.CellReference;
import org.apache.poi.xssf.usermodel.XSSFCell;
import org.apache.poi.xssf.usermodel.XSSFRow;
import org.apache.poi.xssf.usermodel.XSSFSheet;
import org.apache.poi.xssf.usermodel.XSSFTable;
import org.apache.poi.xssf.usermodel.XSSFWorkbook;
import sun.nio.ch.IOUtil;
import umcg.genetica.math.matrix2.DoubleMatrix1dOrder;
import umcg.genetica.math.matrix2.DoubleMatrixDataset;
import umcg.genetica.math.stats.ZScores;

/**
 * @author patri
 */
public class ExcelWriter {

    private CellStyle zscoreStyle;
    private CellStyle largePvalueStyle;
    private CellStyle smallPvalueStyle;
    private CellStyle hlinkStyle;
    private CellStyle boldStyle;

    private String outputBasePath;
    private boolean hlaExcluded;
    private List<String> traits;
    private Depict2Options options;


    public ExcelWriter(List<String> traits, Depict2Options options) {
        this.outputBasePath = options.getOutputBasePath();
        this.hlaExcluded = options.isExcludeHla();
        this.traits = traits;
        this.options = options;
    }

    /**
     * @param pathwayEnrichments
     * @param pathwayDatabases
     * @param outputBasePath
     * @param traits
     * @param hlaExcluded
     * @throws java.io.FileNotFoundException
     */
    public void saveStep2Excel(Depict2Step2Results results) throws FileNotFoundException, IOException {

        final DoubleMatrixDataset<String, String> genePvalues = results.getGenePvalues();
        final List<PathwayEnrichments> pathwayEnrichments = results.getPathwayEnrichments();
        System.setProperty(" java.awt.headless", "true");

        // Each trait gets its own sheet
        for (String trait : traits) {

            Workbook enrichmentWorkbook = new XSSFWorkbook();
            setStyles(enrichmentWorkbook);
            CreationHelper createHelper = enrichmentWorkbook.getCreationHelper();

            // Overview sheet
            populateOverviewSheet(enrichmentWorkbook, trait, pathwayEnrichments, createHelper);

            // Sheet for gene pvalues
            populateGenePvalueSheet(enrichmentWorkbook, trait, genePvalues);

            // Sheet for each pathway database
            for (PathwayEnrichments pathwayEnrichment : pathwayEnrichments) {
                populatePathwaySheet(enrichmentWorkbook, pathwayEnrichment, trait, createHelper);
            }

            // Save the file
            File excelFile = new File(outputBasePath + "_enrichtments" + (traits.size() > 1 ? "_" + trait : "") + (hlaExcluded ? "_exHla.xlsx" : ".xlsx"));
            int nr = 1;
            while (excelFile.exists()) {
                excelFile = new File(outputBasePath + "_enrichtments" + (traits.size() > 1 ? "_" + trait : "") + (hlaExcluded ? "_exHla" : "") + "_" + nr + ".xlsx");
                nr++;
            }
            enrichmentWorkbook.write(new FileOutputStream(excelFile));
        }
    }


    public void saveStep3Excel(Depict2Step2Results step2, Depict2Step3Results step3) throws IOException {

    }

    private void setStyles(Workbook enrichmentWorkbook) {
        DataFormat format = enrichmentWorkbook.createDataFormat();

        zscoreStyle = enrichmentWorkbook.createCellStyle();
        zscoreStyle.setDataFormat(format.getFormat("0.00"));

        largePvalueStyle = enrichmentWorkbook.createCellStyle();
        largePvalueStyle.setDataFormat(format.getFormat("0.0000"));

        smallPvalueStyle = enrichmentWorkbook.createCellStyle();
        smallPvalueStyle.setDataFormat(format.getFormat("0.00E+0"));

        hlinkStyle = enrichmentWorkbook.createCellStyle();
        Font hlinkFont = enrichmentWorkbook.createFont();
        hlinkFont.setUnderline(Font.U_SINGLE);
        hlinkStyle.setFont(hlinkFont);

        boldStyle = enrichmentWorkbook.createCellStyle();
        Font fontBold = enrichmentWorkbook.createFont();
        fontBold.setBold(true);
        boldStyle.setFont(fontBold);
    }

    private void populateOverviewSheet(Workbook enrichmentWorkbook, String trait, Collection<PathwayEnrichments> pathwayEnrichments, CreationHelper createHelper) {
        // -----------------------------------------------------------------------
        // Create overview sheet
        // -----------------------------------------------------------------------
        XSSFSheet overviewSheet = (XSSFSheet) enrichmentWorkbook.createSheet("Overview");

        int r = 0;
        XSSFRow row = overviewSheet.createRow(r++);
        XSSFCell cell = row.createCell(0, CellType.STRING);
        cell.setCellValue("Pathway enrichment analysis for: " + trait);
        cell.setCellStyle(boldStyle);

        row = overviewSheet.createRow(r++);
        cell = row.createCell(0, CellType.STRING);
        cell.setCellValue("Generated using Downstreamer " + Depict2.VERSION);
        cell.setCellStyle(boldStyle);

        overviewSheet.createRow(r++);

        row = overviewSheet.createRow(r++);
        cell = row.createCell(0, CellType.STRING);
        cell.setCellValue("Gene set database");
        cell.setCellStyle(boldStyle);
        cell = row.createCell(1, CellType.STRING);
        cell.setCellValue("Number of sets");
        cell.setCellStyle(boldStyle);

        for (PathwayEnrichments pathwayEnrichment : pathwayEnrichments) {
            row = overviewSheet.createRow(r++);
            cell = row.createCell(0, CellType.STRING);
            cell.setCellValue(pathwayEnrichment.getPathwayDatabase().getName());

            Hyperlink link = createHelper.createHyperlink(HyperlinkType.DOCUMENT);
            link.setAddress(pathwayEnrichment.getPathwayDatabase().getName() + "!A1");
            cell.setHyperlink(link);
            cell.setCellStyle(hlinkStyle);

            row.createCell(1, CellType.NUMERIC).setCellValue(pathwayEnrichment.getNumberOfPathways());
        }

        for (int c = 0; c < 2; ++c) {
            overviewSheet.autoSizeColumn(c);
            overviewSheet.setColumnWidth(c, overviewSheet.getColumnWidth(c) + 500);//compensate for with auto filter and inaccuracies
        }

        overviewSheet.createRow(r++);

        row = overviewSheet.createRow(r++);
        cell = row.createCell(0, CellType.STRING);
        cell.setCellValue("Used settings");
        cell.setCellStyle(boldStyle);

        row = overviewSheet.createRow(r++);
        cell = row.createCell(0, CellType.STRING);
        cell.setCellValue("Number of permutations used for p-values: " + options.getPermutationPathwayEnrichment());

        row = overviewSheet.createRow(r++);
        cell = row.createCell(0, CellType.STRING);
        cell.setCellValue("Number of permutations used for FDR: " + options.getPermutationFDR());

        row = overviewSheet.createRow(r++);
        cell = row.createCell(0, CellType.STRING);
        cell.setCellValue("Gene pruning r: " + options.getGenePruningR());

        row = overviewSheet.createRow(r++);
        cell = row.createCell(0, CellType.STRING);
        cell.setCellValue("Force normal pathway scores: " + options.isForceNormalPathwayPvalues());

        row = overviewSheet.createRow(r++);
        cell = row.createCell(0, CellType.STRING);
        cell.setCellValue("Force normal GWAS gene z-scores: " + options.isForceNormalGenePvalues());

        row = overviewSheet.createRow(r++);
        cell = row.createCell(0, CellType.STRING);
        cell.setCellValue("Regress out gene lengths from GWAS gene z-scores: " + options.isRegressGeneLengths());

        if (options.isIgnoreGeneCorrelations()) {
            row = overviewSheet.createRow(r++);
            cell = row.createCell(0, CellType.STRING);
            cell.setCellValue("Ignoring gene correlations: " + options.isRegressGeneLengths());
        }
    }

    private void populatePathwaySheet(Workbook enrichmentWorkbook, PathwayEnrichments pathwayEnrichment, String trait, CreationHelper createHelper) throws IOException {
        PathwayDatabase pathwayDatabase = pathwayEnrichment.getPathwayDatabase();

        final PathwayAnnotations pathwayAnnotations = new PathwayAnnotations(new File(pathwayDatabase.getLocation() + ".colAnnotations.txt"));
        final int maxAnnotations = pathwayAnnotations.getMaxNumberOfAnnotations();
        final DoubleMatrixDataset<String, String> databaseEnrichmentZscores = pathwayEnrichment.getEnrichmentZscores();
        final DoubleMatrixDataset<String, String> databaseEnrichmentQvalues = pathwayEnrichment.getqValues();
        final ArrayList<String> geneSets = databaseEnrichmentZscores.getRowObjects();
        //final int currentTraitCol = databaseEnrichment.getColIndex(trait);
        final double bonferroniCutoff = 0.05 / pathwayEnrichment.getNumberOfPathways();

        final DoubleMatrix1D traitEnrichment = databaseEnrichmentZscores.getCol(trait);
        final DoubleMatrix1D traitQvalue = databaseEnrichmentQvalues.getCol(trait);
        int[] order = DoubleMatrix1dOrder.sortIndexReverse(traitEnrichment);

        boolean annotateWithGwasData = pathwayDatabase.getName().startsWith("AnnotateGwas");
        int gwasAnnotations = 0;
        if (annotateWithGwasData) {
            gwasAnnotations = 2;
        }

        XSSFSheet databaseSheet = (XSSFSheet) enrichmentWorkbook.createSheet(pathwayDatabase.getName());
        XSSFTable table = databaseSheet.createTable(new AreaReference(new CellReference(0, 0),
                new CellReference(databaseEnrichmentZscores.rows(), 5 + maxAnnotations + gwasAnnotations),
                SpreadsheetVersion.EXCEL2007));

        table.setName(pathwayDatabase.getName() + "_res");
        table.setDisplayName(pathwayDatabase.getName());
        table.setStyleName("TableStyleLight9");
        table.getCTTable().getTableStyleInfo().setShowRowStripes(true);
        table.getCTTable().addNewAutoFilter();

        //databaseSheet.createFreezePane(0, 1);
        XSSFRow headerRow = databaseSheet.createRow(0);
        int hc = 0;
        headerRow.createCell(hc++, CellType.STRING).setCellValue("Gene set");
        for (int i = 0; i < maxAnnotations; ++i) {
            headerRow.createCell(hc++, CellType.STRING).setCellValue("Annotation" + (i + 1));
        }
        headerRow.createCell(hc++, CellType.STRING).setCellValue("Enrichment Z-score");
        headerRow.createCell(hc++, CellType.STRING).setCellValue("Enrichment P-value");
        headerRow.createCell(hc++, CellType.STRING).setCellValue("Enrichment Q-value");
        headerRow.createCell(hc++, CellType.STRING).setCellValue("Bonferroni significant");
        headerRow.createCell(hc++, CellType.STRING).setCellValue("FDR 5% significant");

        // Loop over rows
        for (int r = 0; r < databaseEnrichmentZscores.rows(); ++r) {
            XSSFRow row = databaseSheet.createRow(r + 1); //+1 for header
            String geneSet = geneSets.get(order[r]);
            row.createCell(0, CellType.STRING).setCellValue(geneSet);

            if (maxAnnotations > 0) {
                ArrayList<String> thisPathwayAnnotations = pathwayAnnotations.getAnnotationsForPathway(geneSet);
                if (thisPathwayAnnotations == null) {
                    for (int j = 0; j < maxAnnotations; ++j) {
                        row.createCell(j + 1, CellType.STRING).setCellValue("");
                    }
                } else {
                    for (int j = 0; j < maxAnnotations; ++j) {
                        if (j < thisPathwayAnnotations.size()) {
                            String annotation = thisPathwayAnnotations.get(j);
                            XSSFCell cell = row.createCell(j + 1, CellType.STRING);
                            cell.setCellValue(annotation);

                            if (annotation.startsWith("http")) {
                                Hyperlink link = createHelper.createHyperlink(HyperlinkType.URL);
                                link.setAddress(annotation);
                                cell.setHyperlink(link);
                                cell.setCellStyle(hlinkStyle);
                            }

                        } else {
                            row.createCell(j + 1, CellType.STRING).setCellValue("");
                        }

                    }
                }
            }

            double zscore = traitEnrichment.getQuick(order[r]);
            double qvalue = traitQvalue.getQuick(order[r]);

            XSSFCell zscoreCell = row.createCell(1 + maxAnnotations, CellType.NUMERIC);
            zscoreCell.setCellValue(zscore);
            zscoreCell.setCellStyle(zscoreStyle);

            double pvalue = ZScores.zToP(zscore);

            XSSFCell pvalueCell = row.createCell(2 + maxAnnotations, CellType.NUMERIC);
            pvalueCell.setCellValue(pvalue);
            pvalueCell.setCellStyle(pvalue < 0.001 ? smallPvalueStyle : largePvalueStyle);

            XSSFCell qvalueCell = row.createCell(3 + maxAnnotations, CellType.NUMERIC);
            qvalueCell.setCellValue(qvalue);
            qvalueCell.setCellStyle(qvalue > 0 && qvalue < 0.001 ? smallPvalueStyle : largePvalueStyle);

            XSSFCell bonferroniCell = row.createCell(4 + maxAnnotations, CellType.BOOLEAN);
            bonferroniCell.setCellValue(pvalue <= bonferroniCutoff);

            XSSFCell fdrCell = row.createCell(5 + maxAnnotations, CellType.BOOLEAN);
            fdrCell.setCellValue(qvalue <= 0.05);

        }

        // Auto-scale collumns in sheet
        for (int c = 0; c < (5 + maxAnnotations); ++c) {
            databaseSheet.autoSizeColumn(c);
            databaseSheet.setColumnWidth(c, databaseSheet.getColumnWidth(c) + 1500); //compensate for with auto filter and inaccuracies
        }
    }

    private void populateGenePvalueSheet(Workbook enrichmentWorkbook, String trait, DoubleMatrixDataset<String, String> genePvalues) throws IOException {

        Map<String, Gene> geneInfo = IoUtils.readGenesMap(options.getGeneInfoFile());
        XSSFSheet genePSheet = (XSSFSheet) enrichmentWorkbook.createSheet("GenePvalues");
        XSSFTable table = genePSheet.createTable(new AreaReference(new CellReference(0, 0),
                new CellReference(genePvalues.rows(), 6),
                SpreadsheetVersion.EXCEL2007));

        table.setName("GenePvalues_res");
        table.setDisplayName("GenePvalues");
        table.setStyleName("TableStyleLight9");
        table.getCTTable().getTableStyleInfo().setShowRowStripes(true);
        table.getCTTable().addNewAutoFilter();
        XSSFRow headerRow = genePSheet.createRow(0);

        int hc = 0;
        headerRow.createCell(hc++, CellType.STRING).setCellValue("Gene id");
        headerRow.createCell(hc++, CellType.STRING).setCellValue("P-value");
        headerRow.createCell(hc++, CellType.STRING).setCellValue("Gene name");
        headerRow.createCell(hc++, CellType.STRING).setCellValue("Chromosome");
        headerRow.createCell(hc++, CellType.STRING).setCellValue("Start");
        headerRow.createCell(hc++, CellType.STRING).setCellValue("End");

        for (int r = 0; r < genePvalues.rows(); ++r) {

            XSSFRow row = genePSheet.createRow(r + 1); //+1 for header
            // gene id
            String geneId = genePvalues.getRowObjects().get(r);
            XSSFCell geneIdCell = row.createCell(0, CellType.STRING);
            geneIdCell.setCellValue(geneId);

            // p-value
            double pvalue = genePvalues.getElement(r, genePvalues.getColIndex(trait));
            XSSFCell pvalueCell = row.createCell(1, CellType.NUMERIC);
            pvalueCell.setCellValue(pvalue);
            pvalueCell.setCellStyle(pvalue < 0.001 ? smallPvalueStyle : largePvalueStyle);

            // gene name
            XSSFCell geneNameCell = row.createCell(2, CellType.STRING);
            geneNameCell.setCellValue(geneInfo.get(geneId).getGene());

            // chromosome
            XSSFCell chrCell = row.createCell(3, CellType.NUMERIC);
            chrCell.setCellValue(geneInfo.get(geneId).getSequenceName());

            // start
            XSSFCell startCell = row.createCell(4, CellType.NUMERIC);
            startCell.setCellValue(geneInfo.get(geneId).getStart());

            // end
            XSSFCell endCell = row.createCell(5, CellType.NUMERIC);
            endCell.setCellValue(geneInfo.get(geneId).getEnd());
        }

        // Auto-scale columns in sheet
        for (int c = 0; c < 6; ++c) {
            genePSheet.autoSizeColumn(c);
        }
    }

}


