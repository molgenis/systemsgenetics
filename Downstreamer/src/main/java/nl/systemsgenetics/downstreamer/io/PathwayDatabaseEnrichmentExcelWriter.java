package nl.systemsgenetics.downstreamer.io;

import nl.systemsgenetics.downstreamer.runners.options.DownstreamerOptionsDeprecated;
import nl.systemsgenetics.downstreamer.pathway.PathwayAnnotations;
import nl.systemsgenetics.downstreamer.pathway.PathwayDatabase;
import org.apache.commons.lang3.StringUtils;
import org.apache.commons.lang3.math.NumberUtils;
import org.apache.logging.log4j.Logger;
import org.apache.logging.log4j.LogManager;
import org.apache.poi.common.usermodel.HyperlinkType;
import org.apache.poi.ss.SpreadsheetVersion;
import org.apache.poi.ss.usermodel.CellType;
import org.apache.poi.ss.usermodel.CreationHelper;
import org.apache.poi.ss.usermodel.Hyperlink;
import org.apache.poi.ss.usermodel.Workbook;
import org.apache.poi.ss.util.AreaReference;
import org.apache.poi.ss.util.CellReference;
import org.apache.poi.ss.util.WorkbookUtil;
import org.apache.poi.xssf.usermodel.*;
import nl.systemsgenetics.downstreamer.runners.PathwayDatabaseEnrichments.PathwayDatabaseEnrichmentRecord;

import java.io.File;
import java.io.FileOutputStream;
import java.io.IOException;
import java.util.*;

public class PathwayDatabaseEnrichmentExcelWriter {

    private static final Logger LOGGER = LogManager.getLogger(PathwayDatabaseEnrichmentExcelWriter.class);
    private DownstreamerOptionsDeprecated options;

    public PathwayDatabaseEnrichmentExcelWriter(DownstreamerOptionsDeprecated options) {
        this.options = options;
    }

    public void writeResultsExcel(Map<String, Map<String, List<PathwayDatabaseEnrichmentRecord>>> results, String predictionSource) throws IOException {

        String outputBasePath = options.getOutputBasePath();
        Set<String> gwasTraits = results.keySet();

        System.setProperty("java.awt.headless", "true");

        for (String gwasTrait : gwasTraits) {

            Workbook wb = new XSSFWorkbook();
            ExcelStyles styles = new ExcelStyles(wb);
            CreationHelper createHelper = wb.getCreationHelper();

            //populateOverviewSheet(wb, gwasTrait, geneAnnotationDatabases, createHelper, options, styles, predictionSource);

            // Loop over all the target databases that have been tested against
            for (PathwayDatabase geneAssociations : options.getPathwayDatabases2()) {

                // Retrieve the records to write in this sheet
                List<PathwayDatabaseEnrichmentRecord> curRecords = results.get(gwasTrait).get(geneAssociations.getName());

                // Sort the records on the bonferoni fisher exact pvalue
                curRecords.sort(Comparator.comparing(PathwayDatabaseEnrichmentRecord::getBonfPvalue));

//                // Read the (optional) annotations
                String geneAssociationAnnotationFile = geneAssociations.getLocation() + ".colAnnotations.txt";
                if(!new File(geneAssociationAnnotationFile).canRead()){
                    geneAssociationAnnotationFile+=".gz";
                    if(!new File(geneAssociationAnnotationFile).canRead()){
                        LOGGER.debug("Cannot find file: "+geneAssociations.getLocation() + ".colAnnotations.txt or "+geneAssociations.getLocation() + ".colAnnotations.txt.gz");
                    }
                }

                PathwayAnnotations geneAssociationsAnnotations = new PathwayAnnotations(new File(geneAssociationAnnotationFile));

                int maxAnnotations = geneAssociationsAnnotations.getMaxNumberOfAnnotations();

                // Initialize the sheet
                XSSFSheet sh = (XSSFSheet) wb.createSheet(WorkbookUtil.createSafeSheetName(geneAssociations.getName()));
                XSSFTable table = sh.createTable(new AreaReference(new CellReference(0, 0),
                        new CellReference(curRecords.size(), 11 + maxAnnotations),
                        SpreadsheetVersion.EXCEL2007));

                String tableName = geneAssociations.getName();
                tableName = tableName.replace('-', '_');

                table.setName(tableName + "_enrichment");
                table.setDisplayName(tableName);
                table.setStyleName("TableStyleLight9");
                table.getCTTable().getTableStyleInfo().setShowRowStripes(true);
                table.getCTTable().addNewAutoFilter();

                // Make the header
                XSSFRow headerRow = sh.createRow(0);
                int hc = 0;

                // Header for geneset names
                headerRow.createCell(hc++, CellType.STRING).setCellValue(geneAssociationsAnnotations.getSetName() == null ? "Gene set" : geneAssociationsAnnotations.getSetName());

                // Header names of optional annotations
                for (int i = 0; i < maxAnnotations; ++i) {
                    headerRow.createCell(hc++, CellType.STRING).setCellValue(geneAssociationsAnnotations.getAnnotationHeaders().get(i));
                }

                // Other columns
                headerRow.createCell(hc++, CellType.STRING).setCellValue("# genes in pathway");
                headerRow.createCell(hc++, CellType.STRING).setCellValue("overlap bonf");
                headerRow.createCell(hc++, CellType.STRING).setCellValue("OR bonf");
                headerRow.createCell(hc++, CellType.STRING).setCellValue("P bonf");
                headerRow.createCell(hc++, CellType.STRING).setCellValue("overlap FDR");
                headerRow.createCell(hc++, CellType.STRING).setCellValue("OR FDR");
                headerRow.createCell(hc++, CellType.STRING).setCellValue("P FDR");
                headerRow.createCell(hc++, CellType.STRING).setCellValue("AUC");
                headerRow.createCell(hc++, CellType.STRING).setCellValue("Utest");
                headerRow.createCell(hc++, CellType.STRING).setCellValue("Bonf overlapping genes");
                headerRow.createCell(hc++, CellType.STRING).setCellValue("FDR overlapping genes");

                double v;
                XSSFCell cell;
                int r = 0;

                for (PathwayDatabaseEnrichmentRecord curRecord : curRecords) {

                    int c = 0;

                    XSSFRow row = sh.createRow(r + 1); //+1 for header
                    String geneSet = curRecord.getPathwayName();
                    row.createCell(c++, CellType.STRING).setCellValue(geneSet);

                    // Annotations from .colAnnotations file
                    if (maxAnnotations > 0) {
                        ArrayList<String> thisPathwayAnnotations = geneAssociationsAnnotations.getAnnotationsForPathway(geneSet);
                        if (thisPathwayAnnotations == null) {
                            for (int j = 0; j < maxAnnotations; ++j) {
                                row.createCell(c++, CellType.STRING).setCellValue("");
                            }
                        } else {
                            for (int j = 0; j < maxAnnotations; ++j) {
                                if (j < thisPathwayAnnotations.size()) {

                                    String annotation = thisPathwayAnnotations.get(j);

                                    if (NumberUtils.isCreatable(annotation)) {
                                        cell = row.createCell(c++, CellType.NUMERIC);
                                        cell.setCellValue(Double.parseDouble(annotation));
                                    } else {
                                        cell = row.createCell(c++, CellType.STRING);
                                        cell.setCellValue(annotation);

                                        if (annotation.startsWith("http")) {
                                            Hyperlink link = createHelper.createHyperlink(HyperlinkType.URL);
                                            link.setAddress(annotation);
                                            cell.setHyperlink(link);
                                            cell.setCellStyle(styles.getHlinkStyle());
                                        }
                                    }

                                } else {
                                    row.createCell(c++, CellType.STRING).setCellValue("");
                                }

                            }
                        }
                    }

                    // Genes in the pathway
                    v = curRecord.getBonfSigResult().getGenesInPathway().size();
                    cell = row.createCell(c++, CellType.NUMERIC);
                    cell.setCellValue(v);
                    cell.setCellStyle(styles.getIntStyle());

                    // Bonf overlapping Genes
                    v = curRecord.getBonfSigResult().getSignificantOverlappingGenes().size();
                    cell = row.createCell(c++, CellType.NUMERIC);
                    cell.setCellValue(v);
                    cell.setCellStyle(styles.getIntStyle());

                    // Bonf OR
                    v = curRecord.getBonfSigResult().getOddsRatio();
                    cell = row.createCell(c++, CellType.NUMERIC);
                    cell.setCellValue(v);
                    cell.setCellStyle(styles.getZscoreStyle());

                    // Bonf Pval
                    v = curRecord.getBonfSigResult().getpValue();
                    cell = row.createCell(c++, CellType.NUMERIC);
                    cell.setCellValue(v);
                    cell.setCellStyle(v < 0.001 ? styles.getSmallPvalueStyle() : styles.getLargePvalueStyle());

                    // FDR overlapping Genes
                    v = curRecord.getFdrSigResult().getSignificantOverlappingGenes().size();
                    cell = row.createCell(c++, CellType.NUMERIC);
                    cell.setCellValue(v);
                    cell.setCellStyle(styles.getIntStyle());

                    // Bonf OR
                    v = curRecord.getFdrSigResult().getOddsRatio();
                    cell = row.createCell(c++, CellType.NUMERIC);
                    cell.setCellValue(v);
                    cell.setCellStyle(styles.getZscoreStyle());

                    // Bonf Pval
                    v = curRecord.getFdrSigResult().getpValue();
                    cell = row.createCell(c++, CellType.NUMERIC);
                    cell.setCellValue(v);
                    cell.setCellStyle(v < 0.001 ? styles.getSmallPvalueStyle() : styles.getLargePvalueStyle());

                    // AUC
                    v = curRecord.getAuc();
                    cell = row.createCell(c++, CellType.NUMERIC);
                    cell.setCellValue(v);
                    cell.setCellStyle(styles.getZscoreStyle());

                    // AUC pval
                    v = curRecord.getAucPvalue();
                    cell = row.createCell(c++, CellType.NUMERIC);
                    cell.setCellValue(v);
                    cell.setCellStyle(v < 0.001 ? styles.getSmallPvalueStyle() : styles.getLargePvalueStyle());

                    // Bonferroni signficant genes overlapping the pathway
                    String geneList = StringUtils.join(curRecord.getBonfSigResult().getSignificantOverlappingGenes(), ' ');
                    if (geneList.length() >= 32766 || curRecord.getBonfSigResult().getSignificantOverlappingGenes().size() >= 500) {
                        geneList = "TooManyGenes";
                    }

                    cell = row.createCell(c++, CellType.STRING);
                    cell.setCellValue(geneList);

                    // FDR signficant genes overlapping the pathway
                    geneList = StringUtils.join(curRecord.getFdrSigResult().getSignificantOverlappingGenes(), ' ');
                    cell = row.createCell(c++, CellType.STRING);
                    if (geneList.length() >= 32766 || curRecord.getBonfSigResult().getSignificantOverlappingGenes().size() >= 500) {
                        geneList = "TooManyGenes";
                    }
                    cell.setCellValue(geneList);

                    // Advance row one.
                    r++;
                }

                // Auto-scale columns in sheet
                for (int c = 0; c < hc; ++c) {
                    sh.autoSizeColumn(c);
                    //compensate for with auto filter and inaccuracies;
                    int cw = sh.getColumnWidth(c) + 800;
                    if (cw >= 65280) {
                        cw = 15000;
                    }
                    sh.setColumnWidth(c, cw);
                    if (c >= 1 && sh.getColumnWidth(c) > 15000) {
                        // max col width. Not for first column.
                        sh.setColumnWidth(c, 15000);
                    }
                }

            }

            // Write the excel file
            File excelFile = new File(outputBasePath + "_" + predictionSource + "_Enrichment" + (gwasTraits.size() > 1 ? "_" + gwasTrait : "") + ".xlsx");
            int nr = 1;
            while (excelFile.exists()) {
                excelFile = new File(outputBasePath + "_" + predictionSource + "_Enrichment" + (gwasTraits.size() > 1 ? "_" + gwasTrait : "") + "_" + nr + ".xlsx");
                nr++;
            }
            wb.write(new FileOutputStream(excelFile));

        }
    }


}


