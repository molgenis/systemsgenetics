/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package nl.systemsgenetics.depict2;

import cern.colt.matrix.tdouble.DoubleMatrix1D;
import java.io.File;
import java.io.FileNotFoundException;
import java.io.FileOutputStream;
import java.io.IOException;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;
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
import umcg.genetica.math.matrix2.DoubleMatrix1dOrder;
import umcg.genetica.math.matrix2.DoubleMatrixDataset;
import umcg.genetica.math.stats.ZScores;

/**
 *
 * @author patri
 */
public class ExcelWriter {
	
	/**
	 *
	 * @param pathwayEnrichments
	 * @param pathwayDatabases
	 * @param outputBasePath
	 * @param traits
	 * @param hlaExcluded
	 * @throws java.io.FileNotFoundException
	 */
	public static void saveEnrichmentsToExcel(final List<PathwayEnrichments> pathwayEnrichments, final String outputBasePath, List<String> traits, final boolean hlaExcluded) throws FileNotFoundException, IOException {

		System.setProperty(" java.awt.headless", "true");

		for (String trait : traits) {

			Workbook enrichmentWorkbook = new XSSFWorkbook();
			DataFormat format = enrichmentWorkbook.createDataFormat();
			CreationHelper createHelper = enrichmentWorkbook.getCreationHelper();

			CellStyle zscoreStyle = enrichmentWorkbook.createCellStyle();
			zscoreStyle.setDataFormat(format.getFormat("0.00"));

			CellStyle largePvalueStyle = enrichmentWorkbook.createCellStyle();
			largePvalueStyle.setDataFormat(format.getFormat("0.0000"));

			CellStyle smallPvalueStyle = enrichmentWorkbook.createCellStyle();
			smallPvalueStyle.setDataFormat(format.getFormat("0.00E+0"));

			CellStyle hlinkStyle = enrichmentWorkbook.createCellStyle();
			Font hlinkFont = enrichmentWorkbook.createFont();
			hlinkFont.setUnderline(Font.U_SINGLE);
			hlinkStyle.setFont(hlinkFont);

			CellStyle boldStyle = enrichmentWorkbook.createCellStyle();
			Font fontBold = enrichmentWorkbook.createFont();
			fontBold.setBold(true);
			boldStyle.setFont(fontBold);

			XSSFSheet overviewSheet = (XSSFSheet) enrichmentWorkbook.createSheet("Overview");

			for (PathwayEnrichments pathwayEnrichment : pathwayEnrichments) {
				
				PathwayDatabase pathwayDatabase = pathwayEnrichment.getPathwayDatabase();
				
				final PathwayAnnotations pathwayAnnotations = new PathwayAnnotations(new File(pathwayDatabase.getLocation() + ".colAnnotations.txt"));
				final int maxAnnotations = pathwayAnnotations.getMaxNumberOfAnnotations();
				final DoubleMatrixDataset<String, String> databaseEnrichmentZscores = pathwayEnrichment.getEnrichmentZscores();
				final ArrayList<String> geneSets = databaseEnrichmentZscores.getRowObjects();
				//final int currentTraitCol = databaseEnrichment.getColIndex(trait);
				
				final double bonferroniCutoff = 0.05 / pathwayEnrichment.getNumberOfPathways();

				final DoubleMatrix1D traitEnrichment = databaseEnrichmentZscores.getCol(trait);
				int[] order = DoubleMatrix1dOrder.sortIndexReverse(traitEnrichment);
				XSSFSheet databaseSheet = (XSSFSheet) enrichmentWorkbook.createSheet(pathwayDatabase.getName());

				XSSFTable table = databaseSheet.createTable(new AreaReference(new CellReference(0, 0), new CellReference(databaseEnrichmentZscores.rows(), 3 + maxAnnotations), SpreadsheetVersion.EXCEL2007));
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
				headerRow.createCell(hc++, CellType.STRING).setCellValue("Bonferroni significant");

				for (int r = 0; r < databaseEnrichmentZscores.rows(); ++r) {
					XSSFRow row = databaseSheet.createRow(r + 1);//+1 for header

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

					XSSFCell zscoreCell = row.createCell(1 + maxAnnotations, CellType.NUMERIC);
					zscoreCell.setCellValue(zscore);
					zscoreCell.setCellStyle(zscoreStyle);

					double pvalue = ZScores.zToP(zscore);

					XSSFCell pvalueCell = row.createCell(2 + maxAnnotations, CellType.NUMERIC);
					pvalueCell.setCellValue(pvalue);
					pvalueCell.setCellStyle(pvalue < 0.001 ? smallPvalueStyle : largePvalueStyle);
					
					XSSFCell bonferroniCell = row.createCell(3 + maxAnnotations, CellType.BOOLEAN);
					bonferroniCell.setCellValue(pvalue <= bonferroniCutoff);

				}

				for (int c = 0; c < (3 + maxAnnotations); ++c) {
					databaseSheet.autoSizeColumn(c);
					databaseSheet.setColumnWidth(c, databaseSheet.getColumnWidth(c) + 1500);//compensate for with auto filter and inaccuracies
				}

			}

			File excelFile = new File(outputBasePath + "_enrichtments" + (traits.size() > 1 ? "_" + trait : "") + (hlaExcluded ? "_exHla.xlsx" : ".xlsx"));
			int nr = 1;
			while (excelFile.exists()) {
				excelFile = new File(outputBasePath + "_enrichtments" + (traits.size() > 1 ? "_" + trait : "") + (hlaExcluded ? "_exHla" : "") + "_" + nr + ".xlsx");
				nr++;
			}

			int r = 0;
			XSSFRow row = overviewSheet.createRow(r++);
			XSSFCell cell = row.createCell(0, CellType.STRING);
			cell.setCellValue("Pathway enrichment analysis for: " + trait);
			cell.setCellStyle(boldStyle);

			row = overviewSheet.createRow(r++);
			cell = row.createCell(0, CellType.STRING);
			cell.setCellValue("Generated using DEPICT" + Depict2.VERSION);
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
				overviewSheet.setColumnWidth(c, overviewSheet.getColumnWidth(c) + 1500);//compensate for with auto filter and inaccuracies
			}

			enrichmentWorkbook.write(new FileOutputStream(excelFile));

			System.err.println("WARNING ONLY SAVING FIRST TRAIT TO EXCEL FOR DEBUGING");
			break;

		}
	}
	
}
