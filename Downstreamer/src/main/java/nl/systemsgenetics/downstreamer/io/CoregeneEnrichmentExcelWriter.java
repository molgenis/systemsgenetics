/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package nl.systemsgenetics.downstreamer.io;

import cern.colt.matrix.tdouble.DoubleMatrix1D;
import java.io.File;
import java.io.FileNotFoundException;
import java.io.FileOutputStream;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.HashMap;
import java.util.List;
import nl.systemsgenetics.downstreamer.Downstreamer;
import nl.systemsgenetics.downstreamer.DownstreamerOptions;
import nl.systemsgenetics.downstreamer.pathway.PathwayAnnotations;
import nl.systemsgenetics.downstreamer.pathway.PathwayDatabase;
import org.apache.log4j.Logger;
import org.apache.poi.common.usermodel.HyperlinkType;
import org.apache.poi.ss.SpreadsheetVersion;
import org.apache.poi.ss.usermodel.CellType;
import org.apache.poi.ss.usermodel.CreationHelper;
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

/**
 *
 * @author patri
 */
public class CoregeneEnrichmentExcelWriter {

	private static final Logger LOGGER = Logger.getLogger(CoregeneEnrichmentExcelWriter.class);

	public static void write(
			final DownstreamerOptions options,
			final HashMap<String, DoubleMatrixDataset<String, String>> pathwayDatabase2Auc,
			final HashMap<String, DoubleMatrixDataset<String, String>> pathwayDatabase2Utest,
			final HashMap<String, DoubleMatrixDataset<String, String>> pathwayDatabase2BonfOdds,
			final HashMap<String, DoubleMatrixDataset<String, String>> pathwayDatabase2BonfFisherP,
			final HashMap<String, DoubleMatrixDataset<String, String>> pathwayDatabase2FdrOdds,
			final HashMap<String, DoubleMatrixDataset<String, String>> pathwayDatabase2FdrFisherP,
			final List<String> traits,
			List<PathwayDatabase> geneAnnotationDatabases
	) throws FileNotFoundException, IOException {

		final String outputBasePath = options.getOutputBasePath();
		final boolean hlaExcluded = options.isExcludeHla();

		System.setProperty("java.awt.headless", "true");

		for (String trait : traits) {

			Workbook wb = new XSSFWorkbook();
			ExcelStyles styles = new ExcelStyles(wb);
			CreationHelper createHelper = wb.getCreationHelper();

			populateOverviewSheet(wb, trait, geneAnnotationDatabases, createHelper, options, styles);

			for (PathwayDatabase geneAssociations : geneAnnotationDatabases) {

				LOGGER.debug("Creating sheet for: " + trait + " vs " + geneAssociations.getName());

				final List<String> geneSets = pathwayDatabase2Auc.get(geneAssociations.getName()).getRowObjects();

				final DoubleMatrix1D auc = pathwayDatabase2Auc.get(geneAssociations.getName()).getCol(trait);
				final DoubleMatrix1D uTestP = pathwayDatabase2Utest.get(geneAssociations.getName()).getCol(trait);
				final DoubleMatrix1D bonfOdds = pathwayDatabase2BonfOdds.get(geneAssociations.getName()).getCol(trait);
				final DoubleMatrix1D bonfFisherP = pathwayDatabase2BonfFisherP.get(geneAssociations.getName()).getCol(trait);
				final DoubleMatrix1D fdrOdss = pathwayDatabase2FdrOdds.get(geneAssociations.getName()).getCol(trait);
				final DoubleMatrix1D fdrFisherP = pathwayDatabase2FdrFisherP.get(geneAssociations.getName()).getCol(trait);

				int[] order = DoubleMatrix1dOrder.sortIndex(bonfFisherP);

				PathwayAnnotations geneAssociationsAnnotations = new PathwayAnnotations(new File(geneAssociations.getLocation() + ".colAnnotations.txt"));
				int maxAnnotations = geneAssociationsAnnotations.getMaxNumberOfAnnotations();

				XSSFSheet sh = (XSSFSheet) wb.createSheet(geneAssociations.getName());
				XSSFTable table = sh.createTable(new AreaReference(new CellReference(0, 0),
						new CellReference(geneSets.size(), 6 + maxAnnotations),
						SpreadsheetVersion.EXCEL2007));

				table.setName(geneAssociations.getName() + "_enrichment");
				table.setDisplayName(geneAssociations.getName());
				table.setStyleName("TableStyleLight9");
				table.getCTTable().getTableStyleInfo().setShowRowStripes(true);
				table.getCTTable().addNewAutoFilter();

				//databaseSheet.createFreezePane(0, 1);
				XSSFRow headerRow = sh.createRow(0);
				int hc = 0;
				headerRow.createCell(hc++, CellType.STRING).setCellValue(geneAssociationsAnnotations.getSetName() == null ? "Gene set" : geneAssociationsAnnotations.getSetName());
				for (int i = 0; i < maxAnnotations; ++i) {
					headerRow.createCell(hc++, CellType.STRING).setCellValue(geneAssociationsAnnotations.getAnnotationHeaders().get(i));
				}
				headerRow.createCell(hc++, CellType.STRING).setCellValue("Odds ratio of bonf sig genes");
				headerRow.createCell(hc++, CellType.STRING).setCellValue("Fisher's exact bonf sig genes");
				headerRow.createCell(hc++, CellType.STRING).setCellValue("Odds ratio of FDR 5% sig genes");
				headerRow.createCell(hc++, CellType.STRING).setCellValue("Fisher's exact FDR 5% sig genes");
				headerRow.createCell(hc++, CellType.STRING).setCellValue("AUC of priorization score");
				headerRow.createCell(hc++, CellType.STRING).setCellValue("Utest of priorization score");

				double v;
				XSSFCell cell;

				for (int r = 0; r < geneSets.size(); ++r) {

					XSSFRow row = sh.createRow(r + 1); //+1 for header
					String geneSet = geneSets.get(order[r]);
					row.createCell(0, CellType.STRING).setCellValue(geneSet);

					// Annotations from .colAnnotations file
					if (maxAnnotations > 0) {
						ArrayList<String> thisPathwayAnnotations = geneAssociationsAnnotations.getAnnotationsForPathway(geneSet);
						if (thisPathwayAnnotations == null) {
							for (int j = 0; j < maxAnnotations; ++j) {
								row.createCell(j + 1, CellType.STRING).setCellValue("");
							}
						} else {
							for (int j = 0; j < maxAnnotations; ++j) {
								if (j < thisPathwayAnnotations.size()) {
									String annotation = thisPathwayAnnotations.get(j);
									cell = row.createCell(j + 1, CellType.STRING);
									cell.setCellValue(annotation);

									if (annotation.startsWith("http")) {
										Hyperlink link = createHelper.createHyperlink(HyperlinkType.URL);
										link.setAddress(annotation);
										cell.setHyperlink(link);
										cell.setCellStyle(styles.getHlinkStyle());
									}

								} else {
									row.createCell(j + 1, CellType.STRING).setCellValue("");
								}

							}
						}
					}

					v = bonfOdds.getQuick(order[r]);
					cell = row.createCell(1 + maxAnnotations, CellType.NUMERIC);
					cell.setCellValue(v);
					cell.setCellStyle(styles.getZscoreStyle());

					v = bonfFisherP.getQuick(order[r]);
					cell = row.createCell(2 + maxAnnotations, CellType.NUMERIC);
					cell.setCellValue(v);
					cell.setCellStyle(v < 0.001 ? styles.getSmallPvalueStyle() : styles.getLargePvalueStyle());

					v = fdrOdss.getQuick(order[r]);
					cell = row.createCell(3 + maxAnnotations, CellType.NUMERIC);
					cell.setCellValue(v);
					cell.setCellStyle(styles.getZscoreStyle());

					v = fdrFisherP.getQuick(order[r]);
					cell = row.createCell(4 + maxAnnotations, CellType.NUMERIC);
					cell.setCellValue(v);
					cell.setCellStyle(v < 0.001 ? styles.getSmallPvalueStyle() : styles.getLargePvalueStyle());

					v = auc.getQuick(order[r]);
					cell = row.createCell(5 + maxAnnotations, CellType.NUMERIC);
					cell.setCellValue(v);
					cell.setCellStyle(styles.getZscoreStyle());

					v = uTestP.getQuick(order[r]);
					cell = row.createCell(6 + maxAnnotations, CellType.NUMERIC);
					cell.setCellValue(v);
					cell.setCellStyle(v < 0.001 ? styles.getSmallPvalueStyle() : styles.getLargePvalueStyle());

				}

				// Auto-scale columns in sheet
				for (int c = 0; c < hc; ++c) {
					sh.autoSizeColumn(c);
					sh.setColumnWidth(c, sh.getColumnWidth(c) + 1100); //compensate for with auto filter and inaccuracies
					if (c > 1 && sh.getColumnWidth(c) > 15000) {
						//max col width. Not for first column.
						sh.setColumnWidth(c, 15000);
					}
				}

			}

			File excelFile = new File(outputBasePath + "_prioritizedEnrichment" + (traits.size() > 1 ? "_" + trait : "") + (hlaExcluded ? "_exHla.xlsx" : ".xlsx"));
			int nr = 1;
			while (excelFile.exists()) {
				excelFile = new File(outputBasePath + "_prioritizedEnrichment" + (traits.size() > 1 ? "_" + trait : "") + (hlaExcluded ? "_exHla" : "") + "_" + nr + ".xlsx");
				nr++;
			}
			wb.write(new FileOutputStream(excelFile));

		}

	}

	private static void populateOverviewSheet(final Workbook wb, String trait, final List<PathwayDatabase> geneAnnotationDatabases, final CreationHelper createHelper, final DownstreamerOptions options, final ExcelStyles styles) {
		// -----------------------------------------------------------------------
		// Create overview sheet
		// -----------------------------------------------------------------------
		XSSFSheet overviewSheet = (XSSFSheet) wb.createSheet("Overview");

		int r = 0;
		XSSFRow row = overviewSheet.createRow(r++);
		XSSFCell cell = row.createCell(0, CellType.STRING);
		cell.setCellValue("Enrichments of prioritized genes for: " + trait);
		cell.setCellStyle(styles.getBoldStyle());

		row = overviewSheet.createRow(r++);
		cell = row.createCell(0, CellType.STRING);
		cell.setCellValue("Generated using Downstreamer " + Downstreamer.VERSION);
		cell.setCellStyle(styles.getBoldStyle());

		overviewSheet.createRow(r++);

		row = overviewSheet.createRow(r++);
		cell = row.createCell(0, CellType.STRING);
		cell.setCellValue("Gene set database");
		cell.setCellStyle(styles.getBoldStyle());
		cell = row.createCell(1, CellType.STRING);
		//cell.setCellValue("Number of sets");
		//cell.setCellStyle(styles.getBoldStyle());

		for (PathwayDatabase geneAnno : geneAnnotationDatabases) {
			row = overviewSheet.createRow(r++);
			cell = row.createCell(0, CellType.STRING);
			cell.setCellValue(geneAnno.getName());

			Hyperlink link = createHelper.createHyperlink(HyperlinkType.DOCUMENT);
			link.setAddress(geneAnno.getName() + "!A1");
			cell.setHyperlink(link);
			cell.setCellStyle(styles.getHlinkStyle());

			//row.createCell(1, CellType.NUMERIC).setCellValue(pathwayEnrichment.getNumberOfPathways());
		}

		for (int c = 0; c < 1; ++c) {
			overviewSheet.autoSizeColumn(c);
			overviewSheet.setColumnWidth(c, overviewSheet.getColumnWidth(c) + 1500);//compensate for with auto filter and inaccuracies
		}

		overviewSheet.createRow(r++);

		row = overviewSheet.createRow(r++);
		cell = row.createCell(0, CellType.STRING);
		cell.setCellValue("Used settings");
		cell.setCellStyle(styles.getBoldStyle());

		row = overviewSheet.createRow(r++);
		cell = row.createCell(0, CellType.STRING);
		cell.setCellValue("Number of permutations used for p-values: " + options.getPermutationPathwayEnrichment());

		row = overviewSheet.createRow(r++);
		cell = row.createCell(0, CellType.STRING);
		cell.setCellValue("Number of permutations used for FDR: " + options.getPermutationFDR());

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

}
