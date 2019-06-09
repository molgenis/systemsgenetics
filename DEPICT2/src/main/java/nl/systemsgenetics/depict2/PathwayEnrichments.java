/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package nl.systemsgenetics.depict2;

import cern.colt.matrix.tdouble.DoubleMatrix1D;
import cern.colt.matrix.tdouble.DoubleMatrix2D;
import cern.colt.matrix.tdouble.algo.DenseDoubleAlgebra;
import com.opencsv.CSVWriter;
import edu.emory.mathcs.utils.ConcurrencyUtils;
import java.io.File;
import java.io.FileNotFoundException;
import java.io.FileOutputStream;
import java.io.FileWriter;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Collections;
import java.util.HashMap;
import java.util.HashSet;
import java.util.Iterator;
import java.util.LinkedHashMap;
import java.util.LinkedHashSet;
import java.util.List;
import java.util.Map;
import java.util.Set;
import me.tongfei.progressbar.ProgressBar;
import me.tongfei.progressbar.ProgressBarStyle;
import org.apache.log4j.Logger;
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
import umcg.genetica.math.matrix2.DoubleMatrixDatasetFastSubsetLoader;
import umcg.genetica.math.stats.ZScores;

/**
 *
 * @author patri
 */
public class PathwayEnrichments {

	private static final Logger LOGGER = Logger.getLogger(Depict2Options.class);

	public static HashMap<PathwayDatabase, DoubleMatrixDataset<String, String>> performEnrichmentAnalysis(final DoubleMatrixDataset<String, String> geneZscores, final DoubleMatrixDataset<String, String> geneZscoresNullGwas, final HashSet<String> genesWithPvalue, final List<PathwayDatabase> pathwayDatabases, List<Gene> genes, final String outputBasePath, HashSet<String> hlaGenesToExclude) throws IOException, Exception {

		ConcurrencyUtils.setNumberOfThreads(Depict2Options.getNumberOfThreadsToUse());

		final Set<String> excludeGenes;
		if (hlaGenesToExclude == null) {
			excludeGenes = Collections.emptySet();
		} else {
			excludeGenes = hlaGenesToExclude;
		}

		HashMap<PathwayDatabase, DoubleMatrixDataset<String, String>> enrichmentZscores = new HashMap<>(pathwayDatabases.size());

		//pathwayDatabases.parallelStream().forEach((PathwayDatabase pathwayDatabase) -> {
		for (PathwayDatabase pathwayDatabase : pathwayDatabases) {

			final DoubleMatrixDataset<String, String> genePathwayZscores;
			final DoubleMatrixDatasetFastSubsetLoader pathwayMatrixLoader;
			final LinkedHashSet<String> sharedGenes;

			if (pathwayDatabase.isTextBasedMatrix()) {
				LOGGER.info("WARNING using a debug function!");
				sharedGenes = new LinkedHashSet<>(DoubleMatrixDataset.loadDoubleTextData(pathwayDatabase.getLocation(), '\t').getHashRows().keySet());
				pathwayMatrixLoader = null;

			} else {

				pathwayMatrixLoader = new DoubleMatrixDatasetFastSubsetLoader(pathwayDatabase.getLocation());

				Set<String> pathwayGenes = pathwayMatrixLoader.getOriginalRowMap().keySet();

				sharedGenes = new LinkedHashSet<>();

				for (String gene : genesWithPvalue) {
					if (pathwayGenes.contains(gene) && !excludeGenes.contains(gene)) {
						sharedGenes.add(gene);
					}
				}

			}

			final Map<String, ArrayList<String>> geneChrArmMapping = createChrArmGeneMapping(genes, sharedGenes);

			try (ProgressBar pb = new ProgressBar(pathwayDatabase.getName() + " enrichtment analysis", geneChrArmMapping.size() + 2, ProgressBarStyle.ASCII)) {

				DoubleMatrixDataset<String, String> tmp = geneZscoresNullGwas.viewRowSelection(sharedGenes).duplicate();
				tmp.normalizeColumns();

				HashSet<String> sharedUncorrelatedGenes = findUncorrelatedGenes(tmp, sharedGenes, genes, 0.9);

				LOGGER.info("Number of uncorrelated genes: " + sharedUncorrelatedGenes.size());

				if (pathwayDatabase.isTextBasedMatrix()) {
					genePathwayZscores = DoubleMatrixDataset.loadDoubleTextData(pathwayDatabase.getLocation(), '\t');
				} else {
					genePathwayZscores = pathwayMatrixLoader.loadSubsetOfRowsBinaryDoubleData(sharedUncorrelatedGenes);
				}

				genePathwayZscores.normalizeColumns();

//				final DoubleMatrixDataset<String, String> geneZscoresPathwayMatched = geneZscores.viewRowSelection(sharedGenes).viewDice().createRowForceNormalDuplicate().viewDice();
//				final DoubleMatrixDataset<String, String> geneZscoresNullGwasPathwayMatched = geneZscoresNullGwas.viewRowSelection(sharedGenes).viewDice().createRowForceNormalDuplicate().viewDice();
				final DoubleMatrixDataset<String, String> geneZscoresPathwayMatched = geneZscores.viewRowSelection(sharedUncorrelatedGenes).duplicate();
				final DoubleMatrixDataset<String, String> geneZscoresNullGwasPathwayMatched = geneZscoresNullGwas.viewRowSelection(sharedUncorrelatedGenes).duplicate();

				geneZscoresPathwayMatched.normalizeColumns();
				geneZscoresNullGwasPathwayMatched.normalizeColumns();

//				{
//					DoubleMatrixDataset<String, String> geneZscoresPathwayMatchedTmp = geneZscores.viewRowSelection(sharedGenes).duplicate();
//					DoubleMatrixDataset<String, String> geneZscoresNullGwasPathwayMatchedTmp = geneZscoresNullGwas.viewRowSelection(sharedGenes).duplicate();
//
//					geneZscoresPathwayMatchedTmp.normalizeColumns();
//					geneZscoresPathwayMatched = geneZscoresPathwayMatchedTmp.viewDice().createRowForceNormalDuplicate().viewDice();
//
//					geneZscoresNullGwasPathwayMatchedTmp.normalizeColumns();
//					geneZscoresNullGwasPathwayMatched = geneZscoresNullGwasPathwayMatchedTmp.viewDice().createRowForceNormalDuplicate().viewDice();
//				}
				LinkedHashMap<String, Integer> singleColMap = new LinkedHashMap<>(1);
				singleColMap.put("B1", 0);

				//b1 rows: traits cols: 1
				final DoubleMatrixDataset<String, String> b1 = new DoubleMatrixDataset<>(geneZscores.getHashCols(), singleColMap);
				//b2 rows: traits cols: pathways
				final DoubleMatrixDataset<String, String> b2 = new DoubleMatrixDataset<>(geneZscores.getHashCols(), genePathwayZscores.getHashCols());

				final DoubleMatrixDataset<String, String> b1NullGwas = new DoubleMatrixDataset<>(geneZscoresNullGwas.getHashCols(), singleColMap);
				final DoubleMatrixDataset<String, String> b2NullGwas = new DoubleMatrixDataset<>(geneZscoresNullGwas.getHashCols(), genePathwayZscores.getHashCols());

				pb.step();

				//for (Map.Entry<String, ArrayList<String>> chrArmMappingEntry : geneChrArmMapping.entrySet()) {
				geneChrArmMapping.entrySet().parallelStream().forEach((Map.Entry<String, ArrayList<String>> chrArmMappingEntry) -> {

					try {
						final String chrArm = chrArmMappingEntry.getKey();
						final ArrayList<String> armGenes = chrArmMappingEntry.getValue();

						//invCorMatrixChrArm.setMatrix(DoubleFactory2D.dense.identity(invCorMatrixChrArm.rows()));
						final ArrayList<String> chrArmGenesInPathwayMatrix = new ArrayList<>(armGenes.size());

						for (String armGene : armGenes) {
							if (genePathwayZscores.containsRow(armGene)) {
								chrArmGenesInPathwayMatrix.add(armGene);
							}
						}
						//Now genesInPathwayMatrix will only contain genes that are also in the gene p-value matrix

						LOGGER.debug("Number of genes in chr arm: " + chrArmGenesInPathwayMatrix.size());

						if (!chrArmGenesInPathwayMatrix.isEmpty()) {

							final DoubleMatrixDataset<String, String> geneZscoresSubset = geneZscoresPathwayMatched.viewRowSelection(chrArmGenesInPathwayMatrix);
							final DoubleMatrixDataset<String, String> geneZscoresNullGwasSubset = geneZscoresNullGwasPathwayMatched.viewRowSelection(chrArmGenesInPathwayMatrix);
							final DoubleMatrixDataset<String, String> genePathwayZscoresSubset = genePathwayZscores.viewRowSelection(chrArmGenesInPathwayMatrix);

							DoubleMatrix2D geneInvCorMatrixSubset = new DenseDoubleAlgebra().inverse(geneZscoresNullGwasSubset.viewDice().calculateCorrelationMatrix().getMatrix());

							geneZscoresSubset.save(outputBasePath + "_" + pathwayDatabase.getName() + "_" + chrArm + "_Enrichment_genePvalues.txt");
							geneZscoresNullGwasSubset.save(outputBasePath + "_" + pathwayDatabase.getName() + "_" + chrArm + "_Enrichment_genePvaluesNull.txt");
							//geneInvCorMatrixSubset.save(outputBasePath + "_" + pathwayDatabase.getName() + "_" + chrArm + "_Enrichment_geneInvCor.txt");
							genePathwayZscoresSubset.save(outputBasePath + "_" + pathwayDatabase.getName() + "_" + chrArm + "_Enrichment_pathwayZscores.txt");

							//Inplace update of b1 and b2
							synchronized (b1) {
								glsStep1(geneZscoresSubset, geneInvCorMatrixSubset, genePathwayZscoresSubset, b1, b2);
							}
							synchronized (b1NullGwas) {
								glsStep1(geneZscoresNullGwasSubset, geneInvCorMatrixSubset, genePathwayZscoresSubset, b1NullGwas, b2NullGwas);
							}
//							b1.printMatrix();
//							b2.printMatrix();
//							geneZscoresPerChrArm.put(chrArm, geneZscoresSubset);
//							geneZscoresNullGwasPerChrArm.put(chrArm, geneZscoresNullGwas);
//							
//							//All matrices should now contain the same genes in the same order
//							//Calculate beta's using GLS using the precalculated inverse of the gene correlations
//							DoubleMatrixDataset<String, String> enrichmentChrArm = glsFixedInvCor(geneZscoresSubset, genePathwayZscoresSubset, geneInvCorMatrixSubset);
//							DoubleMatrixDataset<String, String> enrichmentNullChrArm = glsFixedInvCor(geneZscoresNullGwasSubset, genePathwayZscoresSubset, geneInvCorMatrixSubset);
//
//							enrichmentChrArm.save(outputBasePath + "_" + pathwayDatabase.getName() + "_" + chrArm + "_Enrichment_betas.txt");
//							enrichmentNullChrArm.save(outputBasePath + "_" + pathwayDatabase.getName() + "_" + chrArm + "_Enrichment_betasNull.txt");
//							//synchronized (enrichment) {
//							enrichment.getMatrix().assign(enrichmentChrArm.getMatrix(), DoubleFunctions.plus);
//							//}
//							//synchronized (enrichmentNull) {
//							enrichmentNull.getMatrix().assign(enrichmentNullChrArm.getMatrix(), DoubleFunctions.plus);
//							//}
						}

						pb.step();
					} catch (Exception ex) {
						throw new RuntimeException(ex);
					}

				});

				DoubleMatrixDataset<String, String> enrichment = new DoubleMatrixDataset<>(genePathwayZscores.getHashColsCopy(), geneZscores.getHashColsCopy());
				DoubleMatrixDataset<String, String> enrichmentNull = new DoubleMatrixDataset<>(genePathwayZscores.getHashColsCopy(), geneZscoresNullGwas.getHashColsCopy());

				final int numberOfPathways = genePathwayZscores.columns();
				final int numberTraits = geneZscores.columns();

				for (int traitI = 0; traitI < numberTraits; ++traitI) {
					final double b1Trait = b1.getElementQuick(traitI, 0);
					for (int pathwayI = 0; pathwayI < numberOfPathways; ++pathwayI) {

						double beta = b2.getElementQuick(traitI, pathwayI) / b1Trait;
						enrichment.setElementQuick(pathwayI, traitI, beta);

					}
				}

//					b1.printMatrix();
//					b2.printMatrix();
//					enrichment.printMatrix();
				final int numberTraitsNull = geneZscoresNullGwas.columns();

				for (int traitI = 0; traitI < numberTraitsNull; ++traitI) {
					final double b1Trait = b1NullGwas.getElementQuick(traitI, 0);
					for (int pathwayI = 0; pathwayI < numberOfPathways; ++pathwayI) {

						double beta = b2NullGwas.getElementQuick(traitI, pathwayI) / b1Trait;
						enrichmentNull.setElementQuick(pathwayI, traitI, beta);

					}
				}

//					enrichmentNull.printMatrix();
				enrichment.save(outputBasePath + "_" + pathwayDatabase.getName() + "_Enrichment" + (hlaGenesToExclude == null ? "_correlations" : "_correlationsExHla") + ".txt");
				enrichmentNull.save(outputBasePath + "_" + pathwayDatabase.getName() + "_EnrichmentNull" + (hlaGenesToExclude == null ? "_correlations" : "_correlationsExHla") + ".txt");

//					PathwayAnnotations pathwayAnnotations = new PathwayAnnotations(new File(pathwayDatabase.getLocation() + ".colAnnotations.txt"));
				//writeEnrichment(pathwayDatabase, pathwayAnnotations, enrichment, outputBasePath, hlaGenesToExclude == null ? "_correlations" : "_correlationsExHla");
				//writeEnrichment(pathwayDatabase, pathwayAnnotations, enrichmentNull, outputBasePath, hlaGenesToExclude == null ? "_correlations_null" : "_correlationsExHla_null");
				DoubleMatrix2D enrichmentMatrix = enrichment.getMatrix();
				DoubleMatrix2D enrichmentNullMatrix = enrichmentNull.getMatrix();

				final int numberOfPhenotypes = enrichmentMatrix.columns();
				final int numberOfNullGwasPhenotypes = enrichmentNullMatrix.columns();
//					final double numberOfNullGwasPhenotypesPlus1Double = enrichmentNullMatrix.columns() + 1;
//					final double minPvalue = 1 / numberOfNullGwasPhenotypesPlus1Double;
				final double numberOfNullGwasPhenotypesMin1Double = enrichmentNullMatrix.columns() - 1;

				LOGGER.debug("numberOfNullGwasPhenotypes: " + numberOfNullGwasPhenotypes);

				List<String> pathwayNames;
				if (LOGGER.isDebugEnabled()) {
					pathwayNames = enrichmentNull.getRowObjects();
				} else {
					pathwayNames = Collections.emptyList();
				}

//					for (int r = 0; r < numberOfPathways; ++r) {
//
//						for (int c = 0; c < numberOfPhenotypes; ++c) {
//
//							final double corr = Math.abs(enrichmentMatrix.getQuick(r, c));
//
//							int x = 0;
//							for (int p = 0; p < numberOfNullGwasPhenotypes; ++p) {
//								if (corr < Math.abs(enrichmentNullMatrix.getQuick(r, p))) {
//									++x;
//								}
//							}
//							if (x == 0) {
//								enrichmentMatrix.setQuick(r, c, minPvalue);
//							} else {
//								enrichmentMatrix.setQuick(r, c, ((x + 0.5) / numberOfNullGwasPhenotypesPlus1Double));
//							}
//
//						}
//
//					}
				for (int r = 0; r < numberOfPathways; ++r) {

					double meanNull = 0;
					for (int p = 0; p < numberOfNullGwasPhenotypes; ++p) {
						meanNull += enrichmentNullMatrix.getQuick(r, p);
					}
					meanNull /= numberOfNullGwasPhenotypes;

					double x = 0;
					for (int p = 0; p < numberOfNullGwasPhenotypes; ++p) {
						x += (enrichmentNullMatrix.getQuick(r, p) - meanNull) * (enrichmentNullMatrix.getQuick(r, p) - meanNull);
					}
					double sdNull = Math.sqrt(x / numberOfNullGwasPhenotypesMin1Double);

					if (LOGGER.isDebugEnabled()) {
						LOGGER.debug(pathwayNames.get(r) + " mean: " + meanNull + " sd: " + sdNull);
					}

					for (int c = 0; c < numberOfPhenotypes; ++c) {

						final double corr = enrichmentMatrix.getQuick(r, c);

						enrichmentMatrix.setQuick(r, c, (corr - meanNull) / sdNull);

					}
				}

				enrichment.printMatrix();

				//System.out.println(outputBasePath + "_" + pathwayDatabase.getName() + "_Enrichment" + (hlaGenesToExclude == null ? "_zscore" : "_zscoreExHla") + ".txt");
				enrichment.save(outputBasePath + "_" + pathwayDatabase.getName() + "_Enrichment" + (hlaGenesToExclude == null ? "_zscore" : "_zscoreExHla") + ".txt");

				enrichmentZscores.put(pathwayDatabase, enrichment);

				//writeEnrichment(pathwayDatabase, pathwayAnnotations, enrichment, outputBasePath, hlaGenesToExclude == null ? "_zscore" : "_zscoreExHla");
//					for (int r = 0; r < numberOfPathways; ++r) {
//
//						for (int c = 0; c < numberOfPhenotypes; ++c) {
//							
//							enrichmentMatrix.setQuick(r, c, ZScores.zToP(enrichmentMatrix.getQuick(r, c)));
//							
//						}
//					}
//
//					writeEnrichment(pathwayDatabase, pathwayAnnotations, enrichment, outputBasePath, hlaGenesToExclude == null ? "_pvalues" : "_pvaluesExHla");
				LOGGER.debug("Completed " + pathwayDatabase.getName() + " enrichment");

				//This extra step is intentional
				pb.step();

			}
		}

		return enrichmentZscores;

	}

	private static void writeEnrichment(PathwayDatabase pathwayDatabase, PathwayAnnotations pathwayAnnotations, DoubleMatrixDataset<String, String> enrichment, final String outputBasePath, final String nameSuffix) throws IOException {

		final CSVWriter enrichmentWriter = new CSVWriter(new FileWriter(new File(outputBasePath + "_" + pathwayDatabase.getName() + "_Enrichment" + nameSuffix + ".txt")), '\t', '\0', '\0', "\n");
		final String[] outputLine = new String[1 + pathwayAnnotations.getMaxNumberOfAnnotations() + enrichment.columns()];
		int c = 0;
		outputLine[c++] = "Pathway";
		for (int i = 0; i < pathwayAnnotations.getMaxNumberOfAnnotations(); ++i) {
			outputLine[c++] = "PathwayAnnotation" + (i + 1);
		}
		for (String col : enrichment.getHashCols().keySet()) {
			outputLine[c++] = col;
		}
		enrichmentWriter.writeNext(outputLine);

		for (String pathwayKey : enrichment.getHashRows().keySet()) {
			c = 0;
			outputLine[c++] = pathwayKey;

			if (pathwayAnnotations.getMaxNumberOfAnnotations() > 0) {
				ArrayList<String> thisPathwayAnnotations = pathwayAnnotations.getAnnotationsForPathway(pathwayKey);
				if (thisPathwayAnnotations == null) {
					for (int i = 0; i < pathwayAnnotations.getMaxNumberOfAnnotations(); ++i) {
						outputLine[c++] = "";
					}
				} else {
					for (int i = 0; i < pathwayAnnotations.getMaxNumberOfAnnotations(); ++i) {
						outputLine[c++] = i < thisPathwayAnnotations.size() ? thisPathwayAnnotations.get(i) : "";
					}
				}
			}

			DoubleMatrix1D row = enrichment.viewRow(pathwayKey);
			for (int e = 0; e < row.size(); ++e) {
				outputLine[c++] = String.valueOf(row.getQuick(e));
			}

			enrichmentWriter.writeNext(outputLine);
		}
		enrichmentWriter.close();

	}

	/**
	 *
	 * @param pathwayDatabases
	 * @param outputBasePath
	 * @param enrichments
	 * @param traits
	 * @param hlaExcluded
	 * @throws java.io.FileNotFoundException
	 */
	public static void saveEnrichmentsToExcel(final List<PathwayDatabase> pathwayDatabases, final String outputBasePath, HashMap<PathwayDatabase, DoubleMatrixDataset<String, String>> enrichments, List<String> traits, final boolean hlaExcluded) throws FileNotFoundException, IOException {

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

			for (PathwayDatabase pathwayDatabase : pathwayDatabases) {

				final PathwayAnnotations pathwayAnnotations = new PathwayAnnotations(new File(pathwayDatabase.getLocation() + ".colAnnotations.txt"));
				final int maxAnnotations = pathwayAnnotations.getMaxNumberOfAnnotations();
				final DoubleMatrixDataset<String, String> databaseEnrichment = enrichments.get(pathwayDatabase);
				final ArrayList<String> geneSets = databaseEnrichment.getRowObjects();
				//final int currentTraitCol = databaseEnrichment.getColIndex(trait);

				final DoubleMatrix1D traitEnrichment = databaseEnrichment.getCol(trait);
				int[] order = DoubleMatrix1dOrder.sortIndexReverse(traitEnrichment);
				XSSFSheet databaseSheet = (XSSFSheet) enrichmentWorkbook.createSheet(pathwayDatabase.getName());

				XSSFTable table = databaseSheet.createTable(new AreaReference(new CellReference(0, 0), new CellReference(databaseEnrichment.rows(), 2 + maxAnnotations), SpreadsheetVersion.EXCEL2007));
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

				for (int r = 0; r < databaseEnrichment.rows(); ++r) {
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

				}

				for (int c = 0; c < (3 + maxAnnotations); ++c) {
					databaseSheet.autoSizeColumn(c);
					databaseSheet.setColumnWidth(c, databaseSheet.getColumnWidth(c) + 1500);//compensate for with auto filter and inaccuracies
				}

			}

			File excelFile = new File(outputBasePath + "_enrichtments" + (traits.size() > 1 ? "_" + trait : "") + (hlaExcluded ? "_exHla.xlsx" : ".xlsx"));
			int nr = 1;
			while (excelFile.exists()) {
				excelFile = new File(outputBasePath + "_enrichtments" + (traits.size() > 1 ? "_" + trait : "") + (hlaExcluded ? "_exHla" : "" + "_" + nr + ".xlsx"));
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

			XSSFTable table = overviewSheet.createTable(new AreaReference(new CellReference(r, 0), new CellReference(r + pathwayDatabases.size(), 1), SpreadsheetVersion.EXCEL2007));
			table.setName("OverviewTable");
			table.setDisplayName("Overview");

			table.setStyleName("TableStyleLight9");
			table.getCTTable().getTableStyleInfo().setShowRowStripes(true);

			row = overviewSheet.createRow(r++);
			row.createCell(0, CellType.STRING).setCellValue("Gene set database");
			row.createCell(1, CellType.STRING).setCellValue("Number of sets");

			for (PathwayDatabase pathwayDatabase : pathwayDatabases) {
				row = overviewSheet.createRow(r++);
				cell = row.createCell(0, CellType.STRING);
				cell.setCellValue(pathwayDatabase.getName());

				Hyperlink link = createHelper.createHyperlink(HyperlinkType.DOCUMENT);
				link.setAddress(pathwayDatabase.getName() + "!A1");
				cell.setHyperlink(link);
				cell.setCellStyle(hlinkStyle);

				row.createCell(1, CellType.NUMERIC).setCellValue(enrichments.get(pathwayDatabase).rows());
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

//	private static void inplaceCorrectGenePathwayZscoresForWeights(DoubleMatrixDataset<String, String> genePathwayZscores, DoubleMatrixDataset<String, String> geneWeightsSubset) {
//
//		final DoubleMatrix2D genePathwayZscoresMatrix = genePathwayZscores.getMatrix();
//		final DoubleMatrix2D geneWeightsSubsetMatrix = geneWeightsSubset.getMatrix();
//
//		for (int pathway = 0; pathway < genePathwayZscoresMatrix.columns(); ++pathway) {
//
//			final SimpleRegression regression = new SimpleRegression();
//
//			for (int gene = 0; gene < genePathwayZscoresMatrix.rows(); ++gene) {
//
//				regression.addData(geneWeightsSubsetMatrix.getQuick(gene, 0), genePathwayZscoresMatrix.getQuick(gene, pathway));
//
//			}
//
//			final double b0 = regression.getIntercept();
//			final double b1 = regression.getSlope();
//
//			for (int gene = 0; gene < genePathwayZscoresMatrix.rows(); ++gene) {
//
//				//Residual = realValue - b0 + b1 * geneWeight
//				genePathwayZscoresMatrix.setQuick(gene, pathway, genePathwayZscoresMatrix.getQuick(gene, pathway) - (b0 + (b1 * geneWeightsSubsetMatrix.getQuick(gene, 0))));
//
//			}
//
//		}
//
//	}
	/**
	 * Note not a GLS over the full geneZscores matrix but per GWAS
	 *
	 * @param geneZscores
	 * @param genePathwayZscores
	 * @param geneInvCor
	 * @return
	 */
	protected static DoubleMatrixDataset<String, String> glsFixedInvCor(final DoubleMatrixDataset<String, String> geneZscores, final DoubleMatrixDataset<String, String> genePathwayZscores, final DoubleMatrixDataset<String, String> geneInvCor) {

		//Will contain the results rows: pathways, cols: traits
		final DoubleMatrixDataset<String, String> betas = new DoubleMatrixDataset<>(genePathwayZscores.getHashColsCopy(), geneZscores.getHashColsCopy());
		DoubleMatrix2D betasMatrix = betas.getMatrix();

		final int numberOfGenes = geneZscores.rows();
		final int numberTraits = geneZscores.columns();
		final int numberOfPathways = genePathwayZscores.columns();

		final DoubleMatrix2D geneInvCorMatrix = geneInvCor.getMatrix();
		final DoubleMatrix2D genePathwayZscoresMatrix = genePathwayZscores.getMatrix();

		//result of transpose geneZscoresTrait times inv correlation matrix
		DoubleMatrix2D A = geneZscores.getMatrix().like(1, numberOfGenes);
		DoubleMatrix2D B = geneZscores.getMatrix().like(1, numberOfPathways);

		//	final DenseDoubleAlgebra algebra = new DenseDoubleAlgebra();
		for (int traitI = 0; traitI < numberTraits; ++traitI) {

			try {

				DoubleMatrix2D geneZscoresTrait = geneZscores.viewColAsMmatrix(traitI);

				geneZscoresTrait.zMult(geneInvCorMatrix, A, 1, 0, true, false);

				double x = 1d / A.viewRow(0).zDotProduct(geneZscoresTrait.viewColumn(0));

				for (int g = 0; g < numberOfGenes; ++g) {
					A.setQuick(0, g, A.getQuick(0, g) * x);
				}

				A.zMult(genePathwayZscoresMatrix, B);

				betasMatrix.viewColumn(traitI).assign(B.viewRow(0));

			} catch (Exception e) {
				LOGGER.fatal("Number of pathways: " + numberOfPathways);
				LOGGER.fatal("Current trait index: " + traitI);
				LOGGER.fatal("Dim genePathwayZscores: " + genePathwayZscores.rows() + "x" + genePathwayZscores.columns());
				LOGGER.fatal("Dim genePathwayZscores internal: " + genePathwayZscores.getMatrix().rows() + "x" + genePathwayZscores.getMatrix().columns());

				throw (e);
			}
		}

		return betas;

	}

	private static void glsStep1(DoubleMatrixDataset<String, String> geneZscoresSubset, DoubleMatrix2D geneInvCorMatrix, DoubleMatrixDataset<String, String> genePathwayZscoresSubset, DoubleMatrixDataset<String, String> b1, DoubleMatrixDataset<String, String> b2) {

		final int numberOfGenes = geneZscoresSubset.rows();
		final int numberTraits = geneZscoresSubset.columns();
		final int numberOfPathways = genePathwayZscoresSubset.columns();

		//final DoubleMatrix2D geneInvCorMatrix = geneInvCorMatrixSubset.getMatrix();
		final DoubleMatrix2D genePathwayZscoresMatrix = genePathwayZscoresSubset.getMatrix();

		//result of transpose geneZscoresTrait times inv correlation matrix
		DoubleMatrix2D A = geneZscoresSubset.getMatrix().like(1, numberOfGenes);

		for (int traitI = 0; traitI < numberTraits; ++traitI) {

			try {

				DoubleMatrix2D geneZscoresTrait = geneZscoresSubset.viewColAsMmatrix(traitI);
				geneZscoresTrait.zMult(geneInvCorMatrix, A, 1, 0, true, false);

				final double x = A.viewRow(0).zDotProduct(geneZscoresTrait.viewColumn(0));

				//Col order should be the same
				b1.setElementQuick(0, traitI, x + b1.getElementQuick(0, traitI));

				DoubleMatrix2D b2Row = b2.viewRowAsMmatrix(traitI);

				//This does not clear b2 but instead will add and that is what we want now
				A.zMult(genePathwayZscoresMatrix, b2Row, 1, 1, false, false);

			} catch (Exception e) {
				LOGGER.fatal("Number of pathways: " + numberOfPathways);
				LOGGER.fatal("Current trait index: " + traitI);
				LOGGER.fatal("Dim genePathwayZscores: " + genePathwayZscoresSubset.rows() + "x" + genePathwayZscoresSubset.columns());
				LOGGER.fatal("Dim genePathwayZscores internal: " + genePathwayZscoresSubset.getMatrix().rows() + "x" + genePathwayZscoresSubset.getMatrix().columns());

				throw (e);
			}
		}

	}

	private static Map<String, ArrayList<String>> createChrArmGeneMapping(List<Gene> genes, Set<String> includedGenes) {
		Map<String, ArrayList<String>> chrArmToGeneMapping = new HashMap<>(25);
		for (Gene gene : genes) {

			if (includedGenes.contains(gene.getGene())) {

				String chrArm = gene.getChrAndArm();

				ArrayList<String> armGenes = chrArmToGeneMapping.get(chrArm);
				if (armGenes == null) {
					armGenes = new ArrayList<>();
					chrArmToGeneMapping.put(chrArm, armGenes);
				}

				armGenes.add(gene.getGene());

			}

		}
		return chrArmToGeneMapping;
	}

	private static HashSet<String> findUncorrelatedGenes(DoubleMatrixDataset<String, String> geneZscoresNullGwas, HashSet<String> genesWithPvalue, List<Gene> genes, double maxCorrelationBetweenGenes) {

		final Map<String, ArrayList<String>> chrArmToGeneMapping = createChrArmGeneMapping(genes, genesWithPvalue);

		final HashSet<String> includedUncorrelatedGenes = new HashSet<>();

		try (ProgressBar pb = new ProgressBar("Determining uncorrelated genes", chrArmToGeneMapping.size(), ProgressBarStyle.ASCII)) {

			chrArmToGeneMapping.keySet().parallelStream().forEach((String chrArm) -> {

				final ArrayList<String> armGenes = chrArmToGeneMapping.get(chrArm);

				final DoubleMatrixDataset<String, String> geneZscoresNullGwasArm = geneZscoresNullGwas.viewRowSelection(armGenes);

				//final DoubleMatrixDataset<String, String> invCorMatrixArmGenes = invCorMatrix.viewSelection(armGenes, armGenes);
				final DoubleMatrixDataset<String, String> genePvaluesNullGwasArmT = geneZscoresNullGwasArm.viewDice();

				final DoubleMatrixDataset<String, String> genePvaluesNullGwasGeneArmCorrelation = genePvaluesNullGwasArmT.calculateCorrelationMatrix();

				//We need to take the inverse of the correlation matrix. To do that the correlation between genes can't be correlated
				//Simply removing highly correlated genes did not always work, therefor:
				//(1) create correlation matrix of correlations
				//(2) identifie genes that have correlated correlation
				//(3) prune gene correlation matrix
				DoubleMatrixDataset<String, String> correlationOfCorrelations = genePvaluesNullGwasGeneArmCorrelation.calculateCorrelationMatrix();

				ArrayList<String> variantNames = correlationOfCorrelations.getRowObjects();

				rows:
				for (int r = 0; r < correlationOfCorrelations.rows(); ++r) {
					cols:
					for (int c = 0; c < r; ++c) {
						if (Math.abs(correlationOfCorrelations.getElementQuick(r, c)) >= maxCorrelationBetweenGenes && includedUncorrelatedGenes.contains(variantNames.get(c))) {
							continue rows;
						}
					}
					includedUncorrelatedGenes.add(variantNames.get(r));
				}

				pb.step();

			});
		}
		return includedUncorrelatedGenes;
	}
}
