/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package nl.systemsgenetics.depict2;

import cern.colt.matrix.tdouble.DoubleFactory2D;
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
import org.apache.commons.math3.stat.correlation.PearsonsCorrelation;
import org.apache.commons.math3.stat.regression.SimpleRegression;
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

	public static HashMap<PathwayDatabase, DoubleMatrixDataset<String, String>> performEnrichmentAnalysis(final DoubleMatrixDataset<String, String> geneZscores, final DoubleMatrixDataset<String, String> geneZscoresNullGwas, final HashSet<String> genesWithPvalue, final List<PathwayDatabase> pathwayDatabases, List<Gene> genes, final String outputBasePath, final HashSet<String> hlaGenesToExclude, final int nrSampleToUseForCorrelation, final int nrSamplesToUseForNullBetas, final double genePruningR, final boolean ignoreGeneCorrelations, final boolean forceNormalGenePvalues, final boolean forceNormalPathwayPvalues, final int correlationWindow) throws IOException, Exception {

		ConcurrencyUtils.setNumberOfThreads(Depict2Options.getNumberOfThreadsToUse());

		final Set<String> nullGwasRuns = geneZscoresNullGwas.getHashCols().keySet();
		if (nullGwasRuns.size() < (nrSampleToUseForCorrelation + nrSamplesToUseForNullBetas)) {
			throw new Exception("Not enough null gwas runs: " + nullGwasRuns.size() + " < " + nrSampleToUseForCorrelation + " + " + nrSamplesToUseForNullBetas);
		}

		Iterator<String> nullGwasRunIterator = nullGwasRuns.iterator();

		final LinkedHashSet<String> sampleToUseForCorrelation = new LinkedHashSet<>(nrSampleToUseForCorrelation);
		for (int i = 0; i < nrSampleToUseForCorrelation; ++i) {
			sampleToUseForCorrelation.add(nullGwasRunIterator.next());
		}

		final LinkedHashSet<String> samplesToUseForNullBetas = new LinkedHashSet<>(nrSamplesToUseForNullBetas);
		for (int i = 0; i < nrSamplesToUseForNullBetas; ++i) {
			samplesToUseForNullBetas.add(nullGwasRunIterator.next());
		}

		final DoubleMatrixDataset<String, String> geneZscoresNullGwasCorrelation = geneZscoresNullGwas.viewColSelection(sampleToUseForCorrelation);
		final DoubleMatrixDataset<String, String> geneZscoresNullGwasNullBetas = geneZscoresNullGwas.viewColSelection(samplesToUseForNullBetas);

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

			pathwayMatrixLoader = new DoubleMatrixDatasetFastSubsetLoader(pathwayDatabase.getLocation());

			Set<String> pathwayGenes = pathwayMatrixLoader.getOriginalRowMap().keySet();

			sharedGenes = new LinkedHashSet<>();

			for (String gene : genesWithPvalue) {
				if (pathwayGenes.contains(gene) && !excludeGenes.contains(gene)) {
					sharedGenes.add(gene);
				}
			}

			final Map<String, ArrayList<Gene>> geneChrArmMapping = createChrArmGeneMapping(genes, sharedGenes);

			try (ProgressBar pb = new ProgressBar(pathwayDatabase.getName() + " enrichtment analysis", geneChrArmMapping.size() + 2, ProgressBarStyle.ASCII)) {

				DoubleMatrixDataset<String, String> tmp = geneZscoresNullGwasCorrelation.viewRowSelection(sharedGenes).duplicate();
				tmp.normalizeColumns();

				LinkedHashSet<String> sharedUncorrelatedGenes = findUncorrelatedGenes(tmp, sharedGenes, genes, genePruningR);

				LOGGER.debug("Number of uncorrelated genes: " + sharedUncorrelatedGenes.size());

				if (forceNormalPathwayPvalues) {
					genePathwayZscores = pathwayMatrixLoader.loadSubsetOfRowsBinaryDoubleData(sharedUncorrelatedGenes).createColumnForceNormalDuplicate();
				} else {
					genePathwayZscores = pathwayMatrixLoader.loadSubsetOfRowsBinaryDoubleData(sharedUncorrelatedGenes);
				}

				genePathwayZscores.normalizeColumns();

				final DoubleMatrixDataset<String, String> geneZscoresPathwayMatched;
				final DoubleMatrixDataset<String, String> geneZscoresNullGwasCorrelationPathwayMatched;
				final DoubleMatrixDataset<String, String> geneZscoresNullGwasNullBetasPathwayMatched;

				if (forceNormalGenePvalues) {
					geneZscoresPathwayMatched = geneZscores.viewRowSelection(sharedUncorrelatedGenes).createColumnForceNormalDuplicate();
					geneZscoresNullGwasCorrelationPathwayMatched = geneZscoresNullGwasCorrelation.viewRowSelection(sharedUncorrelatedGenes).createColumnForceNormalDuplicate();
					geneZscoresNullGwasNullBetasPathwayMatched = geneZscoresNullGwasNullBetas.viewRowSelection(sharedUncorrelatedGenes).createColumnForceNormalDuplicate();
				} else {
					geneZscoresPathwayMatched = geneZscores.viewRowSelection(sharedUncorrelatedGenes).duplicate();
					geneZscoresNullGwasCorrelationPathwayMatched = geneZscoresNullGwasCorrelation.viewRowSelection(sharedUncorrelatedGenes).duplicate();
					geneZscoresNullGwasNullBetasPathwayMatched = geneZscoresNullGwasNullBetas.viewRowSelection(sharedUncorrelatedGenes).duplicate();
				}

				geneZscoresPathwayMatched.normalizeColumns();
				geneZscoresNullGwasCorrelationPathwayMatched.normalizeColumns();
				geneZscoresNullGwasNullBetasPathwayMatched.normalizeColumns();
				
				genePathwayZscores.save(outputBasePath + "_" + pathwayDatabase.getName() + "_Enrichment_normalizedPathwayScores" + (hlaGenesToExclude == null ? "" : "_ExHla") + ".txt");
				geneZscoresPathwayMatched.save(outputBasePath + "_" + pathwayDatabase.getName() + "_Enrichment_normalizedGwasGeneScores" + (hlaGenesToExclude == null ? "" : "_ExHla") + ".txt");

				LinkedHashMap<String, Integer> singleColMap = new LinkedHashMap<>(1);
				singleColMap.put("B1", 0);

				//b1 rows: traits cols: 1
				final DoubleMatrixDataset<String, String> b1 = new DoubleMatrixDataset<>(geneZscores.getHashCols(), singleColMap);
				//b2 rows: traits cols: pathways
				final DoubleMatrixDataset<String, String> b2 = new DoubleMatrixDataset<>(geneZscores.getHashCols(), genePathwayZscores.getHashCols());

				final DoubleMatrixDataset<String, String> b1NullGwas = new DoubleMatrixDataset<>(geneZscoresNullGwasNullBetas.getHashCols(), singleColMap);
				final DoubleMatrixDataset<String, String> b2NullGwas = new DoubleMatrixDataset<>(geneZscoresNullGwasNullBetas.getHashCols(), genePathwayZscores.getHashCols());

				pb.step();

				//for (Map.Entry<String, ArrayList<String>> chrArmMappingEntry : geneChrArmMapping.entrySet()) {
				geneChrArmMapping.entrySet().parallelStream().forEach((Map.Entry<String, ArrayList<Gene>> chrArmMappingEntry) -> {

					try {
						final String chrArm = chrArmMappingEntry.getKey();
						final ArrayList<Gene> armGenes = chrArmMappingEntry.getValue();

						final ArrayList<String> chrArmGenesInPathwayMatrix = new ArrayList<>(armGenes.size());

						for (Gene armGene : armGenes) {
							if (genePathwayZscores.containsRow(armGene.getGene())) {
								chrArmGenesInPathwayMatrix.add(armGene.getGene());
							}
						}
						//Now genesInPathwayMatrix will only contain genes that are also in the gene p-value matrix

						LOGGER.debug("Number of genes in chr arm: " + chrArmGenesInPathwayMatrix.size());

						if (!chrArmGenesInPathwayMatrix.isEmpty()) {

							final DoubleMatrixDataset<String, String> geneZscoresSubset = geneZscoresPathwayMatched.viewRowSelection(chrArmGenesInPathwayMatrix);
							final DoubleMatrixDataset<String, String> geneZscoresNullGwasCorrelationSubset = geneZscoresNullGwasCorrelationPathwayMatched.viewRowSelection(chrArmGenesInPathwayMatrix);
							final DoubleMatrixDataset<String, String> geneZscoresNullGwasNullBetasSubset = geneZscoresNullGwasNullBetasPathwayMatched.viewRowSelection(chrArmGenesInPathwayMatrix);
							final DoubleMatrixDataset<String, String> genePathwayZscoresSubset = genePathwayZscores.viewRowSelection(chrArmGenesInPathwayMatrix);

							final DoubleMatrixDataset<String, String> geneZscoresNullGwasSubsetGeneCorrelations;
							if(correlationWindow < 0){
								geneZscoresNullGwasSubsetGeneCorrelations = createLocalGeneCorrelation(geneZscoresNullGwasCorrelationSubset, armGenes, correlationWindow); 
							} else {
								geneZscoresNullGwasSubsetGeneCorrelations = geneZscoresNullGwasCorrelationSubset.viewDice().calculateCorrelationMatrix();
							}
								
//							geneZscoresSubset.save(outputBasePath + "_" + pathwayDatabase.getName() + "_" + chrArm + "_Enrichment_genePvalues.txt");
//							geneZscoresNullGwasSubsetGeneCorrelations.save(outputBasePath + "_" + pathwayDatabase.getName() + "_" + chrArm + "_Enrichment_geneCor.txt");
//							genePathwayZscoresSubset.save(outputBasePath + "_" + pathwayDatabase.getName() + "_" + chrArm + "_Enrichment_pathwayZscores.txt");
							final DoubleMatrix2D geneInvCorMatrixSubsetMatrix;
							try {
								if (ignoreGeneCorrelations) {
									geneInvCorMatrixSubsetMatrix = DoubleFactory2D.dense.identity(geneZscoresNullGwasSubsetGeneCorrelations.rows());
								} else {
									geneInvCorMatrixSubsetMatrix = new DenseDoubleAlgebra().inverse(geneZscoresNullGwasSubsetGeneCorrelations.getMatrix());
								}
							} catch (Exception ex) {
								LOGGER.fatal(pathwayDatabase.getName() + " " + chrArm);
								throw ex;
							}

							//final DoubleMatrixDataset<String, String> geneInvCorMatrixSubset = new DoubleMatrixDataset<>(geneInvCorMatrixSubsetMatrix, geneZscoresNullGwasNullBetasSubset.getHashRows(), geneZscoresNullGwasNullBetasSubset.getHashRows());
							//geneInvCorMatrixSubset.save(outputBasePath + "_" + pathwayDatabase.getName() + "_" + chrArm + "_Enrichment_geneInvCor.txt");
							//Inplace update of b1 and b2
							synchronized (b1) {
								glsStep1(geneZscoresSubset, geneInvCorMatrixSubsetMatrix, genePathwayZscoresSubset, b1, b2);
							}
							synchronized (b1NullGwas) {
								glsStep1(geneZscoresNullGwasNullBetasSubset, geneInvCorMatrixSubsetMatrix, genePathwayZscoresSubset, b1NullGwas, b2NullGwas);
							}

						}

						pb.step();
					} catch (Exception ex) {
						throw new RuntimeException(ex);
					}

				});

				DoubleMatrixDataset<String, String> enrichment = new DoubleMatrixDataset<>(genePathwayZscores.getHashColsCopy(), geneZscores.getHashColsCopy());
				DoubleMatrixDataset<String, String> enrichmentNull = new DoubleMatrixDataset<>(genePathwayZscores.getHashColsCopy(), geneZscoresNullGwasNullBetas.getHashColsCopy());

				final int numberOfPathways = genePathwayZscores.columns();
				final int numberTraits = geneZscores.columns();

				for (int traitI = 0; traitI < numberTraits; ++traitI) {
					final double b1Trait = b1.getElementQuick(traitI, 0);
					for (int pathwayI = 0; pathwayI < numberOfPathways; ++pathwayI) {

						double beta = b2.getElementQuick(traitI, pathwayI) / b1Trait;
						enrichment.setElementQuick(pathwayI, traitI, beta);

					}
				}

				final int numberTraitsNull = geneZscoresNullGwasNullBetas.columns();

				for (int traitI = 0; traitI < numberTraitsNull; ++traitI) {
					final double b1Trait = b1NullGwas.getElementQuick(traitI, 0);
					for (int pathwayI = 0; pathwayI < numberOfPathways; ++pathwayI) {

						double beta = b2NullGwas.getElementQuick(traitI, pathwayI) / b1Trait;
						enrichmentNull.setElementQuick(pathwayI, traitI, beta);

					}
				}

				enrichment.save(outputBasePath + "_" + pathwayDatabase.getName() + "_Enrichment" + (hlaGenesToExclude == null ? "_betas" : "_betasExHla") + ".txt");
				enrichmentNull.save(outputBasePath + "_" + pathwayDatabase.getName() + "_EnrichmentNull" + (hlaGenesToExclude == null ? "_betas" : "_betasExHla") + ".txt");

				DoubleMatrix2D enrichmentMatrix = enrichment.getMatrix();
				DoubleMatrix2D enrichmentNullMatrix = enrichmentNull.getMatrix();

				final int numberOfPhenotypes = enrichmentMatrix.columns();
				final int numberOfNullGwasPhenotypes = enrichmentNullMatrix.columns();
				final double numberOfNullGwasPhenotypesMin1Double = enrichmentNullMatrix.columns() - 1;

				LOGGER.debug("numberOfNullGwasPhenotypes: " + numberOfNullGwasPhenotypes);

				List<String> pathwayNames;
				if (LOGGER.isDebugEnabled()) {
					pathwayNames = enrichmentNull.getRowObjects();
				} else {
					pathwayNames = Collections.emptyList();
				}

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

				enrichment.save(outputBasePath + "_" + pathwayDatabase.getName() + "_Enrichment" + (hlaGenesToExclude == null ? "_zscore" : "_zscoreExHla") + ".txt");

				enrichmentZscores.put(pathwayDatabase, enrichment);

				LOGGER.debug("Completed " + pathwayDatabase.getName() + " enrichment");

				//This extra step is intentional
				pb.step();

			}
		}

		return enrichmentZscores;

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

	private static Map<String, ArrayList<Gene>> createChrArmGeneMapping(List<Gene> genes, Set<String> includedGenes) {
		Map<String, ArrayList<Gene>> chrArmToGeneMapping = new HashMap<>(25);
		for (Gene gene : genes) {

			if (includedGenes.contains(gene.getGene())) {

				String chrArm = gene.getChrAndArm();

				ArrayList<Gene> armGenes = chrArmToGeneMapping.get(chrArm);
				if (armGenes == null) {
					armGenes = new ArrayList<>();
					chrArmToGeneMapping.put(chrArm, armGenes);
				}

				armGenes.add(gene);

			}

		}
		return chrArmToGeneMapping;
	}

	private static LinkedHashSet<String> findUncorrelatedGenes(DoubleMatrixDataset<String, String> geneZscoresNullGwas, HashSet<String> genesWithPvalue, List<Gene> genes, double maxCorrelationBetweenGenes) {

		final Map<String, ArrayList<Gene>> chrArmToGeneMapping = createChrArmGeneMapping(genes, genesWithPvalue);

		final LinkedHashSet<String> includedUncorrelatedGenes = new LinkedHashSet<>();

		chrArmToGeneMapping.keySet().parallelStream().forEach((String chrArm) -> {

			final ArrayList<Gene> armGenes = chrArmToGeneMapping.get(chrArm);
			final ArrayList<String> armGenesIds = new ArrayList<>(armGenes.size());

			final HashSet<String> includedUncorrelatedGenesArm = new HashSet<>();

			for (Gene armGene : armGenes) {
				armGenesIds.add(armGene.getGene());
			}

			final DoubleMatrixDataset<String, String> geneZscoresNullGwasArm = geneZscoresNullGwas.viewRowSelection(armGenesIds);

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
					if (Math.abs(correlationOfCorrelations.getElementQuick(r, c)) >= maxCorrelationBetweenGenes && includedUncorrelatedGenesArm.contains(variantNames.get(c))) {
						continue rows;
					}
				}
				includedUncorrelatedGenesArm.add(variantNames.get(r));
			}

			synchronized (includedUncorrelatedGenes) {
				includedUncorrelatedGenes.addAll(includedUncorrelatedGenesArm);
			}

		});

		return includedUncorrelatedGenes;
	}

	private static DoubleMatrixDataset<String, String> createLocalGeneCorrelation(final DoubleMatrixDataset<String, String> geneZscoresNullGwasCorrelationSubset, final ArrayList<Gene> genes, final int correlationWindow) {

		final DoubleMatrixDataset<String, String> correlations = new DoubleMatrixDataset<>(geneZscoresNullGwasCorrelationSubset.getHashRows(), geneZscoresNullGwasCorrelationSubset.getHashRows());
		final DoubleMatrix2D correlationMatrix = correlations.getMatrix();
		final int geneCount = geneZscoresNullGwasCorrelationSubset.rows();
		final int nullGwasCount = geneZscoresNullGwasCorrelationSubset.columns();
		DoubleMatrix2D geneZscoresNullGwasCorrelationSubsetMatrix = geneZscoresNullGwasCorrelationSubset.getMatrix();

		final SimpleRegression regression = new SimpleRegression();

		for (int i = geneCount; --i >= 0;) {
			for (int j = i + 1; --j >= 0;) {
				regression.clear();

				if (i == j) {
					correlationMatrix.setQuick(i, j, 1);
				} else {

					//Genes should be in the same order as the matrix
					Gene geneI = genes.get(i);
					Gene geneJ = genes.get(j);

					//Only look at position because this is done per chromosome arm
					int geneIStart = geneI.getStart();
					int geneIStop = geneI.getStop();

					int geneJStart = geneJ.getStart();
					int geneJStop = geneJ.getStop();

					if (Math.abs(geneIStart - geneJStart) <= correlationWindow
							|| Math.abs(geneIStart - geneJStop) <= correlationWindow
							|| Math.abs(geneIStop - geneJStart) <= correlationWindow
							|| Math.abs(geneIStop - geneJStop) <= correlationWindow) {
						for (int n = 0; n < nullGwasCount; ++n) {
							regression.addData(geneZscoresNullGwasCorrelationSubsetMatrix.getQuick(i, n), geneZscoresNullGwasCorrelationSubsetMatrix.getQuick(j, n));
						}

						double x = regression.getR();

						correlationMatrix.setQuick(i, j, x);
						correlationMatrix.setQuick(j, i, x); // symmetric
					}

				}
			}
		}
		
		return correlations;

	}
}
