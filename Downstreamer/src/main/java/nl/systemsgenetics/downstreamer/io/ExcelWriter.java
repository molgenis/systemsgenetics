/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package nl.systemsgenetics.downstreamer.io;

import nl.systemsgenetics.downstreamer.DownstreamerStep3Results;
import nl.systemsgenetics.downstreamer.DownstreamerStep2Results;
import nl.systemsgenetics.downstreamer.Downstreamer;
import nl.systemsgenetics.downstreamer.DownstreamerOptions;
import cern.colt.matrix.tdouble.DoubleMatrix1D;

import java.io.File;
import java.io.FileOutputStream;
import java.io.IOException;
import java.util.*;
import nl.systemsgenetics.downstreamer.containers.GwasLocus;

import nl.systemsgenetics.downstreamer.gene.Gene;
import nl.systemsgenetics.downstreamer.gene.IndexedDouble;
import nl.systemsgenetics.downstreamer.pathway.PathwayAnnotations;
import nl.systemsgenetics.downstreamer.pathway.PathwayDatabase;
import nl.systemsgenetics.downstreamer.pathway.PathwayEnrichments;
import nl.systemsgenetics.downstreamer.runners.DownstreamerUtilities;
import static nl.systemsgenetics.downstreamer.runners.DownstreamerUtilities.getDistanceGeneToTopCisSnpPerTrait;
import org.apache.log4j.Logger;
import org.apache.poi.common.usermodel.HyperlinkType;
import org.apache.poi.ss.SpreadsheetVersion;
import org.apache.poi.ss.usermodel.*;
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
 * @author patri
 */
public class ExcelWriter {

	private static final Logger LOGGER = Logger.getLogger(ExcelWriter.class);

	private ExcelStyles styles;

	private final String outputBasePath;
	private final boolean hlaExcluded;
	private final List<String> traits;
	private final DownstreamerOptions options;

	public ExcelWriter(List<String> traits, DownstreamerOptions options) {
		this.outputBasePath = options.getOutputBasePath();
		this.hlaExcluded = options.isExcludeHla();
		this.traits = traits;
		this.options = options;
	}

	/**
	 * @throws java.io.FileNotFoundException
	 */
	public void saveStep2Excel(DownstreamerStep2Results results) throws Exception {

		DoubleMatrixDataset<String, String> genePvalues = results.getGenePvalues();
		List<PathwayEnrichments> pathwayEnrichments = results.getPathwayEnrichments();
		System.setProperty("java.awt.headless", "true");

		// Each trait gets its own sheet
		for (String trait : traits) {

			Workbook enrichmentWorkbook = new XSSFWorkbook();
			styles = new ExcelStyles(enrichmentWorkbook);
			CreationHelper createHelper = enrichmentWorkbook.getCreationHelper();

			// Overview sheet
			populateOverviewSheet(enrichmentWorkbook, trait, pathwayEnrichments, createHelper);

			// Sheet for gene pvalues
			//populateGenePvalueSheet(enrichmentWorkbook, trait, genePvalues);
			// Sheet for each pathway database
			for (PathwayEnrichments pathwayEnrichment : pathwayEnrichments) {
				populatePathwaySheet(enrichmentWorkbook, pathwayEnrichment, trait, createHelper, genePvalues);
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

	public void saveStep3Excel(DownstreamerStep2Results step2, DownstreamerStep3Results step3) throws IOException {

		Map<String, List<GwasLocus>> lociPerTrait = step3.getLoci();
		System.setProperty("java.awt.headless", "true");

		String pathwayDatabaseToScore = options.getPathwayDatabasesToAnnotateWithGwas().get(0);

		PathwayEnrichments pathwayEnrichments = null;
		for (PathwayEnrichments pathway : step2.getPathwayEnrichments()) {
			if (pathway.getPathwayDatabase().getName().equals(pathwayDatabaseToScore)) {
				pathwayEnrichments = pathway;
			}
		}

		if (pathwayEnrichments == null) {
			throw new IllegalArgumentException("Provided pathway database name not found. Please check if the pathway database names are correct in -pd and --annotDb");
		}

		// Each trait gets its own sheet
		for (String trait : traits) {

			Workbook enrichmentWorkbook = new XSSFWorkbook();
			styles = new ExcelStyles(enrichmentWorkbook);
			//CreationHelper createHelper = enrichmentWorkbook.getCreationHelper();

			populateCisPrioSheet(enrichmentWorkbook, trait, lociPerTrait.get(trait), pathwayEnrichments);
			// Save the file
			File excelFile = new File(outputBasePath + "_cisPrio" + (traits.size() > 1 ? "_" + trait : "") + (hlaExcluded ? "_exHla.xlsx" : ".xlsx"));
			int nr = 1;
			while (excelFile.exists()) {
				excelFile = new File(outputBasePath + "_cisPrio" + (traits.size() > 1 ? "_" + trait : "") + (hlaExcluded ? "_exHla" : "") + "_" + nr + ".xlsx");
				nr++;
			}
			enrichmentWorkbook.write(new FileOutputStream(excelFile));
		}
	}

	public void saveGenePvalueExcel(DoubleMatrixDataset<String, String> genePvalues) throws IOException {

		System.setProperty("java.awt.headless", "true");
		Workbook enrichmentWorkbook = new XSSFWorkbook();
		styles = new ExcelStyles(enrichmentWorkbook);

		populateGenePvalueSheet(enrichmentWorkbook, genePvalues);

		// Save the file
		File excelFile = new File(outputBasePath + "_genePvalues" + (hlaExcluded ? "_exHla.xlsx" : ".xlsx"));
		int nr = 1;
		while (excelFile.exists()) {
			excelFile = new File(outputBasePath + "_genePvalues" + (hlaExcluded ? "_exHla" : "") + "_" + nr + ".xlsx");
			nr++;
		}
		enrichmentWorkbook.write(new FileOutputStream(excelFile));

	}

	public void savePathwayLoadings(DownstreamerStep2Results step2Results) throws Exception {

		DoubleMatrixDataset<String, String> genePvalues = step2Results.getGenePvalues();
		double bonfSigLevel = 0.05 / genePvalues.rows();
		//TODO: TEMPORARY FOR ALS GWAS
		//LOGGER.warn("Signficiance level for GENES hardcoded to 1e-5");
		//bonfSigLevel = 1e-5;

		for (String trait : traits) {

			System.setProperty("java.awt.headless", "true");
			Workbook enrichmentWorkbook = new XSSFWorkbook();
			styles = new ExcelStyles(enrichmentWorkbook);

			// Determine significant genes
			Set<String> significantGenes = new HashSet<>();

			for (int i = 0; i < genePvalues.rows(); i++) {
				if (genePvalues.getCol(trait).get(i) < bonfSigLevel) {
					significantGenes.add(genePvalues.getRowObjects().get(i));
				}
			}

			// Determine signif pathways and subset matrix
			for (PathwayEnrichments enrichment : step2Results.getPathwayEnrichments()) {

				// Determine significant pathways
				DoubleMatrixDataset<String, String> curZscores = enrichment.getEnrichmentZscores();
				double bonfSigZscore = -ZScores.pToZ(0.05 / enrichment.getNumberOfPathways());

				Set<String> significantPathways = new HashSet<>();

				for (int i = 0; i < curZscores.rows(); i++) {
					if (curZscores.getCol(trait).get(i) > bonfSigZscore) {
						significantPathways.add(curZscores.getRowObjects().get(i));
					}
				}

				// Load the pathway and subset for signif genes
				PathwayDatabase curDb = enrichment.getPathwayDatabase();

				// Intersect all available genes for that pathway
				List<String> availGenes = DoubleMatrixDataset.readDoubleTextDataRowNames(curDb.getLocation() + ".rows.txt", '\t');
				significantGenes.retainAll(availGenes);

				// Make the subset
				DoubleMatrixDataset<String, String> curPathwaySubset = DoubleMatrixDataset.loadSubsetOfRowsBinaryDoubleData(curDb.getLocation(), significantGenes);
				curPathwaySubset = curPathwaySubset.viewColSelection(significantPathways);

				populatePathwayLoadingSheet(enrichmentWorkbook, enrichment, curPathwaySubset);
			}

			// Save the file
			File excelFile = new File(outputBasePath + "_pathwayLoadings" + (traits.size() > 1 ? "_" + trait : "") + (hlaExcluded ? "_exHla.xlsx" : ".xlsx"));
			int nr = 1;
			while (excelFile.exists()) {
				excelFile = new File(outputBasePath + "_pathwayLoadings" + (traits.size() > 1 ? "_" + trait : "") + (hlaExcluded ? "_exHla" : "") + "_" + nr + ".xlsx");
				nr++;
			}
			enrichmentWorkbook.write(new FileOutputStream(excelFile));
		}
	}

	private void populateCisPrioSheet(Workbook enrichmentWorkbook, String trait, List<GwasLocus> loci, PathwayEnrichments databaseForScore) throws IOException {

		int numberOfCols = 10;
		int numberOfRows = 0;
		
		LinkedHashMap<String, Integer> geneHashRow = databaseForScore.getEnrichmentZscores().getHashRows();
		
		for (GwasLocus curLocus : loci) {
			
			if(curLocus.getOverlappingGenes().isEmpty()){
				numberOfRows++;
			} else {
				for(Gene gene : curLocus.getOverlappingGenes()){
					if(geneHashRow.containsKey(gene.getGene())){
						numberOfRows++;
					}
				}
			}

		}
		
		XSSFSheet locusOverview = (XSSFSheet) enrichmentWorkbook.createSheet("LocusOverview");
		XSSFTable table = locusOverview.createTable(new AreaReference(new CellReference(0, 0),
				new CellReference(numberOfRows, numberOfCols),
				SpreadsheetVersion.EXCEL2007));

		table.setName("LocusOverview");
		table.setDisplayName("LocusOverview");
		table.setStyleName("TableStyleLight9");
		table.getCTTable().getTableStyleInfo().setShowRowStripes(true);
		table.getCTTable().addNewAutoFilter();
		XSSFRow headerRow = locusOverview.createRow(0);

		int hc = 0;
		headerRow.createCell(hc++, CellType.STRING).setCellValue("Locus id");
		headerRow.createCell(hc++, CellType.STRING).setCellValue("Locus name");
		headerRow.createCell(hc++, CellType.STRING).setCellValue("Locus chr");
		headerRow.createCell(hc++, CellType.STRING).setCellValue("Locus start");
		headerRow.createCell(hc++, CellType.STRING).setCellValue("Locus end");
		headerRow.createCell(hc++, CellType.STRING).setCellValue("Index SNP");
		headerRow.createCell(hc++, CellType.STRING).setCellValue("Index SNP P-value");
		headerRow.createCell(hc++, CellType.STRING).setCellValue("Gene id");
		headerRow.createCell(hc++, CellType.STRING).setCellValue("Gene name");
		headerRow.createCell(hc++, CellType.STRING).setCellValue("Gene score");
		headerRow.createCell(hc++, CellType.STRING).setCellValue("Distance to index SNP");

		// Sort loci by gene pvalue, smallest first
		loci.sort(Comparator.comparingDouble(GwasLocus::getPvalue));

		int i = 0;
		int r = 1; //+1 for header
		for (GwasLocus curLocus : loci) {
			XSSFRow row = locusOverview.createRow(r);

			// locus id
			// is below
			// locus  name
			XSSFCell locusNameCell = row.createCell(1, CellType.STRING);
			locusNameCell.setCellValue(curLocus.getContig()+ ":" + curLocus.getStart() + "-" + curLocus.getEnd());

			// chromosome
			XSSFCell chrCell = row.createCell(2, CellType.NUMERIC);
			chrCell.setCellValue(curLocus.getContig());

			// start
			XSSFCell startCell = row.createCell(3, CellType.NUMERIC);
			startCell.setCellValue(curLocus.getStart());
			startCell.setCellStyle(styles.getGenomicPositionStyle());

			// end
			XSSFCell endCell = row.createCell(4, CellType.NUMERIC);
			endCell.setCellValue(curLocus.getEnd());
			endCell.setCellStyle(styles.getGenomicPositionStyle());

			// Topsnp id
			XSSFCell snpCell = row.createCell(5, CellType.STRING);
			snpCell.setCellValue(curLocus.getLeadVariant().getVariantId());

			// Topsnp pvalue
			XSSFCell snpPvalCell = row.createCell(6, CellType.NUMERIC);
			double pvalue = curLocus.getPvalue();
			snpPvalCell.setCellValue(pvalue);
			snpPvalCell.setCellStyle(pvalue < 0.001 ? styles.getSmallPvalueStyle() : styles.getLargePvalueStyle());

			// Makes sure loci without genes are reported in the file
			if (curLocus.getOverlappingGenes().size() >= 1) {
				// Sort genes on zscore
				List<IndexedDouble> zscores = new ArrayList<>(curLocus.getOverlappingGenes().size());

				int idx = 0;
				for (Gene curGene : curLocus.getOverlappingGenes()) {
					double zscore;
					if (geneHashRow.containsKey(curGene.getGene())) {
						zscore = databaseForScore.getEnrichmentZscores().getElement(curGene.getGene(), trait);
					} else {
						zscore = Double.NaN;
					}
					zscores.add(new IndexedDouble(zscore, idx));
					idx++;
				}
				zscores.sort(Collections.reverseOrder());

				// Keep track of minimal distance to SNP
				int minimalSnpDistance = Integer.MAX_VALUE;
				int rowWithMinimalSnpDistance = r;

				for (IndexedDouble zscore : zscores) {
					Gene curGene = curLocus.getOverlappingGenes().get(zscore.getIndex());

					if (geneHashRow.containsKey(curGene.getGene())) {
						// Locus id
						XSSFCell locusIdCell = row.createCell(0, CellType.NUMERIC);
						locusIdCell.setCellValue(i + 1);

						// Gene id
						XSSFCell geneIdCell = row.createCell(7, CellType.STRING);
						geneIdCell.setCellValue(curGene.getGene());

						// Gene name
						XSSFCell geneNameCell = row.createCell(8, CellType.STRING);
						geneNameCell.setCellValue(curGene.getGeneSymbol());

						// Prioritzation zscore
						XSSFCell scoreCell = row.createCell(9, CellType.NUMERIC);
						scoreCell.setCellValue(zscore.getValue());
						scoreCell.setCellStyle(styles.getZscoreStyle());

						// Distance to index SNP
						XSSFCell geneDistanceCell = row.createCell(10, CellType.NUMERIC);

						// If SNP overlaps the gene body / intron, set distance to zero
						int dist;
						if (curGene.overlaps(curLocus.getLeadVariant())) {
							dist = 0;
						} else {
							dist = Math.min(Math.abs(curLocus.getLeadVariant().getPos() - curGene.getStart()), Math.abs(curLocus.getLeadVariant().getPos() - curGene.getEnd()));
						}
						geneDistanceCell.setCellValue(dist);
						geneDistanceCell.setCellStyle(styles.getGenomicPositionStyle());

						if (dist < minimalSnpDistance) {
							minimalSnpDistance = dist;
							rowWithMinimalSnpDistance = r;
						}

						r++;
						row = locusOverview.createRow(r);
					}

					/*                    else {
                        XSSFCell scoreCell = row.createCell(9, CellType.STRING);
                        scoreCell.setCellValue("NA");
                        scoreCell.setCellStyle(styles.getRightAlignedText());
                    }*/
				}
				i++;

				// Make the cell with the minimal distance bold
				if (locusOverview.getRow(rowWithMinimalSnpDistance) != null) {
					locusOverview.getRow(rowWithMinimalSnpDistance).getCell(10).setCellStyle(styles.getBoldGenomicPositionStyle());
				}
			} else {
				// Locus id
				XSSFCell locusIdCell = row.createCell(0, CellType.NUMERIC);
				locusIdCell.setCellValue(i + 1);
				i++;
				r++;
				continue;
			}
		}

		for (int c = 0; c < numberOfCols; ++c) {
			locusOverview.autoSizeColumn(c);
			locusOverview.setColumnWidth(c, locusOverview.getColumnWidth(c) + 1500); //compensate for with auto filter and inaccuracies
			if (c > 1 && locusOverview.getColumnWidth(c) > 20000) {
				//max col width. Not for first column.
				locusOverview.setColumnWidth(c, 20000);
			}
		}

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
		cell.setCellValue("Number of sets");
		cell.setCellStyle(styles.getBoldStyle());

		for (PathwayEnrichments pathwayEnrichment : pathwayEnrichments) {
			row = overviewSheet.createRow(r++);
			cell = row.createCell(0, CellType.STRING);
			cell.setCellValue(pathwayEnrichment.getPathwayDatabase().getName());

			Hyperlink link = createHelper.createHyperlink(HyperlinkType.DOCUMENT);
			link.setAddress(pathwayEnrichment.getPathwayDatabase().getName() + "!A1");
			cell.setHyperlink(link);
			cell.setCellStyle(styles.getHlinkStyle());

			row.createCell(1, CellType.NUMERIC).setCellValue(pathwayEnrichment.getNumberOfPathways());
		}

		for (int c = 0; c < 2; ++c) {
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

	private void populatePathwaySheet(Workbook enrichmentWorkbook, PathwayEnrichments pathwayEnrichment, String trait, CreationHelper createHelper, DoubleMatrixDataset<String, String> genePvalues) throws Exception {

		PathwayDatabase pathwayDatabase = pathwayEnrichment.getPathwayDatabase();

		PathwayAnnotations pathwayAnnotations = new PathwayAnnotations(new File(pathwayDatabase.getLocation() + ".colAnnotations.txt"));
		int maxAnnotations = pathwayAnnotations.getMaxNumberOfAnnotations();

		DoubleMatrixDataset<String, String> databaseEnrichmentZscores = pathwayEnrichment.getEnrichmentZscores();
		DoubleMatrixDataset<String, String> databaseEnrichmentQvalues = pathwayEnrichment.getqValues();

		//final DoubleMatrixDataset<String, String> gwasSnpPvalues = DoubleMatrixDataset.loadDoubleBinaryData(options.getGwasZscoreMatrixPath());

		ArrayList<String> geneSets = databaseEnrichmentZscores.getRowObjects();
		double bonferroniCutoff = 0.05 / pathwayEnrichment.getNumberOfPathways();
		DoubleMatrix1D traitEnrichment = databaseEnrichmentZscores.getCol(trait);
		DoubleMatrix1D traitQvalue = databaseEnrichmentQvalues.getCol(trait);
		int[] order = DoubleMatrix1dOrder.sortIndexReverse(traitEnrichment);
		final boolean annotateWithGwasData = options.getPathwayDatabasesToAnnotateWithGwas().contains(pathwayDatabase.getName());
		int gwasAnnotations = 0;
		final int windowExtend = options.getCisWindowExtend();
		final String transLabel = "Trans (>" + (windowExtend >= 1000 ? ((windowExtend / 1000) + " k" ) : (windowExtend + " "))  + "b)";

		if (annotateWithGwasData) {
			gwasAnnotations = 4;
		}

		LOGGER.debug(pathwayDatabase.getName());
		LOGGER.debug(annotateWithGwasData);

		final HashMap<String, DownstreamerUtilities.NearestVariant> distanceGeneToTopCisSnp;
		//Map<String, Gene> genes = null;

		if (annotateWithGwasData) {
			HashMap<String, HashMap<String, DownstreamerUtilities.NearestVariant>> distanceGeneToTopCisSnpPerTrait = getDistanceGeneToTopCisSnpPerTrait(options);
			distanceGeneToTopCisSnp = distanceGeneToTopCisSnpPerTrait.get(trait);
			//genes = IoUtils.readGenesMap(options.getGeneInfoFile());
		} else {
			distanceGeneToTopCisSnp = null;
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

		// Header row
		XSSFRow headerRow = databaseSheet.createRow(0);
		int hc = 0;

		// Use colnames from .colAnnot file if available
		headerRow.createCell(hc++, CellType.STRING).setCellValue(pathwayAnnotations.getSetName() == null ? "Gene set" : pathwayAnnotations.getSetName());
		for (int i = 0; i < maxAnnotations; ++i) {
			headerRow.createCell(hc++, CellType.STRING).setCellValue(pathwayAnnotations.getAnnotationHeaders().get(i));
		}
		headerRow.createCell(hc++, CellType.STRING).setCellValue("Enrichment Z-score");
		headerRow.createCell(hc++, CellType.STRING).setCellValue("Enrichment P-value");
		headerRow.createCell(hc++, CellType.STRING).setCellValue("Enrichment Q-value");
		headerRow.createCell(hc++, CellType.STRING).setCellValue("Bonferroni significant");
		headerRow.createCell(hc++, CellType.STRING).setCellValue("FDR 5% significant");

		if (annotateWithGwasData) {
			headerRow.createCell(hc++, CellType.STRING).setCellValue("Distance to lead GWAS variant");
			headerRow.createCell(hc++, CellType.STRING).setCellValue("GWAS variant ID");
			headerRow.createCell(hc++, CellType.STRING).setCellValue("GWAS variant P-value");
			headerRow.createCell(hc++, CellType.STRING).setCellValue("GWAS gene P-value");
			
			
		}

		// Populate the rows, loop over each row
		for (int r = 0; r < databaseEnrichmentZscores.rows(); ++r) {
			XSSFRow row = databaseSheet.createRow(r + 1); //+1 for header
			String geneSet = geneSets.get(order[r]);
			row.createCell(0, CellType.STRING).setCellValue(geneSet);

			// Annotations from .colAnnotations file
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
								cell.setCellStyle(styles.getHlinkStyle());
							}

						} else {
							row.createCell(j + 1, CellType.STRING).setCellValue("");
						}

					}
				}
			}

			// Fixed fields from enrichment results
			double zscore = traitEnrichment.getQuick(order[r]);
			double qvalue = traitQvalue.getQuick(order[r]);

			// Zscore
			XSSFCell zscoreCell = row.createCell(1 + maxAnnotations, CellType.NUMERIC);
			zscoreCell.setCellValue(zscore);
			zscoreCell.setCellStyle(styles.getZscoreStyle());

			// Pvalue
			double pvalue = ZScores.zToP(zscore);
			XSSFCell pvalueCell = row.createCell(2 + maxAnnotations, CellType.NUMERIC);
			pvalueCell.setCellValue(pvalue);
			pvalueCell.setCellStyle(pvalue < 0.001 ? styles.getSmallPvalueStyle() : styles.getLargePvalueStyle());

			// FDR
			XSSFCell qvalueCell = row.createCell(3 + maxAnnotations, CellType.NUMERIC);
			qvalueCell.setCellValue(qvalue);
			qvalueCell.setCellStyle(qvalue > 0 && qvalue < 0.001 ? styles.getSmallPvalueStyle() : styles.getLargePvalueStyle());

			// Is bonerroni y/n
			XSSFCell bonferroniCell = row.createCell(4 + maxAnnotations, CellType.BOOLEAN);
			bonferroniCell.setCellValue(pvalue <= bonferroniCutoff);

			// Is FDR y/n
			XSSFCell fdrCell = row.createCell(5 + maxAnnotations, CellType.BOOLEAN);
			fdrCell.setCellValue(qvalue <= 0.05);

			// Optional GWAS fields
			if (annotateWithGwasData) {

				String geneId = databaseEnrichmentZscores.getRowObjects().get(order[r]);
				//Gene curGene = genes.get(geneId);

				// Determine the closest independent tophit
				// No overlap = -9
//				int closestDist = -8;
//				SummaryStatisticRecord closestVariant = null;
//				String variantId = "";
				if (distanceGeneToTopCisSnp.containsKey(geneId)) {
					DownstreamerUtilities.NearestVariant nearest = distanceGeneToTopCisSnp.get(geneId);
					int dist = nearest.getDistance();
					if (dist < 0) {
						XSSFCell snpDistCell = row.createCell(6 + maxAnnotations, CellType.STRING);
						
						snpDistCell.setCellValue(transLabel);
						snpDistCell.setCellStyle(styles.getRightAlignedText());
					} else {
						XSSFCell snpDistCell = row.createCell(6 + maxAnnotations, CellType.NUMERIC);
						snpDistCell.setCellValue(dist);
						snpDistCell.setCellStyle(styles.getGenomicPositionStyle());

						// SNP name
						XSSFCell snpNameCell = row.createCell(7 + maxAnnotations, CellType.STRING);
						snpNameCell.setCellValue(nearest.getNearestVariant().getVariantId());

						// SNP p-value
						double snpPvalue = nearest.getNearestVariant().getpValue();
//						if (gwasSnpPvalues.getRowObjects().contains(variantId)) {
//							snpPvalue = ZScores.zToP(gwasSnpPvalues.getElement(variantId, trait));
//						} else {
//							snpPvalue = -9;
//						}

						XSSFCell snpPvalueCell = row.createCell(8 + maxAnnotations, CellType.NUMERIC);
						snpPvalueCell.setCellValue(snpPvalue);
						snpPvalueCell.setCellStyle(snpPvalue < 0.001 ? styles.getSmallPvalueStyle() : styles.getLargePvalueStyle());

						// Gene p-value
						double genePvalue = genePvalues.getElement(geneId, trait);

						if (!Double.isNaN(genePvalue)) {
							XSSFCell genePCell = row.createCell(9 + maxAnnotations, CellType.NUMERIC);
							//genePvalue = ZScores.zToP(genePvalue);
							genePCell.setCellValue(genePvalue);
							genePCell.setCellStyle(genePvalue < 0.001 ? styles.getSmallPvalueStyle() : styles.getLargePvalueStyle());
						} else {
							XSSFCell genePCell = row.createCell(9 + maxAnnotations, CellType.STRING);
							genePCell.setCellValue("NA");
							genePCell.setCellStyle(styles.getRightAlignedText());
						}

					}
				} else {
					XSSFCell snpDistCell = row.createCell(6 + maxAnnotations, CellType.STRING);
					snpDistCell.setCellValue("Missing");
					snpDistCell.setCellStyle(styles.getRightAlignedText());
				}

			}
		}

		// Auto-scale columns in sheet
		for (int c = 0; c < hc; ++c) {
			databaseSheet.autoSizeColumn(c);
			databaseSheet.setColumnWidth(c, databaseSheet.getColumnWidth(c) + 1500); //compensate for with auto filter and inaccuracies
			if (c > 1 && databaseSheet.getColumnWidth(c) > 20000) {
				//max col width. Not for first column.
				databaseSheet.setColumnWidth(c, 20000);
			}
		}

	}

	private void populateGenePvalueSheet(Workbook enrichmentWorkbook, DoubleMatrixDataset<String, String> genePvalues) throws IOException {

		Map<String, Gene> geneInfo = IoUtils.readGenesMap(options.getGeneInfoFile());
		XSSFSheet genePSheet = (XSSFSheet) enrichmentWorkbook.createSheet("GenePvalues");
		XSSFTable table = genePSheet.createTable(new AreaReference(new CellReference(0, 0),
				new CellReference(genePvalues.rows(), 5 + genePvalues.columns()),
				SpreadsheetVersion.EXCEL2007));

		table.setName("GenePvalues_res");
		table.setDisplayName("GenePvalues");
		table.setStyleName("TableStyleLight9");
		table.getCTTable().getTableStyleInfo().setShowRowStripes(true);
		table.getCTTable().addNewAutoFilter();
		XSSFRow headerRow = genePSheet.createRow(0);

		int hc = 0;
		headerRow.createCell(hc++, CellType.STRING).setCellValue("Gene id");
		headerRow.createCell(hc++, CellType.STRING).setCellValue("Gene name");
		headerRow.createCell(hc++, CellType.STRING).setCellValue("Chromosome");
		headerRow.createCell(hc++, CellType.STRING).setCellValue("Start");
		headerRow.createCell(hc++, CellType.STRING).setCellValue("End");

		for (String trait : genePvalues.getColObjects()) {
			// p-value
			headerRow.createCell(hc++, CellType.STRING).setCellValue(trait);
		}

		for (int r = 0; r < genePvalues.rows(); ++r) {
			XSSFRow row = genePSheet.createRow(r + 1); //+1 for header
			// gene id
			String geneId = genePvalues.getRowObjects().get(r);
			XSSFCell geneIdCell = row.createCell(0, CellType.STRING);
			geneIdCell.setCellValue(geneId);

			// gene name
			XSSFCell geneNameCell = row.createCell(1, CellType.STRING);
			geneNameCell.setCellValue(geneInfo.get(geneId).getGeneSymbol());

			// chromosome
			XSSFCell chrCell = row.createCell(2, CellType.NUMERIC);
			chrCell.setCellValue(geneInfo.get(geneId).getChr());

			// start
			XSSFCell startCell = row.createCell(3, CellType.NUMERIC);
			startCell.setCellValue(geneInfo.get(geneId).getStart());
			startCell.setCellStyle(styles.getGenomicPositionStyle());

			// end
			XSSFCell endCell = row.createCell(4, CellType.NUMERIC);
			endCell.setCellValue(geneInfo.get(geneId).getEnd());
			startCell.setCellStyle(styles.getGenomicPositionStyle());

			int x = 1;
			for (String trait : genePvalues.getColObjects()) {
				// p-value
				double pvalue = genePvalues.getElement(r, genePvalues.getColIndex(trait));

				if (!Double.isNaN(pvalue)) {
					XSSFCell pvalueCell = row.createCell(4 + x, CellType.NUMERIC);
					pvalueCell.setCellValue(pvalue);
					pvalueCell.setCellStyle(pvalue < 0.001 ? styles.getSmallPvalueStyle() : styles.getLargePvalueStyle());
				} else {
					XSSFCell pvalueCell = row.createCell(4 + x, CellType.STRING);
					pvalueCell.setCellValue("NA");
					pvalueCell.setCellStyle(styles.getRightAlignedText());

				}

				x++;
			}
		}

		// Auto-scale columns in sheet
		for (int c = 0; c < hc; ++c) {
			genePSheet.autoSizeColumn(c);
			genePSheet.setColumnWidth(c, genePSheet.getColumnWidth(c) + 1500); //compensate for with auto filter and inaccuracies
			if (c > 1 && genePSheet.getColumnWidth(c) > 20000) {
				//max col width. Not for first column.
				genePSheet.setColumnWidth(c, 20000);
			}
		}
	}

	/**
	 * Populates a tab with the pathway loadings for each for the provided
	 * subset of genes and pathways
	 *
	 * @param enrichmentWorkbook
	 * @param pathwayEnrichment
	 * @param pathwayZscores
	 * @throws IOException
	 */
	private void populatePathwayLoadingSheet(Workbook enrichmentWorkbook, PathwayEnrichments pathwayEnrichment, DoubleMatrixDataset<String, String> pathwayZscores) throws IOException {

		PathwayDatabase pathwayDatabase = pathwayEnrichment.getPathwayDatabase();

		Map<String, Gene> geneInfo = IoUtils.readGenesMap(options.getGeneInfoFile());
		XSSFSheet zscoreSheet = (XSSFSheet) enrichmentWorkbook.createSheet(pathwayDatabase.getName());
		XSSFTable table = zscoreSheet.createTable(new AreaReference(new CellReference(0, 0),
				new CellReference(pathwayZscores.rows(), 5 + pathwayZscores.columns()),
				SpreadsheetVersion.EXCEL2007));

		table.setName(pathwayDatabase.getName() + "_res");
		table.setDisplayName(pathwayDatabase.getName());
		table.setStyleName("TableStyleLight9");
		table.getCTTable().getTableStyleInfo().setShowRowStripes(true);
		table.getCTTable().addNewAutoFilter();
		XSSFRow headerRow = zscoreSheet.createRow(0);

		int hc = 0;
		headerRow.createCell(hc++, CellType.STRING).setCellValue("Gene id");
		headerRow.createCell(hc++, CellType.STRING).setCellValue("Gene name");
		headerRow.createCell(hc++, CellType.STRING).setCellValue("Chromosome");
		headerRow.createCell(hc++, CellType.STRING).setCellValue("Start");
		headerRow.createCell(hc++, CellType.STRING).setCellValue("End");

		for (String pathway : pathwayZscores.getColObjects()) {
			// zscore
			headerRow.createCell(hc++, CellType.STRING).setCellValue(pathway);
		}

		for (int r = 0; r < pathwayZscores.rows(); ++r) {
			XSSFRow row = zscoreSheet.createRow(r + 1); //+1 for header
			// gene id
			String geneId = pathwayZscores.getRowObjects().get(r);
			XSSFCell geneIdCell = row.createCell(0, CellType.STRING);
			geneIdCell.setCellValue(geneId);

			// gene name
			XSSFCell geneNameCell = row.createCell(1, CellType.STRING);
			geneNameCell.setCellValue(geneInfo.get(geneId).getGeneSymbol());

			// chromosome
			XSSFCell chrCell = row.createCell(2, CellType.NUMERIC);
			chrCell.setCellValue(geneInfo.get(geneId).getChr());

			// start
			XSSFCell startCell = row.createCell(3, CellType.NUMERIC);
			startCell.setCellValue(geneInfo.get(geneId).getStart());
			startCell.setCellStyle(styles.getGenomicPositionStyle());

			// end
			XSSFCell endCell = row.createCell(4, CellType.NUMERIC);
			endCell.setCellValue(geneInfo.get(geneId).getEnd());
			startCell.setCellStyle(styles.getGenomicPositionStyle());

			int x = 1;
			for (String pathway : pathwayZscores.getColObjects()) {
				double zscore = pathwayZscores.getElement(r, pathwayZscores.getColIndex(pathway));

				if (!Double.isNaN(zscore)) {
					XSSFCell zscoreCell = row.createCell(4 + x, CellType.NUMERIC);
					zscoreCell.setCellValue(zscore);
					zscoreCell.setCellStyle(styles.getZscoreStyle());
				} else {
					XSSFCell zscoreCell = row.createCell(4 + x, CellType.STRING);
					zscoreCell.setCellValue("NA");
					zscoreCell.setCellStyle(styles.getRightAlignedText());
				}

				x++;
			}
		}

		// Auto-scale columns in sheet
		for (int c = 0; c < hc; ++c) {
			zscoreSheet.autoSizeColumn(c);
			zscoreSheet.setColumnWidth(c, zscoreSheet.getColumnWidth(c) + 1500); //compensate for with auto filter and inaccuracies
			if (c > 1 && zscoreSheet.getColumnWidth(c) > 20000) {
				//max col width. Not for first column.
				zscoreSheet.setColumnWidth(c, 20000);
			}
		}
	}

}
