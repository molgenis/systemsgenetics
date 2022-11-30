package nl.systemsgenetics.downstreamer.runners;

import cern.colt.list.tdouble.DoubleArrayList;
import cern.colt.list.tint.IntArrayList;
import cern.colt.matrix.tdouble.DoubleMatrix1D;
import me.tongfei.progressbar.ProgressBar;
import me.tongfei.progressbar.ProgressBarStyle;
import nl.systemsgenetics.downstreamer.io.IoUtils;
import nl.systemsgenetics.downstreamer.runners.options.DownstreamerOptionsDeprecated;
import nl.systemsgenetics.downstreamer.DownstreamerStep2Results;
import nl.systemsgenetics.downstreamer.pathway.PathwayDatabase;
import nl.systemsgenetics.downstreamer.pathway.PathwayEnrichments;
import nl.systemsgenetics.downstreamer.io.PathwayDatabaseEnrichmentExcelWriter;
import org.apache.log4j.Logger;
import umcg.genetica.math.matrix2.DoubleMatrixDataset;
import umcg.genetica.math.matrix2.DoubleMatrixDatasetFastSubsetLoader;
import umcg.genetica.math.stats.FisherExactTest;
import umcg.genetica.math.stats.MannWhitneyUTest2;

import java.util.*;
import java.util.concurrent.ConcurrentHashMap;
import java.util.stream.IntStream;

public class PathwayDatabaseEnrichments {
    private static final Logger LOGGER = Logger.getLogger(PathwayDatabaseEnrichments.class);

    private static final int minimalGeneCountInPathway = 10;

    public static void testPredictionPerformance(DownstreamerOptionsDeprecated options) throws Exception {
        PathwayDatabaseEnrichmentExcelWriter writer = new PathwayDatabaseEnrichmentExcelWriter(options);
        DownstreamerStep2Results step2Results = IoUtils.loadExistingStep2Results(options, true);
        List<PathwayDatabase> targetPathwayDatabases = options.getPathwayDatabases2();

        // Sheet for databases that have been enriched (most likely just GenePriortization)
        for (PathwayEnrichments queryEnrichment : step2Results.getPathwayEnrichments()) {
            // Map is structured as follows: GWAS trait > Pathway database > pathway results
            Map<String, Map<String, List<PathwayDatabaseEnrichmentRecord>>> results = testPredictionPerformanceForGwasTrait(queryEnrichment, targetPathwayDatabases);
            writer.writeResultsExcel(results, queryEnrichment.getPathwayDatabase().getName());
        }

        // Sheet for gene P-values
        PathwayEnrichments genePvalues = PathwayEnrichments.createPathwayEnrichmentsFromGenePvalues(options, step2Results.getGenePvalues());

		Map<String, Map<String, List<PathwayDatabaseEnrichmentRecord>>> results = testPredictionPerformanceForGwasTrait(genePvalues, targetPathwayDatabases);
        writer.writeResultsExcel(results, genePvalues.getPathwayDatabase().getName());
		
		// Sheet for closest gene
		PathwayEnrichments closestGene = PathwayEnrichments.createPathwayEnrichmentsFromClosestGene(options, step2Results.getGenePvalues());
		
        results = testPredictionPerformanceForGwasTrait(closestGene, targetPathwayDatabases);
        writer.writeResultsExcel(results, closestGene.getPathwayDatabase().getName());

    }

    private static Map<String, Map<String, List<PathwayDatabaseEnrichmentRecord>>> testPredictionPerformanceForGwasTrait(PathwayEnrichments query, List<PathwayDatabase> targetPathways) throws Exception {

        // Determine the genes in the query set
        final Set<String> queryGenes = new HashSet<>(query.getEnrichmentZscores().getRowObjects());
        final List<String> gwasTraits = query.getEnrichmentZscores().getColObjects();
        final Map<String, Map<String, List<PathwayDatabaseEnrichmentRecord>>> results = new ConcurrentHashMap<>();

        for (PathwayDatabase curTarget : targetPathways) {

            LOGGER.info("Determining overlapping genes");

            // Determine genes in the pathway
            DoubleMatrixDatasetFastSubsetLoader pathwayMatrixLoader = new DoubleMatrixDatasetFastSubsetLoader(curTarget.getLocation());
            Set<String> pathwayGenes = pathwayMatrixLoader.getOriginalRowMap();

            // Determine the overlapping genes
            Set<String> overlappingGenes = new HashSet<>(queryGenes);
            overlappingGenes.retainAll(pathwayGenes);

            // Load the pathway matrix
            LOGGER.info("Loading " + curTarget.getName());

            final DoubleMatrixDataset<String, String> pathwayMatrix;
            final DoubleMatrixDataset<String, String> pathwayMatrixUnfiltered = pathwayMatrixLoader.loadSubsetOfRowsBinaryDoubleData(overlappingGenes);

            LOGGER.info("Done loading, filtering pathways with fewer then " + minimalGeneCountInPathway + " genes");

            // Filter on pathways with at least 10 genes
            ArrayList<String> allPathways = pathwayMatrixUnfiltered.getColObjects();
            List<String> includedPathways = Collections.synchronizedList(new ArrayList<String>());

            IntStream.range(0, pathwayMatrixUnfiltered.columns()).parallel().forEach(pathwayI -> {
                if (pathwayMatrixUnfiltered.getCol(pathwayI).cardinality() >= minimalGeneCountInPathway) {
                    includedPathways.add(allPathways.get(pathwayI));
                }
            });
            pathwayMatrix = pathwayMatrixUnfiltered.viewColSelection(includedPathways);

            // Determine the bonf and FDR sig geneset
            final Map<String, Set<String>> bonfSigGenesPerGwas = new HashMap<>();
            final Map<String, Set<String>> fdrSigGenesPerGwas = new HashMap<>();

            for (String gwasTrait : gwasTraits) {
                Set<String> bonfSigGenes = getSignficiantGeneIds(query.getpValues(), query.getEnrichmentZscores(), gwasTrait, 0.05 / query.getpValues().rows());
                bonfSigGenes.retainAll(overlappingGenes);
                bonfSigGenesPerGwas.put(gwasTrait, bonfSigGenes);

                // Determine the FDR sig geneset
                Set<String> fdrSigGenes = getSignficiantGeneIds(query.getqValues(), query.getEnrichmentZscores(), gwasTrait, 0.05);
                fdrSigGenes.retainAll(overlappingGenes);
                fdrSigGenesPerGwas.put(gwasTrait, fdrSigGenes);
            }


            try (ProgressBar pb = new ProgressBar(query.getPathwayDatabase().getName() + "_" + curTarget.getName(), pathwayMatrix.columns(), ProgressBarStyle.ASCII)) {

                //TODO: threadpool is not well behaved here, it gobbles up all the threads, perhaps this is going wrong on the
                // options, also sometimes this randomly throws a nullpointer (maybe fixed?)
                //IntStream.range(0, pathwayMatrix.columns()).parallel().forEach(pathwayI -> {
                IntStream.range(0, pathwayMatrix.columns()).parallel().forEach(pathwayI -> {

                    // Pathway name
                    final String pathwayName = pathwayMatrix.getColObjects().get(pathwayI);

                    // Determine which genes are in the pathway
                    final DoubleMatrix1D pathwayAnnotation = pathwayMatrix.getCol(pathwayI);

                    final IntArrayList nonZeroIndices = new IntArrayList();
                    final DoubleArrayList nonZeroValues = new DoubleArrayList();
                    final Set<String> genesInPathway = new HashSet<>();
                    final Set<String> genesNotInPathway = new HashSet<>(pathwayMatrix.getRowObjects());
                    pathwayAnnotation.getNonZeros(nonZeroIndices, nonZeroValues);

                    for (int curIdx : nonZeroIndices.elements()) {
                        genesInPathway.add(pathwayMatrix.getRowObjects().get(curIdx));
                    }
                    genesNotInPathway.removeAll(genesInPathway);

                    // Subset Z-scores for mann whitney U-test
                    final DoubleMatrixDataset<String, String> zscoresPathwayGenes = query.getEnrichmentZscores().viewRowSelection(genesInPathway);
                    final DoubleMatrixDataset<String, String> zscoresNonPathwayGenes = query.getEnrichmentZscores().viewRowSelection(genesNotInPathway);

                    // Loop over the GWAS traits present for current run
                    for (String gwasTrait : gwasTraits) {

                        FisherExactGenesetResult bonfResult = genesetFisherExact(overlappingGenes, genesInPathway, bonfSigGenesPerGwas.get(gwasTrait));
                        FisherExactGenesetResult fdrResult = genesetFisherExact(overlappingGenes, genesInPathway, fdrSigGenesPerGwas.get(gwasTrait));

                        // Mann-whitney U-test
                        final MannWhitneyUTest2 uTest = new MannWhitneyUTest2();
                        uTest.setData(zscoresPathwayGenes.getCol(gwasTrait).toArray(), zscoresNonPathwayGenes.getCol(gwasTrait).toArray());

                        final double auc = uTest.getAuc();
                        final double pval = uTest.getP();

                        PathwayDatabaseEnrichmentRecord result = new PathwayDatabaseEnrichmentRecord(pathwayName, bonfResult, fdrResult, auc, pval);

                        // Add the result to the output
                        results.computeIfAbsent(gwasTrait, k -> new ConcurrentHashMap<>());
                        results.get(gwasTrait).computeIfAbsent(curTarget.getName(), k -> Collections.synchronizedList(new ArrayList<>()));

                        results.get(gwasTrait).get(curTarget.getName()).add(result);

                    }

                    pb.step();

                });
            }
        }

        return results;
    }

    private static Set<String> getSignficiantGeneIds(DoubleMatrixDataset<String, String> queryPvalues, DoubleMatrixDataset<String, String> queryZscores, String column, double threshold) {
        Set<String> significantGenes = new HashSet<>();
        int colIndex = queryPvalues.getColObjects().indexOf(column);

        for (int r = 0; r < queryPvalues.rows(); r++) {
            // Filter on genes that have a positive z-score and are below the significance threshold
            if (queryPvalues.getElementQuick(r, colIndex) < threshold && queryZscores.getElementQuick(r, colIndex) > 0) {
                significantGenes.add(queryPvalues.getRowObjects().get(r));
            }
        }

        return significantGenes;
    }

    private static FisherExactGenesetResult genesetFisherExact(Set<String> allGenes, Set<String> genesInPathway, Set<String> significantGenes) {

        // Determine which genes are overlapping for the contingency table
        Set<String> signifInPathway = new HashSet<>(genesInPathway);
        signifInPathway.retainAll(significantGenes);

        Set<String> notSignifInPathway = new HashSet<>(genesInPathway);
        notSignifInPathway.removeAll(signifInPathway);

        Set<String> signifNotPathway = new HashSet<>(allGenes);
        signifNotPathway.retainAll(significantGenes);
        signifNotPathway.removeAll(genesInPathway);

        Set<String> notSignifNotPathway = new HashSet<>(allGenes);
        notSignifNotPathway.removeAll(significantGenes);
        notSignifNotPathway.removeAll(genesInPathway);

        // Determine the contingency table
        int inPathwaySig = signifInPathway.size();
        int inPathwayNotSig = notSignifInPathway.size();
        int notPathwaySig = signifNotPathway.size();
        int notPathwayNotSig = notSignifNotPathway.size();

        // Fisher exact test
        FisherExactTest ft = new FisherExactTest();

        // First do this before single sided has a value;
        ft.getFisherPValue(inPathwaySig, inPathwayNotSig, notPathwaySig, notPathwayNotSig);
        double pval = ft.getFisherRightTail();
        double or = (double) (inPathwaySig * notPathwayNotSig) / (double) (inPathwayNotSig * notPathwaySig);

        return new FisherExactGenesetResult(signifInPathway, genesInPathway, pval, or);

    }

    /**
     * Store the results of a fisher exact test on a geneset, while preserving the overlapping gene IDs
     */
    public static class FisherExactGenesetResult {

        private final Set<String> significantOverlappingGenes;
        private final Set<String> genesInPathway;
        private final double pValue;
        private final double oddsRatio;

        public FisherExactGenesetResult(Set<String> significantOverlappingGenes, Set<String> genesInPathway, double pValue, double oddsRatio) {
            this.significantOverlappingGenes = significantOverlappingGenes;
            this.genesInPathway = genesInPathway;
            this.pValue = pValue;
            this.oddsRatio = oddsRatio;
        }

        public Set<String> getSignificantOverlappingGenes() {
            return significantOverlappingGenes;
        }

        public Set<String> getGenesInPathway() {
            return genesInPathway;
        }

        public double getpValue() {
            return pValue;
        }

        public double getOddsRatio() {
            return oddsRatio;
        }
    }

    /**
     * Store the results of a pathway enrichment for a single pathway and GWAS
     */
    public static class PathwayDatabaseEnrichmentRecord {
        private final String pathwayName;
        private final FisherExactGenesetResult bonfSigResult;
        private final FisherExactGenesetResult fdrSigResult;
        private final double auc;
        private final double aucPvalue;

        public PathwayDatabaseEnrichmentRecord(String pathwayName, FisherExactGenesetResult bonfSigResult, FisherExactGenesetResult fdrSigResult, double auc, double aucPvalue) {
            this.pathwayName = pathwayName;
            this.bonfSigResult = bonfSigResult;
            this.fdrSigResult = fdrSigResult;
            this.auc = auc;
            this.aucPvalue = aucPvalue;
        }

        public String getPathwayName() {
            return pathwayName;
        }

        public FisherExactGenesetResult getBonfSigResult() {
            return bonfSigResult;
        }

        public FisherExactGenesetResult getFdrSigResult() {
            return fdrSigResult;
        }

        public double getAuc() {
            return auc;
        }

        public double getAucPvalue() {
            return aucPvalue;
        }

        public double getBonfPvalue() {
            return bonfSigResult.getpValue();
        }
    }
}
