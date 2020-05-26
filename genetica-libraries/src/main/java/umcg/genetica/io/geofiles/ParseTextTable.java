/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package umcg.genetica.io.geofiles;

import cern.colt.matrix.tdouble.DoubleMatrix2D;
import cern.colt.matrix.tdouble.impl.DenseDoubleMatrix2D;
import cern.colt.matrix.tdouble.impl.DenseLargeDoubleMatrix2D;
import java.io.IOException;
import java.util.HashSet;
import java.util.LinkedHashMap;
import java.util.LinkedHashSet;
import java.util.logging.Level;
import java.util.logging.Logger;
import java.util.regex.Pattern;
import org.apache.commons.lang.math.NumberUtils;
import umcg.genetica.io.text.TextFile;
import umcg.genetica.math.matrix2.DoubleMatrixDataset;

/**
 *
 * @author MarcJan
 */
public class ParseTextTable {

    private static Pattern SPLIT_ON_TAB = Pattern.compile("\\t");
    protected static final String ENCODING = "ISO-8859-1";
    static final Logger LOGGER = Logger.getLogger(DoubleMatrixDataset.class.getName());

    /**
     * Geo text tables for HT12 v3 / v4 data
     *
     * @param fileInput
     * @return all data out of the TCGA files Methylated / un-methylated and
     * beta-values
     */
    public static DoubleMatrixDataset<String, String> parseGeoTables(String fileInput, boolean debug) throws IOException {
        System.out.println("\tNow parsing file: "+fileInput);
        LinkedHashSet<Integer> desiredColPos = new LinkedHashSet<Integer>();
        
        int columnOffset = 1;
        int rowOffset = 1;
        TextFile in = new TextFile(fileInput, TextFile.R);
        String str = in.readLine(); // header
        while (str.startsWith("#") || str.startsWith("\"#") || str.matches("^\\s*$") || str.equals("") || str.startsWith("This is our raw data.") || str.matches("^GSM[0-9]+.*")
                || str.startsWith("Illumina Inc. GenomeStudio") || str.startsWith("Array Content =") || str.startsWith("Normalization =") || str.startsWith("log")) {
            str = in.readLine();
            rowOffset++;
        }

        String[] headerData = SPLIT_ON_TAB.split(str);
        String str2 = in.readLine();
        String[] nextRowData = SPLIT_ON_TAB.split(str2);

        for (int s = 0; s < headerData.length; s++) {
            if ((headerData[s].toLowerCase().contains("probe") && headerData[s].toLowerCase().contains("id")) || (headerData[s].toLowerCase().contains("ref") && headerData[s].toLowerCase().contains("id")) || (headerData[s].toLowerCase().contains("array") && headerData[s].toLowerCase().contains("address"))) {
                columnOffset = (s + 1);
                break;
            }
        }
        int tmpCols = (headerData.length - columnOffset);

        LinkedHashMap<String, Integer> colMap = new LinkedHashMap<String, Integer>((int) Math.ceil(tmpCols / 0.75));

        int storedCols = 0;
        
        for (int s = 0; s < tmpCols; s++) {
            String colName = headerData[s + columnOffset];
            if (!colMap.containsKey(colName) && !colName.equals("") && !colName.equalsIgnoreCase("Target ID") && !colName.equalsIgnoreCase("TargetID") && !colName.equalsIgnoreCase("Probe ID") && !colName.toLowerCase().contains(".p=") && !colName.toLowerCase().contains("pval") && !colName.toLowerCase().contains("detection") && !colName.toLowerCase().contains("p-val") && !colName.toLowerCase().contains("array") && !colName.toLowerCase().contains("bead")) {
                if((nextRowData.length > (s+columnOffset))){
                    if((isNumeric(nextRowData[ s + columnOffset]))) {
                        colMap.put(colName, storedCols);
                        desiredColPos.add(s+columnOffset);
                        storedCols++;
                    } else if(debug) {
                        System.out.println("In non-numeric, entry: "+nextRowData[s + columnOffset]);
                        System.out.println("###############################");
                    }
                }
            } else if(colMap.containsKey(colName)) {
                LOGGER.warning("Duplicated column name:" + colName + "! In file: " + fileInput);
                throw new IOException("Problem with parsing file");
            } 
            else if (debug) {
                System.out.println("Empthy colname:"+colName.equals(""));
                System.out.println("Colname contains \"target id\":"+colName.equalsIgnoreCase("Target ID"));
                System.out.println("Colname contains \"probe id\":"+colName.equalsIgnoreCase("Probe ID"));
                System.out.println("Colname contains \".p\":"+colName.toLowerCase().contains(".p="));
                System.out.println("Colname contains \"pval\":"+colName.toLowerCase().contains("pval"));
                System.out.println("Colname contains \"p-val\":"+colName.toLowerCase().contains("p-val"));
                System.out.println("Colname contains \"detection\":"+colName.toLowerCase().contains("detection"));
                System.out.println("Colname contains \"array\":"+colName.toLowerCase().contains("array"));
                System.out.println("Colname contains \"bead\":"+colName.toLowerCase().contains("bead"));
                System.out.println("###############################");
            }
        }
        if (colMap.size() == 0) {
            if(debug){
                System.out.println("#Nothing added for this file: "+fileInput+". First two rows:");
                System.out.println("Parsing of values oke? "+ isNumeric(nextRowData[columnOffset]));
                System.out.println(str);
                System.out.println(str2);
                System.out.println(fileInput);
            }
            return null;
        }
        int tmpRows = 1;

        while (in.readLine() != null) {
            tmpRows++;
        }

        in.close();

        double[][] initialMatrix = new double[tmpRows][storedCols];

        in.open();
        String headerRow = null;
        for (int i = 0; i < rowOffset; ++i) {
            headerRow = in.readLine(); // read header
        }
        int row = 0;

        LinkedHashMap<String, Integer> rowMap = new LinkedHashMap<String, Integer>((int) Math.ceil(tmpRows / 0.75));

        boolean correctData = true;
        while ((str = in.readLine()) != null) {
            String[] data = SPLIT_ON_TAB.split(str);
            if(data.length == headerData.length){
               if (!rowMap.containsKey(data[columnOffset - 1])) {
                    rowMap.put(data[columnOffset - 1], row);
                    int columnToPut = 0;
                    for (int s : desiredColPos) {
                        double d;
                        try {
                            d = Double.parseDouble(data[s]);
                        } catch (NumberFormatException e) {
                            correctData = false;
                            d = Double.NaN;
                        }
                        initialMatrix[row][columnToPut] = d;
                        columnToPut++;
                    }
                    row++;
                } else {
                    LOGGER.warning("Duplicated row name: "+data[columnOffset - 1]);
                    System.out.println(str);
                    throw new IOException("Problem in reading file.");
                }
            }
        }
        if (!correctData) {
            LOGGER.warning("Your data contains NaN/unparseable values!");
        }
        in.close();

        DoubleMatrixDataset<String, String> dataset;
        
        DoubleMatrix2D mat;
        
        if ((tmpRows * tmpCols) < (Integer.MAX_VALUE - 2)) {
            mat = new DenseDoubleMatrix2D(initialMatrix);
        } else {
            mat = new DenseLargeDoubleMatrix2D(tmpRows, tmpCols);
            mat.assign(initialMatrix);            
        }
        
        dataset = new DoubleMatrixDataset<String, String>(mat, rowMap, colMap);
        
        LOGGER.log(Level.INFO, "''{0}'' has been loaded, nrRows: {1} nrCols: {2}", new Object[]{fileInput, dataset.rows(), dataset.columns()});
        return dataset;
    }

    public static boolean isNumeric(String str) {
        if (!NumberUtils.isNumber(str)) {
            return false;
        } else {
            return true;
        }
    }
}
