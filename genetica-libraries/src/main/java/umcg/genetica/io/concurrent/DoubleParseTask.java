/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package umcg.genetica.io.concurrent;

import java.util.Set;
import java.util.concurrent.Callable;
import java.util.regex.Pattern;
import umcg.genetica.containers.Triple;
import umcg.genetica.text.Strings;

/**
 *
 * @author harmjan
 */
public class DoubleParseTask implements Callable<Triple<Integer, String, double[]>> {

    String stringData;
    int offset = 0;
    int row = 0;
    private final int[] columnIndex;
    private final Set hashRowsToInclude;
    private int nrColumnsToInclude;
    private Pattern separator;

    public DoubleParseTask(String d, int offset, int row, int[] columnIndex, Set hashProbesToInclude) {
        this.stringData = d;
        this.offset = offset;
        this.row = row;
        this.columnIndex = columnIndex; // converts column position in output to column position 
        this.hashRowsToInclude = hashProbesToInclude;
        if (columnIndex == null) {
            this.nrColumnsToInclude = -1;
        } else {
            this.nrColumnsToInclude = columnIndex.length;
        }
    }
    
    public DoubleParseTask(String d, int offset, int row, int[] columnIndex, Set hashProbesToInclude, Pattern separator) {
        this.stringData = d;
        this.offset = offset;
        this.row = row;
        this.columnIndex = columnIndex; // converts column position in output to column position 
        this.hashRowsToInclude = hashProbesToInclude;
        if (columnIndex == null) {
            this.nrColumnsToInclude = -1;
        } else {
            this.nrColumnsToInclude = columnIndex.length;
        }
        this.separator = separator;
    }

    @Override
    public Triple<Integer, String, double[]> call() throws Exception {
        if (stringData == null) {
            return new Triple<Integer, String, double[]>(-1, null, null);
        }

        if(separator == null){
            separator = Strings.tab;
        }
        
        String[] splitData = separator.split(stringData);
        String rowname = new String(splitData[0].getBytes());

        if (this.nrColumnsToInclude == -1) {
            this.nrColumnsToInclude = splitData.length - 1;
        }

//        if(columnIndex == null){
//            System.out.println("Column index == null");
//        } else {
//            System.out.println("Column index == set");
//        }
        
        if (hashRowsToInclude == null || hashRowsToInclude.contains(rowname)) {
            double[] output = new double[nrColumnsToInclude];
            for (int s = 0; s < nrColumnsToInclude; s++) {
                int columnPositionInFile = s + offset;
                if (columnIndex != null) {
                    columnPositionInFile = columnIndex[s];
                }

                try {
                    output[s] = Double.parseDouble(splitData[columnPositionInFile]);
                } catch (NumberFormatException e) {
                    System.err.println("ERROR! Value is not a double: " + splitData[columnPositionInFile] + "\trow: " + row + "\tcol:" + columnPositionInFile +"\toffset: "+offset);
                    output[s] = Double.NaN;
////                    e.printStackTrace();
//                    System.exit(0);
                }
            }
            this.stringData = null;
            return new Triple<Integer, String, double[]>(row, rowname, output);
        } else {
            this.stringData = null;
            return new Triple<Integer, String, double[]>(-1, null, null);
        }

    }
}
