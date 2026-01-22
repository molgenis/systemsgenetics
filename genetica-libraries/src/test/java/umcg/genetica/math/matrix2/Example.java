package umcg.genetica.math.matrix2;

import java.io.IOException;
import java.text.DateFormat;
import java.text.SimpleDateFormat;
import java.util.*;

public class Example {


    private static void example() throws Exception {


        //Creat link to datg file only reading the row names and column names into memory
        final DoubleMatrixDatasetRowCompressedReader datgReader = new DoubleMatrixDatasetRowCompressedReader("pathToFile");

        //Print the meta data. Dataset name and row and column names are optional.
        final DateFormat DATE_TIME_FORMAT = new SimpleDateFormat("yyyy-MM-dd HH:mm:ss");

        System.out.printf("The dataset is named: '%s'%n", datgReader.getDatasetName());
        System.out.printf("The dataset was created on: %s%n", DATE_TIME_FORMAT.format(new Date(datgReader.getTimestampCreated() * 1000)));
        System.out.printf("Number of rows: %d content on rows: '%s'%n", datgReader.getNumberOfRows(), datgReader.getDataOnRows());
        System.out.printf("Number of columns: %d content on columns: '%s'%n", datgReader.getNumberOfColumns(), datgReader.getDataOnCols());

        /**
         *
         * Manual reading the datg file
         *
         */

        //Get column names
        final Set<String> columnIdentifiers = datgReader.getColumnIdentifiers();

        //Get row names with indices
        final Set<String> rowIdentifiers = datgReader.getRowIdentifiers();

        //Get map with row indices
        final Map<String, Integer> rowMap = datgReader.getRowMap();

        //Load a single row
        double[] rowData = datgReader.loadSingleRow(rowMap.get("Name of the row"));

        //Loop all rows:
        for(int r = 0; r < datgReader.getNumberOfRows(); r++) {
            rowData = datgReader.loadSingleRow(r);
        }

        /**
         *
         *
         * Manual writing a new datg file
         *
         */

        final ArrayList<String> columns = new ArrayList<String>(List.of(new String[]{"Column 1", "Column 2", "Column 3"}));

        //The 3 strings are optional meta data. They can be empty or use constructor with only path and columns: new DoubleMatrixDatasetRowCompressedWriter("FilePath", columns)
        final DoubleMatrixDatasetRowCompressedWriter datgWriter = new DoubleMatrixDatasetRowCompressedWriter("FilePath", columns, "dataset name", "data on row", "data on columns");

        //It is not needed to specify how many rows will be added just keep adding unique names with double[] with values
        double[] dataToWrite = new double[]{0.1, 0.2, 0.3};
        datgWriter.addRow("Row name", dataToWrite);

        //This is important otherwise fill can't be read.
        datgWriter.close();

        /**
         *
         *
         * Easy reading in DoubleMatrixDataset
         *
         */

        //Load full data into memory
        final DoubleMatrixDataset<String, String> dataAsDataset = datgReader.loadFullDataset();

        //Get row, column in dataset
        dataAsDataset.getElementQuick(1,2);

        //View subset of columns in order of columns array
        String[] columnsToView = new String[]{"Column 2", "Column 1"};
        DoubleMatrixDataset<String, String> dataColumnView = dataAsDataset.viewColSelection(columnsToView);

        //The dataColumnView is linked to the original data (IE not a copy)
        dataColumnView.getElementQuick(1,1);

        //Load subset of rows into memory in order of query
        String[] rowsToLoad = new String[]{"Row 2", "Row 1"};
        DoubleMatrixDataset<String, String> subsetOfData = datgReader.loadSubsetOfRows(rowsToLoad);

        //this can also be viewed in specific column order:
        DoubleMatrixDataset<String, String> subsetOfDataColomnView = subsetOfData.viewColSelection(columnsToView);

        subsetOfDataColomnView.getElementQuick(0,1);

        /**
         *
         *
         * Write dataset
         *
         */


        final ArrayList<String> rowNamesNewData = new ArrayList<String>(List.of(new String[]{"Row 1", "Row 2", "Row 3"}));
        final ArrayList<String> colNamesNewData = new ArrayList<String>(List.of(new String[]{"Column 1", "Column 2", "Column 3"}));

        //New doubleMatrixDataset with set row and column names
        final DoubleMatrixDataset<String, String> newData = new DoubleMatrixDataset<>(rowNamesNewData, colNamesNewData);

        newData.setElementQuick(0,0, Double.NaN);

        //Save as datg file, last 3 argument are optional and can also be empty strings.
        newData.saveBinary("File path", "dataset name", "data on row", "data on columns");

    }

}
