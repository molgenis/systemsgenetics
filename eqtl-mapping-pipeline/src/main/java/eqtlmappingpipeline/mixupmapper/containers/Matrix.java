/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */

package eqtlmappingpipeline.mixupmapper.containers;

/**
 *
 * @author harm-jan
 */
public class Matrix {
    private double[][] mat;
    private String[] rowNames;
    private String[] colNames;

    public Matrix(int numRows, int numCols){
        mat = new double[numRows][numCols];
    }

    public void addRow(int rowIndex, double[] row){
        mat[rowIndex] = row;
    }

    public void setRowNames(String[] r){
        rowNames = r;
    }

    public void setColNames(String[] c){
        colNames = c;
    }

    public double[][] getMatrix(){
        return mat;
    }

    public String[] getRowNames(){
        return rowNames;
    }

    public String[] getColNames(){
        return colNames;
    }

    public void setMatrix(double[][] matrix){
        mat = matrix;
    }

    public void scaleAndCenterOverRows(){
        // scale and center per row...
        if (mat.length > 0) {
            for (int row = 0; row < mat.length; row++) {
                double sum = 0.0;
                for (int col = 0; col < mat[row].length; col++) {
                    sum += mat[row][col];
                }

                double mean = sum / mat[row].length;
                double diff = 0.0;

                for (int col = 0; col < mat[row].length; col++) {
                    diff = diff + ((mat[row][col] - mean) * (mat[row][col] - mean));
                }

                double sd = Math.sqrt(diff / mat[row].length);

                for (int col = 0; col < mat[row].length; col++) {
                    mat[row][col] = (mat[row][col] - mean) / sd;
                }

            }
        }
    }

    public void scaleAndCenterOverCols() {
        if (mat.length > 0) {
            // scale and center per column
            for (int i = 0; i < mat[0].length; i++) {
                double[] col = new double[mat.length];
                double sum = 0.0;
                for (int j = 0; j < mat.length; j++) {
                    col[j] = mat[j][i];
                    sum += col[j];
                }

                double mean = sum / mat.length;
                double diff = 0.0;
                for (int j = 0; j < mat.length; j++) {
                    diff = diff + ((col[j] - mean) * (col[j] - mean));
                }

                double sd = Math.sqrt(diff / mat.length);

                // scale and center in one operation :D
                for (int j = 0; j < mat.length; j++) {
                    mat[j][i] = (col[j] - mean) / sd;
                }
            }

            
        } 
    }

    public double[] asArray() {
	double[] output = new double[mat.length*mat[0].length];
	int combos = 0;
	for(int i=0;i<mat.length; i++){
	    for(int j=0; j<mat[i].length; j++){
		output[combos] = mat[i][j];
		combos++;
	    }
	}
	return output;
    }
}
