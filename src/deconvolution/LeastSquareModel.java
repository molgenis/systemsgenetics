package deconvolution;

import java.util.ArrayList;
import java.util.List;

/** 
* TODO: Now doing everything 3x for the 3 genotypes, should be merged to 1 single function that takes as input which genotype it is
*/
public class LeastSquareModel {
	// save the expression values for each dosage separatly
	private List<Double> expressionValuesGt0 = new ArrayList<Double>();
	private List<Double> expressionValuesGt1 = new ArrayList<Double>();
	private List<Double> expressionValuesGt2 = new ArrayList<Double>();
	private List<List<Double>> cellCountsGt0 = new ArrayList<List<Double>>();
	private List<List<Double>> cellCountsGt1 = new ArrayList<List<Double>>();
	private List<List<Double>> cellCountsGt2 = new ArrayList<List<Double>>();
	private String qtlName;
	private List<String> celltypes = new ArrayList<String>();
	private Boolean noIntercept;
	LeastSquareModel(){};

	public void addExpressionValuesGt0(double expressionValue){
	    this.expressionValuesGt0.add(expressionValue);
	  }
	
	public void addExpressionValuesGt1(double expressionValue){
		this.expressionValuesGt1.add(expressionValue);
	  }
	
	public void addExpressionValuesGt2(double expressionValue){
		this.expressionValuesGt2.add(expressionValue);
	  }
	
	public double[] getExpessionValuesGt0() throws IllegalAccessException{
		if(this.expressionValuesGt0 == null){
			throw new IllegalAccessException("expressionValues not set for this model");
		}
		// need double[] for regression model
		 double[] expressionVector = new double[this.expressionValuesGt0.size()];
		 for (int i = 0; i < expressionVector.length; i++) {
		    expressionVector[i] = this.expressionValuesGt0.get(i);                // java 1.5+ style (outboxing)
		 }
	    return(expressionVector);
	  }
	
	public double[] getExpessionValuesGt1() throws IllegalAccessException{
		if(this.expressionValuesGt1 == null){
			throw new IllegalAccessException("expressionValues not set for this model");
		}
		// need double[] for regression model
		 double[] expressionVector = new double[this.expressionValuesGt1.size()];
		 for (int i = 0; i < expressionVector.length; i++) {
		    expressionVector[i] = this.expressionValuesGt1.get(i);                // java 1.5+ style (outboxing)
		 }
	    return(expressionVector);
	  }
	
	public double[] getExpessionValuesGt2() throws IllegalAccessException{
		if(this.expressionValuesGt2 == null){
			throw new IllegalAccessException("expressionValues not set for this model");
		}
		// need double[] for regression model
		 double[] expressionVector = new double[this.expressionValuesGt2.size()];
		 for (int i = 0; i < expressionVector.length; i++) {
		    expressionVector[i] = this.expressionValuesGt2.get(i);                // java 1.5+ style (outboxing)
		 }
	    return(expressionVector);
	  }
	
	public void addCellCountsGt0(int celltypeIndex, double cellcount){
	    this.cellCountsGt0.get(celltypeIndex).add(cellcount);
	  }
	
	public void addCellCountsGt1(int celltypeIndex, double cellcount){
		this.cellCountsGt1.get(celltypeIndex).add(cellcount);
	  }
	
	public void addCellCountsGt2(int celltypeIndex, double cellcount){
		this.cellCountsGt2.get(celltypeIndex).add(cellcount);
	  }
	
	public double[][] getCellCountsGt0() throws IllegalAccessException{
		if(this.cellCountsGt0 == null){
			throw new IllegalAccessException("expressionValues not set for this model");
		}
		// need double[][] for regression model
		 double[][] cellCountMatrix = new double[this.cellCountsGt0.size()][this.cellCountsGt0.get(0).size()];
		 for (int i = 0; i < cellCountMatrix.length; i++) {
			 for (int j = 0; j < cellCountMatrix[0].length; j++){
			 cellCountMatrix[i][j] = this.cellCountsGt0.get(i).get(j);                // java 1.5+ style (outboxing)
			 }
		}
	    return(cellCountMatrix);
	  }

	public double[][] getCellCountsGt1() throws IllegalAccessException{
		if(this.cellCountsGt1 == null){
			throw new IllegalAccessException("expressionValues not set for this model");
		}
		// need double[][] for regression model
		 double[][] cellCountMatrix = new double[this.cellCountsGt1.size()][this.cellCountsGt1.get(0).size()];
		 for (int i = 0; i < cellCountMatrix.length; i++) {
			 for (int j = 0; j < cellCountMatrix[0].length; j++){
			 cellCountMatrix[i][j] = this.cellCountsGt1.get(i).get(j);                // java 1.5+ style (outboxing)
			 }
		}
	    return(cellCountMatrix);
	  }
	

	public double[][] getCellCountsGt2() throws IllegalAccessException{
		if(this.cellCountsGt2 == null){
			throw new IllegalAccessException("expressionValues not set for this model");
		}
		// need double[][] for regression model
		 double[][] cellCountMatrix = new double[this.cellCountsGt2.size()][this.cellCountsGt2.get(0).size()];
		 for (int i = 0; i < cellCountMatrix.length; i++) {
			 for (int j = 0; j < cellCountMatrix[0].length; j++){
			 cellCountMatrix[i][j] = this.cellCountsGt2.get(i).get(j);                // java 1.5+ style (outboxing)
			 }
		}
	    return(cellCountMatrix);
	  }
	
	public void addCelltype(String celltype) throws IllegalAccessException{
		this.celltypes.add(celltype);
	}
	
	/** 
	 * Get a list of all the celltypes given as input
	 */
	public List<String> getCelltypes() throws IllegalAccessException{
		if(this.celltypes == null){
			throw new IllegalAccessException("celltypes not set for this model");
		}
		return(this.celltypes);
	}

	public void SetQtlName(String qtlName){
	    this.qtlName = qtlName;
	  }
	
	public String getQtlName() throws IllegalAccessException{
		if(this.qtlName == null){
			throw new IllegalAccessException("QTL name not set for this model");
		}
	    return(this.qtlName);
	  }		
	
	/** 
	 * Set No Intercept for the model. If true, the intercept will be removed
	 */
	public void setNoIntercept(Boolean noIntercept){
	    this.noIntercept = noIntercept;
	  }
	
	/** 
	 * Get the No Intercept value for the model. If true, the intercept will be removed
	 */
	public Boolean getNoIntercept() throws IllegalAccessException{
		if(this.noIntercept == null){
			throw new IllegalAccessException("noIntercept not set for this model");
		}
	    return(this.noIntercept);
	  }
}
