package deconvolution;

import java.util.ArrayList;
import java.util.List;

public class LeastSquareModel {
	/* independentVariables = the names of the independent variables, e.g. neut%, mono%, neut%:GT */
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

	public void AddExpressionValuesGt0(double expressionValue){
		/* 
		 * Set the expression values (y) of the model. 
		 */
	    this.expressionValuesGt0.add(expressionValue);
	  }
	
	public void AddExpressionValuesGt1(double expressionValue){
		/* 
		 * Set the expression values (y) of the model. 
		 */
		this.expressionValuesGt1.add(expressionValue);
	  }
	
	public void AddExpressionValuesGt2(double expressionValue){
		/* 
		 * Set the expression values (y) of the model. 
		 */
		this.expressionValuesGt2.add(expressionValue);
	  }
	
	public double[] GetExpessionValuesGt0() throws IllegalAccessException{
		/* 
		 * Get the expression values (y) of the model. 
		 */
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
	
	public double[] GetExpessionValuesGt1() throws IllegalAccessException{
		/* 
		 * Get the expression values (y) of the model. 
		 */
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
	
	public double[] GetExpessionValuesGt2() throws IllegalAccessException{
		/* 
		 * Get the expression values (y) of the model. 
		 */
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
	
	public void AddCellCountsGt0(int celltypeIndex, double cellcount){
		/* 
		 * Set the expression values (y) of the model. 
		 */
	    this.cellCountsGt0.get(celltypeIndex).add(cellcount);
	  }
	
	public void AddCellCountsGt1(int celltypeIndex, double cellcount){
		/* 
		 * Set the expression values (y) of the model. 
		 */
		this.cellCountsGt1.get(celltypeIndex).add(cellcount);
	  }
	
	public void AddCellCountsGt2(int celltypeIndex, double cellcount){
		/* 
		 * Set the expression values (y) of the model. 
		 */
		this.cellCountsGt2.get(celltypeIndex).add(cellcount);
	  }
	
	public double[][] GetCellCountsGt0() throws IllegalAccessException{
		/* 
		 * Get the expression values (y) of the model. 
		 */
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
	
	public double[][] GetCellCountsGt1() throws IllegalAccessException{
		/* 
		 * Get the expression values (y) of the model. 
		 */
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
	
	public double[][] GetCellCountsGt2() throws IllegalAccessException{
		/* 
		 * Get the expression values (y) of the model. 
		 */
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
	
		public void AddCelltype(String celltype) throws IllegalAccessException{
		/* 
		 * Get a list of all the celltypes given as input
		 */
		this.celltypes.add(celltype);
	}
	
	public List<String> GetCelltypes() throws IllegalAccessException{
		/* 
		 * Get a list of all the celltypes given as input
		 */
		if(this.celltypes == null){
			throw new IllegalAccessException("celltypes not set for this model");
		}
		return(this.celltypes);
	}

	public void SetQtlName(String qtlName){
	    this.qtlName = qtlName;
	  }
	
	public String GetQtlName() throws IllegalAccessException{
		if(this.qtlName == null){
			throw new IllegalAccessException("QTL name not set for this model");
		}
	    return(this.qtlName);
	  }		
	
	public void SetNoIntercept(Boolean noIntercept){
		/* 
		 * Set No Intercept for the model. If true, the intercept will be removed
		 */
	    this.noIntercept = noIntercept;
	  }
	
	public Boolean GetNoIntercept() throws IllegalAccessException{
		/* 
		 * Get the No Intercept value for the model. If true, the intercept will be removed
		 */
		if(this.noIntercept == null){
			throw new IllegalAccessException("noIntercept not set for this model");
		}
	    return(this.noIntercept);
	  }
}
