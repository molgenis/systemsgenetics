package deconvolution;

import java.util.ArrayList;
import java.util.List;

public class InteractionModel {
	private List<String> celltypes = new ArrayList<String>();
	private double[][] observedValues;
	private double[] expressionValues;
	private Boolean noIntercept;
	private String qtlName;
	private String modelName;
	InteractionModel(){};
	public InteractionModel( List<String> celltypes, double[][] observedValues){
		    this.celltypes = celltypes;
		    this.observedValues = observedValues;
		  }
	
	public void SetObservedValues( double[][] observedValues){
	    this.observedValues = observedValues;
	  }
	
	public double[][] GetObservedValues(){
	    return(this.observedValues);
	  }
	
	public void AddCelltype(String celltype){
	    this.celltypes.add(celltype);
	  }
	
	public void SetCelltypes(List<String> celltypes){
	    this.celltypes = celltypes;
	  }
	
	public List<String> GetCelltypes(){
	    return(this.celltypes);
	  }
	
	public void SetExpressionValues(double[] expressionValues){
	    this.expressionValues = expressionValues;
	  }
	
	public double[] GetExpessionValues(){
	    return(this.expressionValues);
	  }
	
	public void SetNoIntercept(Boolean noIntercept){
	    this.noIntercept = noIntercept;
	  }
	
	public Boolean GetNoIntercept(){
	    return(this.noIntercept);
	  }

	public void SetQtlName(String qtlName){
	    this.qtlName = qtlName;
	  }
	
	public String GetQtlName(){
	    return(this.qtlName);
	  }
	
	public void SetModelName(String modelName){
	    this.modelName = modelName;
	  }
	
	public String GetModelName(){
	    return(this.modelName);
	  }
		
}
