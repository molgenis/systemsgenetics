package deconvolution;

import java.util.ArrayList;
import java.util.List;

public class InteractionModel {
	private List<String> independentVariables = new ArrayList<String>();
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
	
	public double[][] GetObservedValues() throws IllegalAccessException{
		if(this.observedValues == null){
			throw new IllegalAccessException("observedValues not set for this model");
		}
	    return(this.observedValues);
	  }
	
	public void AddIndependentVariable(String independentVariables){
	    this.independentVariables.add(independentVariables);
	  }
	
	public void AddIndependentVariable(int index, String independentVariables){
	    this.independentVariables.add(index, independentVariables);
	  }
	
	public void SetIndependentVariable(List<String> independentVariables){
	    this.independentVariables = independentVariables;
	  }
	
	public List<String> GetIndependentVariables() throws IllegalAccessException{
		if(this.independentVariables == null){
			throw new IllegalAccessException("celltypes not set for this model");
		}
	    return(this.independentVariables);
	  }
	
	public void SetExpressionValues(double[] expressionValues){
	    this.expressionValues = expressionValues;
	  }
	
	public double[] GetExpessionValues() throws IllegalAccessException{
		if(this.expressionValues == null){
			throw new IllegalAccessException("expressionValues not set for this model");
		}
	    return(this.expressionValues);
	  }
	
	public void SetNoIntercept(Boolean noIntercept){
	    this.noIntercept = noIntercept;
	  }
	
	public Boolean GetNoIntercept() throws IllegalAccessException{
		if(this.noIntercept == null){
			throw new IllegalAccessException("noIntercept not set for this model");
		}
	    return(this.noIntercept);
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
	
	public void SetModelName(String modelName){
	    this.modelName = modelName;
	  }
	
	public String GetModelName() throws IllegalAccessException{
		if(this.modelName == null){
			throw new IllegalAccessException("modelName name not set for this model");
		}
	    return(this.modelName);
	  }
		
}
