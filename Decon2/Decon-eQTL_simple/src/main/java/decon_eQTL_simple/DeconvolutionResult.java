package main.java.decon_eQTL_simple;

import java.util.ArrayList;
import java.util.List;

public class DeconvolutionResult {
	private String qtlName;
	private List<Double> pvalues = new ArrayList<Double>();
	private InteractionModelCollection interactionModelCollection;
	
	public DeconvolutionResult(){};
	
	/**
	 * Set the deconvolutionResult with the InteractionModels.
	 * 
	 * @param interactionModelCollection Collection of the interaction models used to get the deconvolution result
	 * 
	 * @throws IllegalAccessException	QTL name or p-value can not be retrieved from interactionModelCollection
	 */
	public DeconvolutionResult( InteractionModelCollection interactionModelCollection) throws IllegalAccessException{
		this.qtlName = interactionModelCollection.getQtlName();
		for (String modelName : interactionModelCollection.getModelNames()){
			Double pvalue = interactionModelCollection.getPvalue(modelName);
			this.pvalues.add(pvalue);
		}
		this.interactionModelCollection = interactionModelCollection;
	}

	public String getQtlName() throws IllegalAccessException{
		if(this.qtlName == null){
			throw new IllegalAccessException("QTL name not set for this model");
		}
		return this.qtlName;
	}
	
	public List<Double> getPvalues() throws IllegalAccessException{
		if(this.pvalues == null){
			throw new IllegalAccessException("pvalues not set for this model");
		}
		return this.pvalues;
	}
	
	public InteractionModelCollection getInteractionModelCollection() throws IllegalAccessException{
		if(this.interactionModelCollection == null){
			throw new IllegalAccessException("interactionModelCollection not set");
		}
		return this.interactionModelCollection;
	}
}
