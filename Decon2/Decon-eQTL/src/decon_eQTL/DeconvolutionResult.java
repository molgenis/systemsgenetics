package decon_eQTL;

import java.util.ArrayList;
import java.util.List;

public class DeconvolutionResult {
	private List<String> celltypes;
	private String qtlName;
	private List<Double> pvalues = new ArrayList<Double>();
	private InteractionModelCollection interactionModelCollection;
	private double wholeBloodQTL;
	private double wholeBloodQTLpvalue;
	
	public DeconvolutionResult(){};
	
	/**
	 * Set the deconvolutionResult with the InteractionModels.
	 * 
	 * @param interactionModelCollection Collection of the interaction models used to get the deconvolution result
	 * 
	 * @param wholeBloodQTL Spearman correlation of genotypes and expression levels
	 * 
	 * @param wholeBloodQTLpvalue p-value of the spearman correlation of genotypes and expression levels
	 * 
	 * @throws IllegalAccessException	QTL name or p-value can not be retrieved from interactionModelCollection
	 */
	public DeconvolutionResult( InteractionModelCollection interactionModelCollection, double wholeBloodQTL, 
								double wholeBloodQTLpvalue) throws IllegalAccessException{

		celltypes = interactionModelCollection.getAllCelltypes();
		this.qtlName = interactionModelCollection.getQtlName();
		for (int i = 0; i < celltypes.size(); i++){
			String modelName = celltypes.get(i);
			Double pvalue = interactionModelCollection.getPvalue(modelName);
			this.pvalues.add(pvalue);
		}
		this.interactionModelCollection = interactionModelCollection;
		this.wholeBloodQTL = wholeBloodQTL;
		this.wholeBloodQTLpvalue = wholeBloodQTLpvalue;
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

	public double  getWholeBloodQTL() throws IllegalAccessException{
		return this.wholeBloodQTL;
	}
	public double  getWholeBloodQTLpvalue() throws IllegalAccessException{
		return this.wholeBloodQTLpvalue;
	}
}
