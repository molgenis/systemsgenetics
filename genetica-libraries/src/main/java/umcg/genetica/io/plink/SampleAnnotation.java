package umcg.genetica.io.plink;

import cern.colt.matrix.tdouble.DoubleMatrix2D;
import umcg.genetica.enums.DiseaseStatus;
import umcg.genetica.enums.Gender;
import umcg.genetica.individuals.Family;
import umcg.genetica.individuals.Individual;

import java.util.ArrayList;
import java.util.HashMap;
import java.util.HashSet;

/**
 * Created by hwestra on 9/7/16.
 */
public class SampleAnnotation {

	private ArrayList<String> covariateNames;
	private Individual[] individuals;
	private ArrayList<Family> families;
	private String[] individualNames;
	private Gender[] individualGenders;
	private DiseaseStatus[][] sampleDiseaseStatus; // format: [samples][diseases]
	private DoubleMatrix2D covariates;

	public DoubleMatrix2D getCovariates() {
		return covariates;
	}

	public void setCovariates(DoubleMatrix2D covariates) {
		this.covariates = covariates;
	}

	public DiseaseStatus[][] getSampleDiseaseStatus() {
		if (sampleDiseaseStatus == null) {
			int nrDiseases = individuals[0].getDiseaseStatuses().length;
			sampleDiseaseStatus = new DiseaseStatus[individuals.length][nrDiseases];
			for (int i = 0; i < individuals.length; i++) {
				for (int j = 0; j < nrDiseases; j++) {
					sampleDiseaseStatus[i][j] = individuals[i].getDiseaseStatus(j);
				}
			}
		}
		return sampleDiseaseStatus;
	}

//	public void setSampleDiseaseStatus(DiseaseStatus[][] sampleDiseaseStatus) {
//		this.sampleDiseaseStatus = sampleDiseaseStatus;
//	}

//	public void setSampleDiseaseStatus(DiseaseStatus[] sampleDiseaseStatus) {
//		this.sampleDiseaseStatus = new DiseaseStatus[sampleDiseaseStatus.length][1];
//		for (int q = 0; q < sampleDiseaseStatus.length; q++) {
//			this.sampleDiseaseStatus[q][0] = sampleDiseaseStatus[q];
//		}
//	}

	public Gender[] getIndividualGender() {
		if (individualGenders == null) {
			individualGenders = new Gender[individuals.length];
			for (int i = 0; i < individuals.length; i++) {
				individualGenders[i] = individuals[i].getGender();
			}
		}
		return individualGenders;
	}

	public String[] getSampleName() {
		if (individualNames == null) {
			individualNames = new String[individuals.length];
			for (int i = 0; i < individuals.length; i++) {
				individualNames[i] = individuals[i].getName();
			}
		}
		return individualNames;
	}

	public void reorder(ArrayList<String> order, boolean pruneMissing) {

		System.out.println("Reordering sample annotation..");
		System.out.println(individuals.length + " individuals available.");
		// clear precalculated arrays
		individualNames = null;
		individualGenders = null;
		sampleDiseaseStatus = null;

		HashSet<String> set1 = new HashSet<String>();
		set1.addAll(order);

		System.out.println(set1.size() + " individuals requested");
		HashSet<String> intersect = new HashSet<>();
		for (int i = 0; i < individuals.length; i++) {
			if (set1.contains(individuals[i].getName())) {
				intersect.add(individuals[i].getName());
			}
		}
		System.out.println(intersect.size() + " overlapping individuals.");

		HashMap<String, Integer> orderMap = new HashMap<String, Integer>();
		int ctr = 0;
		for (int i = 0; i < order.size(); i++) {
			if (pruneMissing) {
				// don't add sample to order if not present in data
				boolean inData = intersect.contains(order.get(i));
				if (inData) {
					orderMap.put(order.get(i), ctr);
					ctr++;
				}
			} else {
				// add sample to order always
				orderMap.put(order.get(i), ctr);
				ctr++;
			}
		}

		// reorder samples
		Individual[] tmpInds = new Individual[orderMap.size()];
		for (int i = 0; i < individuals.length; i++) {
			Integer id = orderMap.get(individuals[i].getName());
			if (id != null) {
				tmpInds[id] = individuals[i];
			}
		}

		individuals = tmpInds;
		System.out.println("Final size of sample annotation: " + individuals.length);
		// reorder and prune covariates

//
//		private DiseaseStatus[][] sampleDiseaseStatus;
//		private ArrayList<String> covariateNames;
//		private ArrayList<Individual> individuals;
//		private ArrayList<Family> families;


	}

	public void setCovariateNames(ArrayList<String> covariateNames) {
		this.covariateNames = covariateNames;
	}

	public void setIndividuals(ArrayList<Individual> inds) {
		this.individuals = inds.toArray(new Individual[0]);
	}

	public void setFamilies(ArrayList<Family> families) {
		this.families = families;
	}
}
