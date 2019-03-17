package umcg.genetica.individuals;

import java.util.ArrayList;

/**
 * Created by hwestra on 6/28/16.
 */
public class Family {

	String name;
	ArrayList<Individual> individuals = new ArrayList<>();

	public String getName() {
		return name;
	}

	public void setName(String name) {
		this.name = name;
	}

	public ArrayList<Individual> getIndividuals() {
		return individuals;
	}

	public void setIndividuals(ArrayList<Individual> individuals) {
		this.individuals = individuals;
	}

	public Family(String name) {
		this.name = name;
		individuals = new ArrayList<>();
	}

	public Family(String name, ArrayList<Individual> individuals) {
		this.name = name;
		this.individuals = individuals;
	}

	public void addIndividual(Individual i) {
		individuals.add(i);
	}

}
