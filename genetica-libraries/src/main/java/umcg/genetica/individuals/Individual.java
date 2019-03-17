package umcg.genetica.individuals;


import umcg.genetica.enums.DiseaseStatus;
import umcg.genetica.enums.Gender;

import java.util.ArrayList;

/**
 * Created by hwestra on 6/28/16.
 */
public class Individual {

	String name;
	Individual father;
	Individual mother;
	ArrayList<Individual> children = new ArrayList<>();
	Family family;
	Gender gender;
	DiseaseStatus[] diseaseStatus;

	public Individual(String name, Individual father, Individual mother, Family family, Gender gender, DiseaseStatus diseaseStatus) {
		this.name = name;
		this.father = father;
		this.mother = mother;
		this.family = family;
		this.gender = gender;
		this.diseaseStatus = new DiseaseStatus[]{diseaseStatus};
	}

	public Individual(String name, Gender gender, DiseaseStatus diseaseStatus) {
		this.name = name;
		this.gender = gender;
		this.diseaseStatus = new DiseaseStatus[]{diseaseStatus};
	}

	public Individual(String name, Gender gender, DiseaseStatus[] diseaseStatus) {
		this.name = name;
		this.gender = gender;
		this.diseaseStatus = diseaseStatus;
	}

	public Individual(String name) {
		this.name = name;
	}

	public String getName() {
		return name;
	}

	public void setName(String name) {
		this.name = name;
	}

	public Individual getFather() {
		return father;
	}

	public void setFather(Individual father) {
		this.father = father;
	}

	public Individual getMother() {
		return mother;
	}

	public void setMother(Individual mother) {
		this.mother = mother;
	}

	public ArrayList<Individual> getChildren() {
		return children;
	}

	public void setChildren(ArrayList<Individual> children) {
		this.children = children;
	}

	public Family getFamily() {
		return family;
	}

	public void setFamily(Family family) {
		this.family = family;
	}

	public Gender getGender() {
		return gender;
	}

	public void setGender(Gender gender) {
		this.gender = gender;
	}

	public DiseaseStatus[] getDiseaseStatuses() {
		return diseaseStatus;
	}

	public DiseaseStatus getDiseaseStatus() {
		return diseaseStatus[0];
	}

	public DiseaseStatus getDiseaseStatus(int disease) {
		return diseaseStatus[disease];
	}

	public void setDiseaseStatus(DiseaseStatus diseaseStatus) {
		this.diseaseStatus = new DiseaseStatus[]{diseaseStatus};
	}

	public void setDiseaseStatus(DiseaseStatus[] diseaseStatus) {
		this.diseaseStatus = diseaseStatus;
	}

	public void addChild(Individual i) {
		children.add(i);

	}


}
