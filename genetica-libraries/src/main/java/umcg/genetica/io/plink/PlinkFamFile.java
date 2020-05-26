package umcg.genetica.io.plink;

import umcg.genetica.enums.DiseaseStatus;
import umcg.genetica.enums.Gender;
import umcg.genetica.individuals.Family;
import umcg.genetica.individuals.Individual;
import umcg.genetica.io.text.TextFile;
import umcg.genetica.text.Strings;

import java.io.IOException;
import java.util.ArrayList;
import java.util.HashMap;

/**
 * Created by hwestra on 6/28/16.
 */
public class PlinkFamFile {

	ArrayList<Family> families;
	ArrayList<Individual> samples;

	public PlinkFamFile(String file) throws IOException {

		System.out.println("Parsing FAM path: " + file);


		HashMap<String, Individual> strToInd = new HashMap<String, Individual>();
		HashMap<String, Family> strToFam = new HashMap<String, Family>();
		samples = new ArrayList<>();
		families = new ArrayList<>();

		TextFile tf = new TextFile(file, TextFile.R);
		String[] elems = tf.readLineElems(Strings.whitespace);
		while (elems != null) {
			if (!elems[0].startsWith("#") && elems.length > 5) {
				String family = elems[0];
				String sample = elems[1];
				Family fam = strToFam.get(family);
				if (fam == null) {
					fam = new Family(family);
					families.add(fam);
					strToFam.put(family, fam);
				}

				Individual ind = strToInd.get(sample);
				if (ind == null) {
					ind = new Individual(sample);
					samples.add(ind);
					strToInd.put(sample, ind);
				}
				ind.setFamily(fam);
			}
			elems = tf.readLineElems(Strings.whitespace);
		}
		tf.close();

		tf.open();

		int cases = 0;
		int controls = 0;
		int unknowndisease = 0;
		int male = 0;
		int female = 0;
		int unknowngender = 0;
		int malecases = 0;
		int femalecases = 0;
		int malecontrols = 0;
		int femalecontrols = 0;

		elems = tf.readLineElems(Strings.whitespace);
		while (elems != null) {
			if (!elems[0].startsWith("#") && elems.length > 5) {
				String family = elems[0];
				String sample = elems[1];
				String fathid = elems[2];
				String mothid = elems[3];
				Gender gender = Gender.parseStatus(elems[4]);
				DiseaseStatus diseas = DiseaseStatus.parseStatus(elems[5]);

				if (gender.equals(Gender.MALE)) {
					male++;
					if (diseas.equals(DiseaseStatus.CASE)) {
						malecases++;
					} else if (diseas.equals(DiseaseStatus.CONTROL)) {
						malecontrols++;
					}
				} else if (gender.equals(Gender.FEMALE)) {
					female++;
					if (diseas.equals(DiseaseStatus.CASE)) {
						femalecases++;
					} else if (diseas.equals(DiseaseStatus.CONTROL)) {
						femalecontrols++;
					}
				} else {
					unknowngender++;
				}

				if (diseas.equals(DiseaseStatus.CASE)) {
					cases++;
				} else if (diseas.equals(DiseaseStatus.CONTROL)) {
					controls++;
				} else {
					unknowndisease++;
				}

				Individual father = null;
				Individual mother = null;
				if (!mothid.equals("0")) {
					Individual ind = strToInd.get(mothid);
					mother = ind;
				}
				if (!fathid.equals("0")) {
					Individual ind = strToInd.get(fathid);
					father = ind;
				}

				Family fam = strToFam.get(family);
				if (fam == null) {
					fam = new Family(family);
					families.add(fam);
					strToFam.put(family, fam);
				}

				Individual ind = strToInd.get(sample);
				if (ind == null) {
					ind = new Individual(sample);
					samples.add(ind);
					strToInd.put(sample, ind);
				}

				ind.setGender(gender);
				ind.setDiseaseStatus(diseas);
				ind.setFamily(fam);
				fam.addIndividual(ind);

				if (father != null) {
					ind.setFather(father);
					father.addChild(ind);
				}
				if (mother != null) {
					ind.setMother(mother);
					mother.addChild(ind);
				}
			}
			elems = tf.readLineElems(Strings.whitespace);
		}
		tf.close();

		System.out.println("Loaded: " + samples.size() + " samples in " + families.size() + " families ");
		System.out.println(cases + "\tcases and\t" + controls + "\tcontrols,\t" + unknowndisease + " unknown");
		System.out.println(male + "\tmales and\t" + female + "\tfemales,\t" + unknowngender + " unknown");
		System.out.println(malecases + "\tmale cases and\t" + malecontrols + "\tmale controls");
		System.out.println(femalecases + "\tfemale cases and\t" + femalecontrols + "\tfemale controls");
	}

	public ArrayList<Family> getFamilies() {
		return families;
	}

	public ArrayList<Individual> getSamples() {
		return samples;
	}

	public SampleAnnotation getSampleAnnotation() {
		SampleAnnotation s = new SampleAnnotation();
		s.setIndividuals(samples);
		s.setFamilies(families);
		return s;
	}
}
