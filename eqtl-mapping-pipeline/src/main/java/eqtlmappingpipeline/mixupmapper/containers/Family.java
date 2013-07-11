/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */

package eqtlmappingpipeline.mixupmapper.containers;

import java.util.HashMap;
import java.util.Vector;

/**
 *
 * @author harmjan
 */
public class Family {
    public String id;
    public Vector<Individual> individuals = new Vector<Individual>();
    public HashMap<String, Individual> sampleToIndividual = new HashMap<String, Individual>();
    public HashMap<String, Trio> sampleToTrio = new HashMap<String, Trio>();

    public Vector<Trio> trios = new Vector<Trio>();

    void makeTrios() {
	for(int i=0; i<individuals.size(); i++){
	    Individual ind = individuals.get(i);
	    if(ind.fatherId.equals("0") && ind.motherId.equals("0")){
		// this is a parent

	    } else {
		// this is a child
		Trio a = new Trio();


		a.child = ind;
		a.parent1 = sampleToIndividual.get(ind.fatherId);
		a.parent2 = sampleToIndividual.get(ind.motherId);
		System.out.println(ind.sampleName+" - fa "+ind.fatherId +" - mo "+ind.motherId);
		trios.add(a);
		sampleToTrio.put(ind.sampleName, a);
		sampleToTrio.put(a.parent1.sampleName, a);
		sampleToTrio.put(a.parent2.sampleName, a);
	    }
	}
    }

}
