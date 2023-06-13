/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package nl.systemsgenetics.downstreamer.pathway;

import java.util.ArrayList;
import java.util.List;

/**
 *
 * @author patri
 */
public interface PathwayAnnotations {

	List<String> getAnnotationHeaders();

	List<String> getAnnotationsForPathway(String pathway);

	int getMaxNumberOfAnnotations();

	String getSetName();
	
}
