/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package nl.umcg.deelenp.genotypeharmonizer;

/**
 *
 * @author Patrick Deelen
 */
public class GenotypeAlignmentException extends Exception {

	/**
	 * Creates a new instance of
	 * <code>GenotypeAlignmentException</code> without detail message.
	 */
	public GenotypeAlignmentException() {
	}

	/**
	 * Constructs an instance of
	 * <code>GenotypeAlignmentException</code> with the specified detail
	 * message.
	 *
	 * @param msg the detail message.
	 */
	public GenotypeAlignmentException(String msg) {
		super(msg);
	}
}
