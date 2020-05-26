/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package eqtlmappingpipeline.ase;

/**
 *
 * @author Patrick Deelen
 */
public class AseException extends Exception {

	/**
	 * Creates a new instance of
	 * <code>AseException</code> without detail message.
	 */
	public AseException() {
	}

	/**
	 * Constructs an instance of
	 * <code>AseException</code> with the specified detail message.
	 *
	 * @param msg the detail message.
	 */
	public AseException(String msg) {
		super(msg);
	}
}
