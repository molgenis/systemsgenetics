/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package eqtlmappingpipeline.ase.pileup;

/**
 *
 * @author Patrick Deelen
 */
public class PileupParseException extends Exception {

	/**
	 * Creates a new instance of
	 * <code>PileupParseException</code> without detail message.
	 */
	public PileupParseException() {
	}

	/**
	 * Constructs an instance of
	 * <code>PileupParseException</code> with the specified detail message.
	 *
	 * @param msg the detail message.
	 */
	public PileupParseException(String msg) {
		super(msg);
	}
}
