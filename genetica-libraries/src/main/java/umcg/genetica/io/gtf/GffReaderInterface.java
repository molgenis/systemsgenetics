/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package umcg.genetica.io.gtf;

import java.io.IOException;

/**
 *
 * @author Patrick Deelen
 */
public interface GffReaderInterface extends Iterable<GffElement>{
	
	/**
	 * 
	 * 
	 * @return next element null when done
	 */
	public GffElement nextElement() throws IOException;
	
	
	public void close() throws IOException;
	
}
