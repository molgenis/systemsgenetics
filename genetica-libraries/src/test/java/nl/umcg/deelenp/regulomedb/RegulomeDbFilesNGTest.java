/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package nl.umcg.deelenp.regulomedb;

import umcg.genetica.io.regulomedb.RegulomeDbEntry;
import umcg.genetica.io.regulomedb.RegulomeDbFile;
import umcg.genetica.io.regulomedb.RegulomeDbFiles;
import java.io.File;
import java.util.ArrayList;
import java.util.Iterator;
import static org.testng.Assert.*;
import org.testng.annotations.BeforeMethod;
import org.testng.annotations.Test;

/**
 *
 * @author Patrick Deelen
 */
public class RegulomeDbFilesNGTest {
	
	RegulomeDbFiles instance;
	
	public RegulomeDbFilesNGTest() {
	}

	@BeforeMethod
	public void setUpMethod() throws Exception {
		ArrayList<File> files = new ArrayList<File>();
		files.add(new File(this.getClass().getResource("/regulomeDbTestFile.txt").toURI()));
		files.add(new File(this.getClass().getResource("/regulomeDbTestFileEmpty.txt").toURI()));
		instance = new RegulomeDbFiles(files);
		instance.addRegulomeDbFile(new RegulomeDbFile(new File(this.getClass().getResource("/regulomeDbTestFile2.txt").toURI())));
	}

	/**
	 * Test of iterator method, of class RegulomeDbFiles.
	 */
	@Test
	public void testIterator() throws Exception {
		System.out.println("iterator");
		Iterator<RegulomeDbEntry> regulomeDbFileIterator = instance.iterator();
		
		RegulomeDbEntry currentEntry;
		RegulomeDbEntry expectedEntry;
		
		currentEntry = regulomeDbFileIterator.next();
		expectedEntry = new RegulomeDbEntry("chr1	20207	rs4030268	Chromatin_Structure|DNase-seq|Sknmc	5");
		assertEquals(currentEntry, expectedEntry);
		
		currentEntry = regulomeDbFileIterator.next();
		expectedEntry = new RegulomeDbEntry("chr1	29499	rs3971598	Protein_Binding|ChIP-seq|TAF1, Protein_Binding|ChIP-seq|TAF1, Protein_Binding|ChIP-seq|HEY1	5");
		assertEquals(currentEntry, expectedEntry);
		
		currentEntry = regulomeDbFileIterator.next();
		expectedEntry = new RegulomeDbEntry("chr1	11002	rs79537094	Motifs|PWM|Tcfap2c, Motifs|Footprinting|Sp3, Motifs|PWM|Sp3	6");
		assertEquals(currentEntry, expectedEntry);
		
		currentEntry = regulomeDbFileIterator.next();
		expectedEntry = new RegulomeDbEntry("chr1	11457	rs80029109	.	7");
		assertEquals(currentEntry, expectedEntry);
		
		currentEntry = regulomeDbFileIterator.next();
		expectedEntry = new RegulomeDbEntry("chr1	752566	rs3094315	Single_Nucleotides|eQTL|FLJ22639, Motifs|Footprinting|NF-kappaB, Motifs|Footprinting|NFKB1, Motifs|Footprinting|, Chromatin_Structure|DNase-seq|Ag09319, Chromatin_Structure|DNase-seq|Hipe, Chromatin_Structure|DNase-seq|Hmveclly, Chromatin_Structure|DNase-seq|Hconf, Chromatin_Structure|DNase-seq|Hmvecdblneo, Chromatin_Structure|DNase-seq|Hmvecdlyad, Chromatin_Structure|DNase-seq|Hmvecdneo, Chromatin_Structure|DNase-seq|Huvec, Chromatin_Structure|DNase-seq|Nhdfneo, Chromatin_Structure|DNase-seq|Sknmc, Chromatin_Structure|DNase-seq|Nhdfad, Chromatin_Structure|DNase-seq|Hcm, Chromatin_Structure|DNase-seq|Lncap, Chromatin_Structure|DNase-seq|Hae, Chromatin_Structure|DNase-seq|Hcf, Chromatin_Structure|DNase-seq|Hmvecdlyneo, Chromatin_Structure|DNase-seq|Ag09309, Chromatin_Structure|DNase-seq|Hac, Chromatin_Structure|DNase-seq|Hrgec, Chromatin_Structure|DNase-seq|Hgf, Chromatin_Structure|DNase-seq|Hmec, Chromatin_Structure|DNase-seq|Hasp	1f");
		assertEquals(currentEntry, expectedEntry);
		
		currentEntry = regulomeDbFileIterator.next();
		expectedEntry = new RegulomeDbEntry("chr1	96560	rs3871782	Chromatin_Structure|DNase-seq|Ag09309	5");
		assertEquals(currentEntry, expectedEntry);
		
		assertEquals(regulomeDbFileIterator.hasNext(), false);
		
	}
}
