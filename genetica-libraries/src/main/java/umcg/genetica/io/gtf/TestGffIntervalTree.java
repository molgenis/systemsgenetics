package umcg.genetica.io.gtf;

import java.io.File;
import java.util.List;
import org.apache.log4j.BasicConfigurator;
import umcg.genetica.collections.intervaltree.PerChrIntervalTree;

/**
 *
 * @author Patrick Deelen
 */
public class TestGffIntervalTree {

    /**
     * @param args the command line arguments
     */
    public static void main(String[] args) throws Exception {
        
		BasicConfigurator.configure();
		
		GtfReader reader = new GtfReader(new File("D:\\UMCG\\Genetica\\Projects\\RnaSeqEqtl\\Homo_sapiens.GRCh37.71.cut.sorted.gtf"));
		
		PerChrIntervalTree<GffElement> perChrIntervalTree = reader.createIntervalTree();
		List<GffElement> elements = perChrIntervalTree.searchPosition("1", 861393);
		
		for(GffElement element : elements){
			System.out.println(element.toString());
		}
		
    }

}
