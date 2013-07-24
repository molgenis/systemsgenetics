package umcg.genetica.io.gmt;

import java.io.IOException;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.HashSet;
import java.util.List;
import java.util.Map;
import java.util.Set;
import umcg.genetica.io.text.TextFile;

/**
 *
 * @author harmjan
 *
 * $LastChangedDate: 2012-08-24 11:50:07 +0200 (Fri, 24 Aug 2012) $ 
 * $LastChangedRevision: 43 $
 */
public class GMTFile {

    List<String> pathways = new ArrayList<String>();
    Map<String, Set<String>> genesInPathways = new HashMap<String, Set<String>>();

    public GMTFile(String file) throws IOException {
        if (file == null) {
            throw new IllegalArgumentException("Error loading GMT file: no file specified.");
        } else {
            read(file);
        }
    }

    private void read(String gmtfile) throws IOException {
        System.out.println("Reading GMT file: " + gmtfile);
        TextFile tfpathway = new TextFile(gmtfile, TextFile.R);
        String[] pwelems;
        while ((pwelems = tfpathway.readLineElems(TextFile.tab)) != null) {
            if (pwelems.length > 2) {
                pathways.add(pwelems[0]);
//            System.out.println("Loading pathway: " + pwelems[0]);
                HashSet<String> genes = new HashSet<String>();
                for (int i = 2; i < pwelems.length; i++) {
                    genes.add(pwelems[i]);
                }
                genesInPathways.put(pwelems[0], genes);
            }

        }
        tfpathway.close();
        System.out.println(pathways.size() + " pathways loaded.");
    }

    public List<String> getPathways() {
        return pathways;
    }

    public Set<String> getGenesForPathway(String pathway) {
        return genesInPathways.get(pathway);
    }
    
    public Map<String, Set<String>> getPathwayToGenes() {
        return genesInPathways;
    }
}
