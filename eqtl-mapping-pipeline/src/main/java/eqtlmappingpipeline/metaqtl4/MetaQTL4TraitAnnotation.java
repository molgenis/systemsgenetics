/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package eqtlmappingpipeline.metaqtl4;

import java.io.File;
import java.io.IOException;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.HashSet;
import java.util.Set;
import java.util.TreeSet;
import umcg.genetica.io.text.TextFile;

/**
 *
 * @author harmjan
 */
class MetaQTL4TraitAnnotation {

    private final String[] platforms;
    private final ArrayList<HashMap<String, MetaQTL4MetaTrait>> traitHashPerPlatform; // maps each probe on a specific platform to a specific meta-analysis object
    private final HashMap<String, Integer> platformToId; // maps a platform to an id
    private final MetaQTL4MetaTraitTreeSet metatraits; // lists all meta-analysis objects
    private final HashMap<String, MetaQTL4MetaTrait> metaTraitNameToObj; // maps a meta-analsysis object name to a specific meta-analysis object

    public MetaQTL4TraitAnnotation(File probeAnnotationFile, Set<String> platformsToInclude) throws IOException {
        TextFile tf = new TextFile(probeAnnotationFile, TextFile.R);
        int nrLines = tf.countLines();
        String[] header = tf.readLineElems(TextFile.tab);
        boolean[] colsToInclude = new boolean[header.length];
        int nrPlatforms = 0;
        HashSet<String> visitedPlatforms = new HashSet<String>();
        for (int i = 0; i < header.length; i++) {
            if (platformsToInclude.contains(header[i])) {
                if (!visitedPlatforms.contains(header[i])) {
                    colsToInclude[i] = true;
                    visitedPlatforms.add(header[i]);
                    nrPlatforms++;
                } else {
                    System.err.println("ERROR: your probe annotation file contains duplicate platform identifiers!");
                }

            }
        }

        metatraits = new MetaQTL4MetaTraitTreeSet();
        metaTraitNameToObj = new HashMap<String, MetaQTL4MetaTrait>();
        platformToId = new HashMap<String, Integer>();
        platforms = new String[nrPlatforms];
        traitHashPerPlatform = new ArrayList<HashMap<String, MetaQTL4MetaTrait>>();

        int platformNr = 0;
        for (int i = 0; i < header.length; i++) {
            if (colsToInclude[i]) {
                platforms[platformNr] = header[i];
                platformToId.put(header[i], platformNr);
                HashMap<String, MetaQTL4MetaTrait> probeToId = new HashMap<String, MetaQTL4MetaTrait>();
                traitHashPerPlatform.add(probeToId);
                platformNr++;
            }
        }

        for (String platform : platformsToInclude) {
            if (!visitedPlatforms.contains(platform)) {
                System.err.println("WARNING: no annotation will be loaded for platform: " + platform);
            }
        }


        int probeCounter = 0;
        for (String[] elems : tf.readLineElemsIterable(TextFile.tab)) {

            String metaTraitName = new String(elems[0].getBytes("UTF-8")).intern();
            String chr = new String(elems[2].getBytes("UTF-8")).intern();
            String chrpos = elems[3];
            String[] chrposElems = chrpos.split("-");
            int chrstartpos = -1;
            int chrendpos = -1;
            if (chrposElems.length >= 1) {
                try {
                    chrstartpos = Integer.parseInt(chrposElems[0]);
                } catch (NumberFormatException e) {
                }
                try {
                    chrendpos = Integer.parseInt(chrposElems[chrposElems.length - 1]);
                } catch (NumberFormatException e) {
                }
            }

            String hugo = new String(elems[4].getBytes("UTF-8")).intern();
            String[] platformIds = new String[nrPlatforms];
            // int metaTraitId, String metaTraitName, String chr, int chrStart, int chrEnd, String annotation, String[] platformIds
            MetaQTL4MetaTrait metaTraitObj = new MetaQTL4MetaTrait(probeCounter, metaTraitName, chr, chrstartpos, chrendpos, hugo, platformIds);

            for (int i = 5; i < elems.length; i++) {
                platformNr = 0;
                if (colsToInclude[i]) {
                    platformIds[platformNr] = new String(elems[i].getBytes("UTF-8")).intern();
                    HashMap<String, MetaQTL4MetaTrait> probeToId = traitHashPerPlatform.get(platformNr);
                    probeToId.put(elems[i], metaTraitObj);
                    platformNr++;
                }
            }
            probeCounter++;
            metatraits.add(metaTraitObj);
            metaTraitNameToObj.put(metaTraitName, metaTraitObj);
        }
        System.out.println(tf.getFileName() + " has annotation for " + visitedPlatforms.size() + " platforms and " + metatraits.size() + " traits.");
        tf.close();
    }

    public Integer getPlatformId(String platform) {
        return platformToId.get(platform);
    }

    public MetaQTL4MetaTrait getTraitForPlatformId(Integer platformId, String platformTrait) {
        return traitHashPerPlatform.get(platformId).get(platformTrait);
    }

    public String[] getPlatforms() {
        return platforms;
    }

    public ArrayList<HashMap<String, MetaQTL4MetaTrait>> getTraitHashPerPlatform() {
        return traitHashPerPlatform;
    }

    public HashMap<String, Integer> getPlatformToId() {
        return platformToId;
    }

    public TreeSet<MetaQTL4MetaTrait> getMetatraits() {
        return metatraits;
    }

    public HashMap<String, MetaQTL4MetaTrait> getMetaTraitNameToObj() {
        return metaTraitNameToObj;
    }
}
