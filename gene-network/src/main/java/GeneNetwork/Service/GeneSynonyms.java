/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package GeneNetwork.Service;

import GeneNetwork.Service.exceptions.EnsemblIdNotFoundException;
import GeneNetwork.Service.exceptions.GeneNotFoundException;
import java.io.BufferedReader;
import java.io.File;
import java.io.FileNotFoundException;
import java.io.FileReader;
import java.io.IOException;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;
import java.util.logging.Level;
import java.util.logging.Logger;

/**
 *
 * @author juha
 */
public class GeneSynonyms {

    // gene name (any synonym) to a list of unique ids
    private HashMap<String, ArrayList<Integer>> d_synonymMap;
    // gene unique id to hugo id
    private HashMap<Integer, String> d_hugoMap;
    // gene unique id to ensembl id
    private HashMap<Integer, String> d_ensemblMap;

    /**
     * 
     * Constructor calls readFile().
     * 
     * @throws IOException
     * @throws FileNotFoundException 
     */
    public GeneSynonyms(String fileName) throws IOException, FileNotFoundException {
        readFile(fileName);
    }

    /**
     * 
     * Reads the given tab-delimited file and parses gene synonyms to maps from it.
     * 
     * @param fileName file to be read
     * @throws IOException
     * @throws FileNotFoundException 
     */
    public void readFile(String fileName) throws IOException, FileNotFoundException {

        Logger.getLogger(GeneNetworkService.class.getName()).log(Level.INFO, "reading gene synonyms from file {0}", fileName);
        
        d_synonymMap = new HashMap();
        d_hugoMap = new HashMap();
        d_ensemblMap = new HashMap();

        BufferedReader br = br = new BufferedReader(new FileReader(new File(fileName)));

        String line = null;
        String[] cols = null;
        String[] synonyms = null;

        while ((line = br.readLine()) != null) {
            cols = line.split("\t", -1);
            int id = Integer.parseInt(cols[1]);
            String hugoId = cols[2].trim().toUpperCase();
            String earlierId = d_hugoMap.put(id, hugoId);
            if (earlierId != null) {
                System.err.println(id + " was already assigned to id " + earlierId + ", replaced with " + cols[2]);
            }
            addSynonyms(new String[]{hugoId}, id);
            synonyms = cols[4].split("\\|", -1);
            addSynonyms(synonyms, id);
            synonyms = cols[5].split("\\|", -1);
            addSynonyms(synonyms, id);
        }

        Logger.getLogger(GeneNetworkService.class.getName()).log(Level.INFO, "{0} synonyms read to {1} unique ids", new Object[]{d_synonymMap.size(), d_hugoMap.size()});
        Logger.getLogger(GeneNetworkService.class.getName()).log(Level.INFO, "{0} ensembl ids", d_ensemblMap.size());
    }

    /**
     * 
     * Adds the given array of synonyms and given id to a synonym-to-id map. If an Ensembl synonym exists in the array, add it to an id-to-Ensemble-id map.
     * 
     * @param synonyms synonyms to be added
     * @param id gene id
     */
    private void addSynonyms(String[] synonyms, int id) {
        for (String synonym : synonyms) {
            String trimmed = synonym.trim().toUpperCase();
            if (!"-".equals(trimmed)) {
                if (trimmed.startsWith("ENSEMBL:")) {
                    trimmed = trimmed.replace("ENSEMBL:", "");
                    d_ensemblMap.put(id, trimmed);
                }
                if (d_synonymMap.containsKey(trimmed)) {
                    d_synonymMap.get(trimmed).add(id);
                } else {
                    ArrayList<Integer> al = new ArrayList();
                    al.add(id);
                    d_synonymMap.put(trimmed, al);
                }
            }
        }
    }

    /**
     * 
     * Gets HUGO ids for a given gene synonym (there can be many).
     * 
     * @param synonym gene synonym to look for
     * @return list of HUGO ids for the synonym
     * @throws GeneNotFoundException 
     */
    public List<String> getHugoIds(String synonym) throws GeneNotFoundException {
        String trimmed = synonym.trim().toUpperCase();
        if (d_synonymMap.get(trimmed) == null) {
            throw new GeneNotFoundException("gene " + synonym + " not found");
        }
        List<Integer> ids = d_synonymMap.get(trimmed);
        List<String> synonyms = new ArrayList();
        for (int id : ids) {
            synonyms.add(d_hugoMap.get(id));
        }
        return synonyms;
    }

    /**
     * 
     * Gets Ensembl ids for a given gene synonym (there can be many).
     * 
     * @param synonym gene synonym to look for
     * @return list of Ensembl ids for the synonym
     * @throws EnsemblIdNotFoundException no Ensembl ids found for the gene
     * @throws GeneNotFoundException 
     */
    public List<String> getEnsemblIds(String synonym) throws EnsemblIdNotFoundException, GeneNotFoundException {
        String trimmed = synonym.trim().toUpperCase();
        if (d_synonymMap.get(trimmed) == null) {
            throw new GeneNotFoundException("gene " + synonym + " not found");
        }
        List<Integer> ids = d_synonymMap.get(trimmed);
        List<String> synonyms = new ArrayList();
        for (int id : ids) {
            String get = d_ensemblMap.get(id);
            if (get != null) {
                synonyms.add(d_ensemblMap.get(id));
            }
        }
        if (synonyms.isEmpty()) {
            throw new EnsemblIdNotFoundException("no ensembl ids for " + synonym + " found");
        }
        return synonyms;
    }
}
