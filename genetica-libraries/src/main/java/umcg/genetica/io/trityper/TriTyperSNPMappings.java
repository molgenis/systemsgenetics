/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package umcg.genetica.io.trityper;

import java.io.IOException;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.HashSet;
import java.util.Set;
import umcg.genetica.io.text.TextFile;
import umcg.genetica.io.trityper.util.ChrAnnotation;

/**
 *
 * @author harmjan
 */
public class TriTyperSNPMappings {

    private String[] snps;
    private Byte[] chr;
    private Integer[] chrPos;
    private TextFile snpFile;
    private HashMap<String, Integer> snpIndex = new HashMap<String, Integer>();

    public TriTyperSNPMappings(String loc) throws IOException {
        System.out.println("SNP annotation loading: "+loc);
        snpFile = new TextFile(loc, TextFile.R);
        parse(null);
        snpFile.close();
    }

    public TriTyperSNPMappings(String loc, String querySNPFile) throws IOException {
        System.out.println("SNP annotation loading: "+loc);
        snpFile = new TextFile(loc, TextFile.R);
        HashSet<String> querySNPs = null;
        TextFile tf = new TextFile(querySNPFile, TextFile.R);
        parse(tf.readAsSet(0, TextFile.tab));
        tf.close();
        snpFile.close();
    }

    public TriTyperSNPMappings(String loc, Set<String> querySNPs) throws IOException {
        System.out.println("SNP annotation loading: "+loc);                
        snpFile = new TextFile(loc, TextFile.R);
        parse(querySNPs);
        snpFile.close();
    }

    private void parse(Set<String> query) throws IOException {
        String[] elems = snpFile.readLineElems(TextFile.tab);
        int snpCtr = 0;
        ArrayList<String> tmpSNPs = new ArrayList<String>();
        ArrayList<Byte> tmpChr = new ArrayList<Byte>();
        ArrayList<Integer> tmpChrPos = new ArrayList<Integer>();

        while (elems != null) {
            String snp = elems[2];
            if (query == null || query.contains(snp)) {

                tmpSNPs.add(snp);
                tmpChr.add(ChrAnnotation.parseChr(elems[0]));
                tmpChrPos.add(Integer.parseInt(elems[1]));




                snpIndex.put(snp, snpCtr);
                snpCtr++;
            }
            elems = snpFile.readLineElems(TextFile.tab);
        }

        snps = tmpSNPs.toArray(new String[0]);
        chr = tmpChr.toArray(new Byte[0]);
        chrPos = tmpChrPos.toArray(new Integer[0]);
    }

    public Byte getChr(String snp) {
        Integer id = snpIndex.get(snp);
        if (id != null) {
            return chr[id];
        } else {
            return null;
        }
    }

    public Integer getChrPos(String snp) {
        Integer id = snpIndex.get(snp);
        if (id != null) {
            return chrPos[id];
        } else {
            return null;
        }
    }
    
    public Byte getChr(Integer id) {
        
        if (id != null) {
            return chr[id];
        } else {
            return null;
        }
    }

    public Integer getChrPos(Integer id) {
        if (id != null) {
            return chrPos[id];
        } else {
            return null;
        }
    }
    
    public String getSNP(Integer id) {
        if (id != null) {
            return snps[id];
        } else {
            return null;
        }
    }

    public Integer getSNPId(String snp) {
        Integer id = snpIndex.get(snp);
        return id;
    }

    public String[] getSNPs() {
        return snps;
    }

    public int getNumSNPs() {
        return snps.length;
    }
}
