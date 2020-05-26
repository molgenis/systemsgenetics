/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package org.molgenis.genotype.variantFilter;

import java.util.HashSet;
import org.molgenis.genotype.variant.GeneticVariant;

/**
 *
 * @author Patrick Deelen
 */
public class VariantFilterSeq implements VariantFilter {

    private final HashSet<String> toIncludeSeqs;

    public VariantFilterSeq() {
        this.toIncludeSeqs = new HashSet<String>();
    }

    public VariantFilterSeq(String... seq) {
        this.toIncludeSeqs = new HashSet<String>();
        for (String s : seq) {
            toIncludeSeqs.add(s);
        }
    }

    public VariantFilterSeq(HashSet<String> toIncludeSeqs) {
        this.toIncludeSeqs = toIncludeSeqs;
    }

    public void addSeq(String seq) {
        toIncludeSeqs.add(seq);
    }

    @Override
    public boolean doesVariantPassFilter(GeneticVariant variant) {
        return toIncludeSeqs.contains(variant.getSequenceName());
    }

    @Override
    public boolean doesIdPassFilter(String id) {
        return true;
    }

}
