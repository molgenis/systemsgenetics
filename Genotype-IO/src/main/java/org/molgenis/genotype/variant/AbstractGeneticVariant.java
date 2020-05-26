package org.molgenis.genotype.variant;

import org.molgenis.genotype.Allele;
import org.molgenis.genotype.Alleles;
import org.molgenis.genotype.GenotypeDataException;
import org.molgenis.genotype.util.ChromosomeComparator;
import org.molgenis.genotype.util.Ld;
import org.molgenis.genotype.util.LdCalculator;
import org.molgenis.genotype.util.LdCalculatorException;
import org.molgenis.genotype.util.MachR2Calculator;

abstract public class AbstractGeneticVariant implements GeneticVariant {

    private static final ChromosomeComparator chrComparator = new ChromosomeComparator();
	
    @Override
    public final int compareTo(GeneticVariant other) {
        if (other == null) {
            return 0;
        }

        if (this == other) {
            return 0;
        }

        if (this.equals(other)) {
            return 0;
        }

        if (!this.getSequenceName().equals(other.getSequenceName())) {
            return chrComparator.compare(this.getSequenceName(), other.getSequenceName());
        } else {
            // same sequence
            if (this.getStartPos() != other.getStartPos()) {
                return this.getStartPos() < other.getStartPos() ? -1 : (this.getStartPos() == other.getStartPos() ? 0 : 1);
            } else {

                // same sequence and same start

                if (ReadOnlyGeneticVariantTriTyper.class.isInstance(this)) {

                    if (!ReadOnlyGeneticVariantTriTyper.class.isInstance(other)) {
                        return 1;
                    } else {
                        //Both are ReadOnlyGeneticVariantTriTyper
                        return this.getSampleVariantsProvider().getSampleVariantProviderUniqueId()
                                - other.getSampleVariantsProvider().getSampleVariantProviderUniqueId();
                    }

                }

                if (ReadOnlyGeneticVariantTriTyper.class.isInstance(other)) {
                    return -1;

                }

                // System.out.println(this.getClass().toString() +
                // this.getSequenceName() + ":" + this.getStartPos() + " "
                // + this.getVariantAlleles() + " vs " +
                // other.getClass().toString() + other.getSequenceName()
                // + ":" + other.getStartPos() + " " +
                // other.getVariantAlleles());

                if (!this.getVariantAlleles().equals(other.getVariantAlleles())) {
                    // Alleles are different

                    return this.getVariantAlleles().compareTo(other.getVariantAlleles());
                } else {

                    return this.getSampleVariantsProvider().getSampleVariantProviderUniqueId()
                            - other.getSampleVariantsProvider().getSampleVariantProviderUniqueId();

                }

            }
        }
    }

    /*
     * (non-Javadoc)
     * 
     * @see java.lang.Object#hashCode()
     */
    @Override
    public int hashCode() {
        final int prime = 31;
        int result = 0;
        result = prime * result + ((getSampleVariantsProvider() == null) ? 0 : getSampleVariantsProvider().hashCode());
        result = prime * result + ((getSequenceName() == null) ? 0 : getSequenceName().hashCode());
        result = prime * result + getStartPos();
        return result;
    }

    /*
     * (non-Javadoc)
     * 
     * @see java.lang.Object#equals(java.lang.Object)
     */
    @Override
    public boolean equals(Object obj) {
        if (this == obj) {
            return true;
        }
        if (!(obj instanceof GeneticVariant)) {
            return false;
        }
        GeneticVariant other = (GeneticVariant) obj;
        if (getSequenceName() == null) {
            if (other.getSequenceName() != null) {
                return false;
            }
        } else if (!getSequenceName().equals(other.getSequenceName())) {
            return false;
        }
        if (getStartPos() != other.getStartPos()) {
            return false;
        }

        // If we get here pos and sequence are identical

        if (!other.getClass().equals(obj.getClass())) {
            return false;
        }

        if (getVariantAlleles() == null) {
            if (other.getVariantAlleles() != null) {
                return false;
            }
        } else if (!getVariantAlleles().equals(other.getVariantAlleles())) {
            return false;
        }
        if (getSampleVariantsProvider() == null) {
            if (other.getSampleVariantsProvider() != null) {
                return false;
            }
        } else if (!getSampleVariantsProvider().equals(other.getSampleVariantsProvider())) {
            return false;
        }
        return true;
    }

    @Override
    public boolean isMapped() {
        return !(this.getSequenceName().equals("0") || this.getStartPos() < 0);
    }

    @Override
    public double getCallRate() {

        int total = 0;
        int call = 0;

        for (Alleles sampleAlleles : this.getSampleVariants()) {
            ++total;

            if (!sampleAlleles.contains(Allele.ZERO)) {
                ++call;
            }

        }

        return (double) call / total;
    }

    @Override
    public double getHwePvalue() {

        if (!this.isBiallelic()) {
            return Double.NaN;
        }

        Alleles hom1 = Alleles.createAlleles(this.getVariantAlleles().get(0), this.getVariantAlleles().get(0));
        Alleles het1 = Alleles.createAlleles(this.getVariantAlleles().get(1), this.getVariantAlleles().get(0));
        Alleles het2 = Alleles.createAlleles(this.getVariantAlleles().get(0), this.getVariantAlleles().get(1));
        Alleles hom2 = Alleles.createAlleles(this.getVariantAlleles().get(1), this.getVariantAlleles().get(1));

        int obs_hom1 = 0;
        int obs_hets = 0;
        int obs_hom2 = 0;

        for (Alleles sampleAlleles : this.getSampleVariants()) {

            if (sampleAlleles == hom1) {
                ++obs_hom1;
            } else if (sampleAlleles == het1 || sampleAlleles == het2) {
                ++obs_hets;
            } else if (sampleAlleles == hom2) {
                ++obs_hom2;
            }

        }

        //Code below from genetica lib

        int obs_homc = obs_hom1 < obs_hom2 ? obs_hom2 : obs_hom1;
        int obs_homr = obs_hom1 < obs_hom2 ? obs_hom1 : obs_hom2;

        int rare_copies = 2 * obs_homr + obs_hets;
        int l_genotypes = obs_hets + obs_homc + obs_homr;

        if (l_genotypes == 0) {
            return Double.NaN;
        }

        double[] het_probs = new double[rare_copies + 1];

        int i;
        for (i = 0; i <= rare_copies; i++) {
            het_probs[i] = 0.0;
        }

        /* start at midpoint */
        int mid = rare_copies * (2 * l_genotypes - rare_copies) / (2 * l_genotypes);

        /* check to ensure that midpoint and rare alleles have same parity */
        if (mid % 2 != rare_copies % 2) {
            mid++;
        }

        int curr_hets = mid;
        int curr_homr = (rare_copies - mid) / 2;
        int curr_homc = l_genotypes - curr_hets - curr_homr;

        het_probs[mid] = 1.0;
        double sum = het_probs[mid];
        for (curr_hets = mid; curr_hets > 1; curr_hets -= 2) {
            het_probs[curr_hets - 2] = het_probs[curr_hets] * curr_hets * (curr_hets - 1.0) / (4.0 * (curr_homr + 1.0) * (curr_homc + 1.0));
            sum += het_probs[curr_hets - 2];
            /* 2 fewer heterozygotes for next iteration -> add one rare, one common homozygote */
            curr_homr++;
            curr_homc++;
        }

        curr_hets = mid;
        curr_homr = (rare_copies - mid) / 2;
        curr_homc = l_genotypes - curr_hets - curr_homr;
        for (curr_hets = mid; curr_hets <= rare_copies - 2; curr_hets += 2) {
            het_probs[curr_hets + 2] = het_probs[curr_hets] * 4.0 * curr_homr * curr_homc / ((curr_hets + 2.0) * (curr_hets + 1.0));
            sum += het_probs[curr_hets + 2];
            curr_homr--;
            curr_homc--;
        }

        for (i = 0; i <= rare_copies; i++) {
            het_probs[i] /= sum;
        }

        double p_hwe = 0.0;
        for (i = 0; i <= rare_copies; i++) {
            if (het_probs[i] <= het_probs[obs_hets]) {
                p_hwe += het_probs[i];
            }
        }

        p_hwe = p_hwe > 1.0 ? 1.0 : p_hwe;

        return p_hwe;

    }

    @Override
    public boolean isSnp() {
        return this.getVariantAlleles().isSnp();
    }

    @Override
    public boolean isAtOrGcSnp() {
        return this.getVariantAlleles().isAtOrGcSnp();
    }

    @Override
    public Ld calculateLd(GeneticVariant other) throws LdCalculatorException {
        return LdCalculator.calculateLd(this, other);
    }

    @Override
    public boolean isBiallelic() {
        return this.getVariantAlleles().getAlleleCount() == 2;
    }
	
	@Override
	public double getMachR2() {
		return MachR2Calculator.calculateMachR2(this.getSampleGenotypeProbilities());
	}
}
