package org.molgenis.genotype.variant;

import org.molgenis.genotype.util.ChromosomeComparator;
import org.molgenis.genotype.util.GeneticVariantTreeSet;

abstract public class AbstractGeneticVariant implements GeneticVariant
{

	private static final ChromosomeComparator chrComparator = new ChromosomeComparator();

	@Override
	public int compareTo(GeneticVariant other)
	{
		if (other == null)
		{
			return 0;
		}

		if (this == other)
		{
			return 0;
		}

		if (this.equals(other))
		{
			return 0;
		}

		if (!this.getSequenceName().equals(other.getSequenceName()))
		{
			return chrComparator.compare(this.getSequenceName(), other.getSequenceName());
		}
		else
		{
			// same sequence
			if (this.getStartPos() != other.getStartPos())
			{
				return this.getStartPos() - other.getStartPos();
			}
			else
			{

				// same sequence and same start

				if (GeneticVariantTreeSet.getDummyGeneticVariantClass().isInstance(this))
				{

					if (!GeneticVariantTreeSet.getDummyGeneticVariantClass().isInstance(other))
					{
						return 1;
					}
					else
					{
						return 0;
					}

				}

				if (GeneticVariantTreeSet.getDummyGeneticVariantClass().isInstance(other))
				{
					return -1;

				}

				// System.out.println(this.getClass().toString() +
				// this.getSequenceName() + ":" + this.getStartPos() + " "
				// + this.getVariantAlleles() + " vs " +
				// other.getClass().toString() + other.getSequenceName()
				// + ":" + other.getStartPos() + " " +
				// other.getVariantAlleles());

				if (!this.getVariantAlleles().equals(other.getVariantAlleles()))
				{
					// Alleles are different

					return this.getVariantAlleles().compareTo(other.getVariantAlleles());
				}
				else
				{

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
	public int hashCode()
	{
		final int prime = 31;
		int result = 0;
		result = prime * result + ((getVariantAlleles() == null) ? 0 : getVariantAlleles().hashCode());
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
	public boolean equals(Object obj)
	{
		if (this == obj) return true;
		if (!(obj instanceof GeneticVariant)) return false;
		GeneticVariant other = (GeneticVariant) obj;
		if (getSequenceName() == null)
		{
			if (other.getSequenceName() != null) return false;
		}
		else if (!getSequenceName().equals(other.getSequenceName())) return false;
		if (getStartPos() != other.getStartPos()) return false;

		// If we get here pos and sequence are identical

		if (getVariantAlleles() == null)
		{
			if (other.getVariantAlleles() != null) return false;
		}
		else if (!getVariantAlleles().equals(other.getVariantAlleles())) return false;
		if (getSampleVariantsProvider() == null)
		{
			if (other.getSampleVariantsProvider() != null) return false;
		}
		else if (!getSampleVariantsProvider().equals(other.getSampleVariantsProvider())) return false;
		return true;
	}

	@Override
	public boolean isMapped()
	{
		return !(this.getSequenceName().equals("0") && this.getStartPos() == 0);
	}
}
