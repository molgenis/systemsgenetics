package eqtlmappingpipeline.ase;

import cern.colt.list.tdouble.DoubleArrayList;
import cern.colt.list.tint.IntArrayList;
import java.util.ArrayList;
import org.molgenis.genotype.Allele;
import org.molgenis.genotype.variant.id.GeneticVariantId;

/**
 *
 * @author Patrick Deelen
 */
public interface AseVariant extends Comparable<AseVariant> {

	Allele getA1();

	IntArrayList getA1Counts();

	Allele getA2();

	IntArrayList getA2Counts();

	String getChr();

	double getCountPearsonR();

	GeneticVariantId getId();

	double getMetaPvalue();

	double getMetaZscore();

	AseMle getMle();

	DoubleArrayList getPValues();

	int getPos();

	int getSampleCount();

	ArrayList<String> getSampleIds();

}
