/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package umcg.genetica.io.gwascatalog;

import umcg.genetica.io.text.TextFile;

import java.io.IOException;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.HashMap;
import java.util.HashSet;

/**
 * @author harmjan
 */
public class GWASCatalog {
	
	private HashSet<GWASLocus> loci = new HashSet<GWASLocus>();
	private HashSet<GWASSNP> snps = new HashSet<GWASSNP>();
	private HashSet<GWASTrait> traits = new HashSet<GWASTrait>();
	private GWASTrait[] traitArray = null;
	private HashMap<String, GWASPublication> publicationToObj = new HashMap<String, GWASPublication>();
	private HashMap<String, GWASSNP> snpToObj = new HashMap<String, GWASSNP>();
	private HashMap<String, GWASLocus> locusToObj = new HashMap<String, GWASLocus>();
	private HashMap<String, GWASTrait> traitToObj = new HashMap<String, GWASTrait>();
	private HashMap<String, GWASTrait> cleanTraitToObj = new HashMap<String, GWASTrait>();
	
	public GWASCatalog() {
	}
	
	public GWASCatalog(String gwasCatalogLoc) throws IOException {
		this.read(gwasCatalogLoc);
	}
	
	public GWASCatalog(String gwasCatalogLoc, double pvaluethreshold) throws IOException {
		this.read(gwasCatalogLoc);
	}
	
	public void read(String calatogloc) throws IOException {
		TextFile tf = new TextFile(calatogloc, TextFile.R);
		String[] headerelems = tf.readLineElemsReturnReference(TextFile.tab);
		
		int dateAddedCol = -1;
		int pubMedidCol = -1;
		int firstAuthorCol = -1;
		int publishDateCol = -1;
		int journalCol = -1;
		int studyCol = -1;
		
		int diseaseCol = -1;
		int samplesizeCol = -1;
		int samplesizeReplicationCol = -1;
		
		int topSNPCol = -1;
		int snpCol = -1;
		
		int pvalCol = -1;
		int chrCol = -1;
		int chrPosCol = -1;
		
		int reportedGeneCol = -1;
		int mappedGeneCol = -1;
		
		int col = 0;
		for (String e : headerelems) {
			
			e = e.toLowerCase();
			
			if (e.equals("Date Added to Catalog".toLowerCase())) {
				dateAddedCol = col;
			} else if (e.equals("PUBMEDID".toLowerCase())) {
				pubMedidCol = col;
			} else if (e.equals("First Author".toLowerCase())) {
				firstAuthorCol = col;
			} else if (e.equals("Date".toLowerCase())) {
				publishDateCol = col;
			} else if (e.equals("Journal".toLowerCase())) {
				journalCol = col;
			} else if (e.equals("Study".toLowerCase())) {
				studyCol = col;
			} else if (e.equals("Disease/Trait".toLowerCase())) {
				diseaseCol = col;
			} else if (e.equals("Initial Sample Size".toLowerCase())) {
				samplesizeCol = col;
			} else if (e.equals("Replication Sample Size".toLowerCase())) {
				samplesizeReplicationCol = col;
			} else if (e.equals("Strongest SNP-Risk Allele".toLowerCase())) {
				topSNPCol = col;
			} else if (e.equals("SNPs".toLowerCase())) {
				snpCol = col;
			} else if (e.equals("p-Value".toLowerCase())) {
				pvalCol = col;
			} else if (e.equals("Chr_id".toLowerCase())) {
				chrCol = col;
			} else if (e.equals("Chr_pos".toLowerCase())) {
				chrPosCol = col;
			} else if (e.equals("Reported Gene(s)".toLowerCase())) {
				reportedGeneCol = col;
			} else if (e.equals("Mapped_gene".toLowerCase())) {
				mappedGeneCol = col;
			} else if (e.equals("LINK")) {
				// below is new stuff that is added by EBI (might be useful later on?)
			} else if (e.equals("STUDY")) {
			} else if (e.equals("REGION")) {
			} else if (e.equals("UPSTREAM_GENE_ID")) {
			} else if (e.equals("DOWNSTREAM_GENE_ID")) {
			} else if (e.equals("SNP_GENE_IDS")) {
			} else if (e.equals("UPSTREAM_GENE_DISTANCE")) {
			} else if (e.equals("DOWNSTREAM_GENE_DISTANCE")) {
			} else if (e.equals("MERGED")) {
			} else if (e.equals("SNP_ID_CURRENT")) {
			} else if (e.equals("CONTEXT")) {
			} else if (e.equals("INTERGENIC")) {
			} else if (e.equals("RISK ALLELE FREQUENCY")) {
			} else if (e.equals("PVALUE_MLOG")) {
			} else if (e.equals("P-VALUE (TEXT)")) {
			} else if (e.equals("OR or BETA")) {
			} else if (e.equals("95% CI (TEXT)")) {
			} else if (e.equals("PLATFORM [SNPS PASSING QC]")) {
			} else if (e.equals("CNV")) {
			
			}
			
			
			// v2
            /*
            DATE ADDED TO CATALOG
PUBMEDID
FIRST AUTHOR
DATE
JOURNAL
LINK
STUDY
DISEASE/TRAIT
INITIAL SAMPLE SIZE
REPLICATION SAMPLE SIZE
REGION
CHR_ID
CHR_POS
REPORTED GENE(S)
MAPPED_GENE
UPSTREAM_GENE_ID
DOWNSTREAM_GENE_ID
SNP_GENE_IDS
UPSTREAM_GENE_DISTANCE
DOWNSTREAM_GENE_DISTANCE
STRONGEST SNP-RISK ALLELE
SNPS
MERGED
SNP_ID_CURRENT
CONTEXT
INTERGENIC
RISK ALLELE FREQUENCY
P-VALUE
PVALUE_MLOG
P-VALUE (TEXT)
OR or BETA
95% CI (TEXT)
PLATFORM [SNPS PASSING QC]
CNV

             */
			
			col++;
		}
		
		String[] elems = tf.readLineElemsReturnReference(TextFile.tab);
		int numtraits = 0;
		int numsnps = 0;
		int numpubs = 0;
		while (elems != null) {
			if (elems.length > 11) {
				String pubname = elems[pubMedidCol] + "; " + elems[firstAuthorCol] + "; " + elems[publishDateCol] + "; " + elems[journalCol] + "; " + elems[studyCol];
//                String pubname = elems[pubMedidCol];
//	    int studysize = Integer.parseInt(elems[samplesizeCol]);
//	    int studySizeReplication = Integer.parseInt(elems[samplesizeReplicationCol]);
				
				String trait = elems[diseaseCol].trim();
				String cleanedTrait = trait.replaceAll(" ", "_").replaceAll("[^a-zA-Z0-9\\-_]+", "");
				String otherSNPs = elems[snpCol].trim();
				String[] topSNPElems = elems[topSNPCol].split("-");
				String riskallele = null;
				
				String mGene = elems[mappedGeneCol];
				String rGene = elems[reportedGeneCol];
				
				HashSet<String> mappedGenes = new HashSet();
				HashSet<String> reportedGenes = new HashSet();
				if (!mGene.equals("NR") && !mGene.equals("Intergenic") && !mGene.equals(" - ")) {
					if (mGene.contains(" - ")) {
						String[] mGenes = mGene.split(" - ");
						mappedGenes.addAll(Arrays.asList(mGenes));
					} else if (mGene.contains(";")) {
						String[] mGenes = mGene.split(";");
						mappedGenes.addAll(Arrays.asList(mGenes));
					} else {
						mappedGenes.add(mGene);
					}
				}
				if (!rGene.equals("NR") && !rGene.equals("Intergenic") && !rGene.equals(" - ")) {
					if (rGene.contains(" - ")) {
						String[] rGenes = rGene.split(" - ");
						reportedGenes.addAll(Arrays.asList(rGenes));
					} else if (rGene.contains(";")) {
						String[] rGenes = rGene.split(";");
						reportedGenes.addAll(Arrays.asList(rGenes));
					} else {
						reportedGenes.add(rGene);
					}
				}
				
				byte chr = -1;
				int chrPos = -1;
				try {
					chr = Byte.parseByte(elems[chrCol]);
					chrPos = Integer.parseInt(elems[chrPosCol]);
				} catch (NumberFormatException ex) {
					//System.out.println("Chromosome and/or position unparseable for trait: " + trait + " associated with SNP " + snp + ": chr: " + elems[chrCol] + ", pos: " + elems[chrPosCol]);
				}
				if (topSNPElems.length > 1) {
					riskallele = topSNPElems[1];
					if (riskallele.equals("?")) {
						riskallele = null;
					}
				}
				
				GWASPublication pub = publicationToObj.get(pubname);
				if (pub == null) {
					pub = new GWASPublication();
					pub.id = numpubs;
					pub.name = pubname;
					publicationToObj.put(pubname, pub);
					numpubs++;
				}
				
				GWASTrait gwasTraitObj = traitToObj.get(trait);
				if (gwasTraitObj == null) {
					gwasTraitObj = new GWASTrait();
					gwasTraitObj.name = trait;
					gwasTraitObj.cleanName = cleanedTrait;
					gwasTraitObj.id = numtraits;
					gwasTraitObj.setMappedGenes(mappedGenes);
					gwasTraitObj.setReportedGenes(reportedGenes);
					traitToObj.put(trait, gwasTraitObj);
					cleanTraitToObj.put(cleanedTrait, gwasTraitObj);
					traits.add(gwasTraitObj);
					numtraits++;
				}
				
				// parse the top SNP: remove whitespace..
				String topSNP = topSNPElems[0];
				topSNP = topSNP.trim();
				while (topSNP.startsWith(" ")) {
					topSNP = topSNP.substring(1);
				}
				
				GWASSNP gwasTopSNPObj = snpToObj.get(topSNP);
				if (gwasTopSNPObj == null) {
					gwasTopSNPObj = new GWASSNP();
					gwasTopSNPObj.setName(topSNP);
					gwasTopSNPObj.setId(numsnps);
					gwasTopSNPObj.setChr(chr);
					gwasTopSNPObj.setPosition(chrPos);
					snpToObj.put(topSNP, gwasTopSNPObj);
					snps.add(gwasTopSNPObj);
					numsnps++;
				}
				
				Double topSNPAssocPVal = null;
				try {
					topSNPAssocPVal = Double.parseDouble(elems[pvalCol]);
				} catch (NumberFormatException e) {
					// Sometimes the pvalue is unreported...
					// System.out.println("P-value unparseable for trait: " + gwasTraitObj.getName() + " associated with SNP " + gwasSNPObj.getName() + ": " + elems[pvalCol]);
				}
				gwasTopSNPObj.getAssociatedTraits().add(gwasTraitObj);
				gwasTraitObj.addTopSNP(gwasTopSNPObj);
				
				if (topSNPAssocPVal != null) {
					Double previousP = gwasTopSNPObj.getPValueAssociatedWithTrait(gwasTraitObj);
					if (previousP == null || previousP > topSNPAssocPVal) {
						gwasTopSNPObj.setPValueAssociatedWithTrait(gwasTraitObj, topSNPAssocPVal);
						gwasTopSNPObj.getRiskAllele().put(gwasTraitObj, riskallele);
					}
				}
				// parse all the other reported SNPs..
				String[] otherSNPElems = otherSNPs.split(",");
				for (int s = 0; s < otherSNPElems.length; s++) {
					
					String snpname = otherSNPElems[s].trim();
					while (snpname.startsWith(" ")) {
						snpname = snpname.substring(1);
					}
					
					GWASSNP gwasSNPObj = snpToObj.get(snpname);
					if (gwasSNPObj == null) {
						gwasSNPObj = new GWASSNP();
						gwasSNPObj.setName(snpname);
						gwasSNPObj.setId(numsnps);
						gwasSNPObj.setChr(chr);
						gwasSNPObj.setPosition(chrPos);
						snpToObj.put(snpname, gwasSNPObj);
						snps.add(gwasSNPObj);
						numsnps++;
					}
					
					// The GWAS Catalog often only publishes a single p-value for a couple of SNPs.
					// We'll assume that all the reported SNPs have an LD ~ 1.0
					Double pval = null;
					try {
						pval = Double.parseDouble(elems[pvalCol]);
					} catch (NumberFormatException e) {
						//System.out.println("P-value unparseable for trait: " + gwasTraitObj.getName() + " associated with SNP " + gwasSNPObj.getName() + ": " + elems[pvalCol]);
					}
					gwasSNPObj.getAssociatedTraits().add(gwasTraitObj);
					
					if (pval != null) {
						Double previousP = gwasSNPObj.getPValueAssociatedWithTrait(gwasTraitObj);
						if (previousP == null || previousP > pval) {
							gwasSNPObj.setPValueAssociatedWithTrait(gwasTraitObj, pval);
							gwasSNPObj.getRiskAllele().put(gwasTraitObj, riskallele);
						}
					}
					
					gwasTraitObj.snps.add(gwasSNPObj);
					pub.snps.add(gwasSNPObj);
					pub.setPValueAssociatedWithTrait(gwasSNPObj, gwasTraitObj, pval);
					gwasSNPObj.getPublishedIn().add(pub);
				}
				
				gwasTraitObj.appendMappedGenes(mappedGenes);
				gwasTraitObj.appendReportedGenes(reportedGenes);
				gwasTraitObj.publishedIn.add(pub);
				pub.traits.add(gwasTraitObj);
				
			}
			elems = tf.readLineElemsReturnReference(TextFile.tab);
		}
		
		System.out.println(numpubs + " pubs, " + numsnps + " snps, " + numtraits + " traits");
		tf.close();
	}
	
	public GWASTrait[] getTraits() {
		if (traitArray == null) {
			traitArray = new GWASTrait[traits.size()];
			traits.toArray(traitArray);
		}
		return traitArray;
	}
	
	/**
	 * @return the loci
	 */
	public HashSet<GWASLocus> getLoci() {
		return loci;
	}
	
	/**
	 * @param loci the loci to set
	 */
	public void setLoci(HashSet<GWASLocus> loci) {
		this.loci = loci;
	}
	
	/**
	 * @return the snps
	 */
	public HashSet<GWASSNP> getSnps() {
		return snps;
	}
	
	/**
	 * @param snps the snps to set
	 */
	public void setSnps(HashSet<GWASSNP> snps) {
		this.snps = snps;
	}
	
	/**
	 * @param traits the traits to set
	 */
	public void setTraits(HashSet<GWASTrait> traits) {
		this.traits = traits;
	}
	
	/**
	 * @return the publicationToObj
	 */
	public HashMap<String, GWASPublication> getPublicationToObj() {
		return publicationToObj;
	}
	
	/**
	 * @param publicationToObj the publicationToObj to set
	 */
	public void setPublicationToObj(HashMap<String, GWASPublication> publicationToObj) {
		this.publicationToObj = publicationToObj;
	}
	
	/**
	 * @return the snpToObj
	 */
	public HashMap<String, GWASSNP> getSnpToObj() {
		return snpToObj;
	}
	
	/**
	 * @param snpToObj the snpToObj to set
	 */
	public void setSnpToObj(HashMap<String, GWASSNP> snpToObj) {
		this.snpToObj = snpToObj;
	}
	
	/**
	 * @return the locusToObj
	 */
	public HashMap<String, GWASLocus> getLocusToObj() {
		return locusToObj;
	}
	
	/**
	 * @param locusToObj the locusToObj to set
	 */
	public void setLocusToObj(HashMap<String, GWASLocus> locusToObj) {
		this.locusToObj = locusToObj;
	}
	
	/**
	 * @return the traitToObj
	 */
	public HashMap<String, GWASTrait> getTraitToObj() {
		return traitToObj;
	}
	
	/**
	 * @param traitToObj the traitToObj to set
	 */
	public void setTraitToObj(HashMap<String, GWASTrait> traitToObj) {
		this.traitToObj = traitToObj;
	}
	
	public GWASSNP[] getSnpsArray() {
		GWASSNP[] snpsr = new GWASSNP[snps.size()];
		snpsr = snps.toArray(snpsr);
		return snpsr;
	}
	
	public GWASSNP[] getSNPsForTraitContainingKey(String key) {
		System.out.println("Looking for " + key + " snps");
		HashSet<GWASSNP> s = new HashSet<GWASSNP>();
		key = key.toLowerCase();
		for (GWASTrait t : traits) {
			if (t.getName().toLowerCase().contains(key)) {
				System.out.println("Found trait: " + t.getName());
				GWASSNP[] traitsnps = t.getSNPs();
				s.addAll(Arrays.asList(traitsnps));
			}
		}
		
		return s.toArray(new GWASSNP[s.size()]);
	}
	
	public GWASTrait[] getTraitsForCertainKey(String key) {
		key = key.toLowerCase();
		ArrayList<GWASTrait> selected = new ArrayList<GWASTrait>();
		for (GWASTrait t : traits) {
			if (t.getName().toLowerCase().contains(key)) {
				selected.add(t);
			}
		}
		
		return selected.toArray(new GWASTrait[selected.size()]);
	}
	
	public GWASLocus[] getLociForCertainKey(String key) {
		System.out.println("Looking for " + key + " snps");
		HashSet<GWASLocus> s = new HashSet<GWASLocus>();
		key = key.toLowerCase();
		for (GWASTrait t : traits) {
			if (t.getName().toLowerCase().contains(key)) {
				System.out.println("Found trait: " + t.getName());
				s.addAll(t.loci);
			}
		}
		
		return s.toArray(new GWASLocus[s.size()]);
	}
	
	public HashSet<String> getReportedGenesForCertainKey(String key) {
		System.out.println("Looking for " + key + " snps");
		HashSet<String> s = new HashSet<String>();
		key = key.toLowerCase();
		for (GWASTrait t : traits) {
			if (t.getName().toLowerCase().contains(key)) {
				System.out.println("Found trait: " + t.getName());
				s.addAll(t.getReportedGenes());
			}
		}
		
		return s;
	}
	
	public HashSet<String> getMappedGenesForCertainKey(String key) {
		System.out.println("Looking for " + key + " snps");
		HashSet<String> s = new HashSet<String>();
		key = key.toLowerCase();
		for (GWASTrait t : traits) {
			if (t.getName().toLowerCase().contains(key)) {
				System.out.println("Found trait: " + t.getName());
				s.addAll(t.getMappedGenes());
			}
		}
		
		return s;
	}
	
	public HashSet<String> getTraitsForCertainSnps(String key) {
		HashSet<String> m = new HashSet<String>();
		key = key.toLowerCase();
		for (GWASSNP s : snps) {
			if (s.getName().equalsIgnoreCase(key)) {
				HashSet<GWASTrait> t = s.getAssociatedTraits();
				for (GWASTrait tmp : t) {
					m.add(tmp.cleanName);
				}
			}
		}
		
		return m;
	}
	
	public HashSet<GWASTrait> getFullTraitsForCertainSnps(String key) {
		HashSet<GWASTrait> m = new HashSet<GWASTrait>();
		key = key.toLowerCase();
		for (GWASSNP s : snps) {
			if (s.getName().equalsIgnoreCase(key)) {
				m.addAll(s.getAssociatedTraits());
			}
		}
		
		return m;
	}
}


/*
Date Added to Catalog
PUBMEDID
First Author
Date
Journal
Link
Study
Disease/Trait
Initial Sample Size
Replication Sample Size
Region
Chr_id
Chr_pos
Reported Gene(s)
Mapped_gene
Upstream_gene_id
Downstream_gene_id
Snp_gene_ids
Upstream_gene_distance
Downstream_gene_distance
Strongest SNP-Risk Allele
SNPs
Merged
Snp_id_current
Context
Intergenic
Risk Allele Frequency
p-Value
Pvalue_mlog
p-Value (text)
OR or beta
95% CI (text)
Platform [SNPs passing QC]
CNV
 */
