package umcg.genetica.graphics;

import com.itextpdf.text.DocumentException;
import umcg.genetica.enums.Chromosome;
import umcg.genetica.features.Feature;
import umcg.genetica.features.FeatureComparator;
import umcg.genetica.features.GFFFile;
import umcg.genetica.graphics.themes.ComplementaryColor;
import umcg.genetica.graphics.themes.DefaultTheme;
import umcg.genetica.io.text.TextFile;


import java.awt.*;
import java.awt.geom.Arc2D;
import java.io.FileNotFoundException;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Collections;
import java.util.HashSet;

/**
 * Created by hwestra on 6/6/15.
 */
public class ChromosomePlot extends DefaultGraphics {

	private int margin;

	public ChromosomePlot(String name, int pageWidth, int pageHeight) throws FileNotFoundException, DocumentException {
		super(name, pageWidth, pageHeight);
		figureWidth = pageWidth;
	}

	public void setMargin(int margin) {
		this.margin = margin;
	}

	private ArrayList<Feature> getCytobands(String cytobandfile, Chromosome chr) throws IOException {

		ArrayList<Feature> output = new ArrayList<Feature>();
		TextFile tf = new TextFile(cytobandfile, TextFile.R);
		String[] elems = tf.readLineElems(TextFile.tab);
		while (elems != null) {

			Chromosome c = Chromosome.parseChr(elems[0]);
			if (c.equals(chr)) {
				Feature f = new Feature();
				f.setStart(Integer.parseInt(elems[1]));
				f.setStop(Integer.parseInt(elems[2]));
				f.setName(elems[4]);
				output.add(f);
			}

			elems = tf.readLineElems(TextFile.tab);
		}
		tf.close();
		return output;
	}

	private ArrayList<Feature> readRegionFile(String f) throws IOException {
		TextFile tf = new TextFile(f, TextFile.R);

		System.out.println("reading regions: " + f);
		ArrayList<Feature> output = new ArrayList<Feature>();
		String[] elems = tf.readLineElems(TextFile.tab);
		while (elems != null) {
			if (elems.length > 2) {
				Chromosome chr = Chromosome.parseChr(elems[0]);
				Feature feat = new Feature();
				feat.setChromosome(chr);
				feat.setStart(Integer.parseInt(elems[1]));
				feat.setStop(Integer.parseInt(elems[2]));
				output.add(feat);
			}
			elems = tf.readLineElems(TextFile.tab);
		}

		tf.close();
		return output;
	}

	public void plot(String cytobandfile, String[] gffLocusFiles, boolean onlySuggestedLoci, String sequencedRegionFile) throws IOException {

		// plot the chromosomes
		int maxChrHeight = 400;
		int chrWidth = 30;
		int nrRows = 2;

		Chromosome[] allChr = Chromosome.values();

		int ctr = 0;
		int maxSize = 0;
		int nrChr = 0;

		int betweenChrMarginX = 250;
		int betweenChrMarginY = 300;

		int fontsize = 15;

		for (Chromosome chr : allChr) {
			if (chr.getLengthB37() > 1) {
				nrChr++;
			}
			if (chr.getLengthB37() > maxSize) {
				maxSize = chr.getLengthB37();

			}
		}

		ArrayList<Feature> sequencedRegions = null;
		if (sequencedRegionFile != null) {
			sequencedRegions = readRegionFile(sequencedRegionFile);
		}

		ArrayList<ArrayList<Feature>> associatedLociPerDisease = new ArrayList<ArrayList<Feature>>();
		GFFFile gff = new GFFFile();
		for (int f = 0; f < gffLocusFiles.length; f++) {
			String gfffile = gffLocusFiles[f];
			associatedLociPerDisease.add(gff.readGFF(gfffile, onlySuggestedLoci));
		}


		int nrPerRow = nrChr / nrRows;

		Color[] colors = new Color[]{new Color(70, 67, 58), new Color(174, 164, 140)};


		g2d.setFont(new Font("Helvetica", Font.PLAIN, fontsize));


		double pixelsPerBp = (double) maxChrHeight / maxSize;
		System.out.println(pixelsPerBp);
		for (Chromosome chr : allChr) {
			// get cytobands


			if (chr.getLengthB37() > 1) {


				int row = (int) Math.floor(ctr / nrPerRow);
				int colNr = ctr % nrPerRow;

				// get the centromere
				ArrayList<Feature> allBands = getCytobands(cytobandfile, chr);
				// draw a circle
				int starty = margin + (maxChrHeight * row) + (betweenChrMarginY * row);
				int startx = margin + (betweenChrMarginX * colNr) + (chrWidth * colNr);
				int chrHeight = (int) Math.ceil(pixelsPerBp * chr.getLengthB37());
				System.out.println(chr.getName() + "\t" + row + "\t" + colNr + "\t" + startx + "\t" + starty + "\t" + chrHeight);

				boolean centromereseen = false;

				int bottomcentromerestop = 0;
				for (Feature f : allBands) {
					int bandStart = starty + (int) Math.ceil(f.getStart() * pixelsPerBp);
					int bandStop = starty + (int) Math.ceil(f.getStop() * pixelsPerBp);
					String name = f.getName();
					boolean draw = false;

					if (name.equals("gneg")) {

					} else if (name.equals("gpos25")) {
						g2d.setColor(new Color(70, 67, 58, 100));
						draw = true;
					} else if (name.equals("gpos50")) {
						g2d.setColor(new Color(70, 67, 58, 75));
						draw = true;
					} else if (name.equals("gpos75")) {
						g2d.setColor(new Color(70, 67, 58, 50));
						draw = true;
					} else if (name.equals("gpos100")) {

						g2d.setColor(new Color(70, 67, 58, 25));
						draw = true;
					} else if (name.equals("acen")) {
						int[] xpoints = new int[3];
						int[] ypoints = new int[3];
						g2d.setColor(colors[0]);
						xpoints[0] = startx;
						xpoints[1] = startx + chrWidth / 2;
						xpoints[2] = startx + chrWidth;

						if (!centromereseen) {
							g2d.drawLine(startx, starty, startx, bandStart);
							g2d.drawLine(startx + chrWidth, starty, startx + chrWidth, bandStart);
							ypoints[0] = bandStart;
							ypoints[1] = bandStop;
							ypoints[2] = bandStart;
							centromereseen = true;
						} else {
							bottomcentromerestop = bandStop;
							ypoints[0] = bandStop;
							ypoints[1] = bandStart;
							ypoints[2] = bandStop;
						}
						Polygon p = new Polygon(xpoints, ypoints, 3);

						g2d.fillPolygon(p);

					} else {
						System.out.println("unknown dinges: " + name);
					}

					if (draw) {

//					System.out.println(f.getStart() + "\t" + f.getStop() + "\t" + bandStart + "\t" + bandStop);
						g2d.fillRect(startx, bandStart, chrWidth, (bandStop - bandStart));

					}


				}
				ctr++;


				g2d.setColor(colors[1]);
				g2d.setFont(new Font("Helvetica", Font.BOLD, fontsize));
				FontMetrics metrics = g2d.getFontMetrics();
				g2d.drawString(chr.getName(), startx + (chrWidth / 2) - (metrics.stringWidth(chr.getName()) / 2), starty - 40);

				g2d.setColor(colors[0]);
				g2d.setFont(new Font("Helvetica", Font.PLAIN, fontsize));

				g2d.drawLine(startx, bottomcentromerestop, startx, bottomcentromerestop + (chrHeight - (bottomcentromerestop - starty)));
				g2d.drawLine(startx + chrWidth, bottomcentromerestop, startx + chrWidth, bottomcentromerestop + (chrHeight - (bottomcentromerestop - starty)));

				g2d.draw(new Arc2D.Double(startx, starty - (chrWidth / 2), chrWidth, chrWidth, 0, 180, Arc2D.OPEN));

				g2d.draw(new Arc2D.Double(startx, starty + chrHeight - (chrWidth / 2), chrWidth, chrWidth, 180, 180, Arc2D.OPEN));


				// print the loci...
				g2d.setColor(colors[0]);

				// this only works for two datasets..
				if (associatedLociPerDisease.size() <= 2) {
					for (int dataset = 0; dataset < associatedLociPerDisease.size(); dataset++) {

						ArrayList<Feature> lociOnChr = new ArrayList<Feature>();
						ArrayList<Feature> lociForDs = associatedLociPerDisease.get(dataset);

						for (Feature feature : lociForDs) {

							if (feature.getChromosome() == null) {
								System.err.println(feature.getName() + " null chr");
							}
							if (feature.getName().equals("ACOXL")) {
								System.out.println("found it");
							}
							if (feature.getChromosome().equals(chr)) {
								if (sequencedRegions != null) {
									boolean includeregion = false;
									for (Feature seqregion : sequencedRegions) {
										if (seqregion.getChromosome().equals(chr)) {
											if (feature.overlaps(seqregion)) {
												includeregion = true;
											}
										}
									}
									if (includeregion) {
										lociOnChr.add(feature);
									}
								} else {
									lociOnChr.add(feature);
								}
							}
						}

						// assign pixel locations
						for (int i = 0; i < lociOnChr.size(); i++) {
							Feature feature = lociOnChr.get(i);
							// starty + (int) Math.ceil(f.getStart() * pixelsPerBp)

							int geneY1 = starty + (int) Math.ceil(feature.getStart() * pixelsPerBp);
							System.out.println(feature.getName() + "\t" + feature.getStart() + "\t" + geneY1);
							int geneY2 = starty + (int) Math.ceil(feature.getStop() * pixelsPerBp);
							feature.setStart(geneY1);
							feature.setStop(geneY2);

						}


						System.out.println(lociOnChr.size() + " loci on chr " + chr.getName());
						Collections.sort(lociOnChr, new FeatureComparator(false));

						// plot the gene names
						int lastGeneY = starty;
						int geneMarginY = 10;

						metrics = g2d.getFontMetrics();
						int fontheight = metrics.getHeight();
						for (int i = 0; i < lociOnChr.size(); i++) {
							Feature f2 = lociOnChr.get(i);
							int start = f2.getStart();
							System.out.println(f2.getName() + "\t" + f2.getStart() + "\t" + start);
							if (start <= lastGeneY + geneMarginY + fontheight) {
								start = lastGeneY + geneMarginY + fontheight;
							}

							lastGeneY = start;
							if (dataset == 0) {
								int geneStartX = startx + chrWidth + chrWidth + 10;
								System.out.println(f2.getStart() + "\t" + start);
								g2d.drawLine(startx + chrWidth + 5, f2.getStart() - 5, startx + chrWidth + 30, start - 5);

								g2d.drawString(f2.getName(), geneStartX, start);
							} else if (dataset == 1) {
								int geneStartX = startx - chrWidth - 10;
								System.out.println(f2.getStart() + "\t" + start);
								int geneStrwidth = metrics.stringWidth(f2.getName());
								g2d.drawLine(startx - 5, f2.getStart() - 5, startx - 30, start - 5);

								g2d.drawString(f2.getName(), geneStartX - geneStrwidth, start);
							}
						}
					}

					DefaultTheme theme = new DefaultTheme();

					if (sequencedRegions != null) {
						// plot the regions as well
						for (Feature feature : sequencedRegions) {
							if (feature.getChromosome().equals(chr)) {
								starty = margin + (maxChrHeight * row) + (betweenChrMarginY * row);
								int geneY1 = starty + (int) Math.ceil(feature.getStart() * pixelsPerBp);
								int geneY2 = starty + (int) Math.ceil(feature.getStop() * pixelsPerBp);
								int height = (geneY2 - geneY1);
								if (height < 2) {
									height = 2;
								}

								g2d.setColor(theme.getColor(1));
								g2d.fillRect(startx, geneY1 - 5, chrWidth, height);
							}
						}
					}
				} else {

					HashSet<Feature> lociOnChrHash = new HashSet<Feature>();
					ArrayList<HashSet<Feature>> lociOnChrPerDisease = new ArrayList<HashSet<Feature>>();
					for (int dataset = 0; dataset < associatedLociPerDisease.size(); dataset++) {
						ArrayList<Feature> f = associatedLociPerDisease.get(dataset);
						HashSet<Feature> f2 = new HashSet<Feature>();
						for (int i = 0; i < f.size(); i++) {
							Feature feat = f.get(i);
							if (feat.getChromosome().equals(chr)) {
								lociOnChrHash.add(f.get(i));
								f2.add(f.get(i));
							}
						}
						lociOnChrPerDisease.add(f2);
					}

					ArrayList<Feature> lociOnChr = new ArrayList<Feature>();
					lociOnChr.addAll(lociOnChrHash);
					System.out.println(lociOnChr.size() + " unique loci on chr");

					// assign pixel locations
//					ArrayList<Feature> lociOnChrPixelPos = new ArrayList<Feature>();
//					for (int i = 0; i < lociOnChr.size(); i++) {
//						Feature feature = lociOnChr.get(i);
//						// starty + (int) Math.ceil(f.getStart() * pixelsPerBp)
//						Feature copy = new Feature();
//						copy.setStrand(feature.getStrand());
//						copy.setChromosome(feature.getChromosome());
//
//
//						int geneY1 = starty + (int) Math.ceil(feature.getStart() * pixelsPerBp);
//						// System.out.println(feature.getName() + "\t" + feature.getStart() + "\t" + geneY1);
//						int geneY2 = starty + (int) Math.ceil(feature.getStop() * pixelsPerBp);
//						copy.setStart(geneY1);
//						copy.setStop(geneY2);
//						copy.setName(feature.getName());
//						lociOnChrPixelPos.add(copy);
//
//					}

					System.out.println(lociOnChr.size() + " loci on chr " + chr.getName());
					Collections.sort(lociOnChr, new FeatureComparator(false));


					// plot the gene names
					int lastGeneY = starty;
					int geneMarginY = 10;

					metrics = g2d.getFontMetrics();
					int fontheight = metrics.getHeight();

					ComplementaryColor complement = new ComplementaryColor();

					for (int i = 0; i < lociOnChr.size(); i++) {
						Feature feature = lociOnChr.get(i);

						int start = starty + (int) Math.ceil(feature.getStart() * pixelsPerBp);
						int origStart = start;
						//System.out.println(f2.getName() + "\t" + f2.getStart() + "\t" + start);
						if (start <= lastGeneY + geneMarginY + fontheight) {
							start = lastGeneY + geneMarginY + fontheight;
						}

						lastGeneY = start;
//						if (dataset == 0) {
						int geneStartX = startx + chrWidth + (chrWidth) + 10;
						System.out.println(origStart + "\t" + start);
						g2d.setColor(colors[1]);
						g2d.drawLine(startx + chrWidth + 5, origStart - 5, startx + chrWidth + 30, start - 5);

						int circlemargin = 5;
						int totalWidthOfDiseaseColumns = (lociOnChrPerDisease.size() * circlemargin) + (lociOnChrPerDisease.size() * fontsize);
						g2d.setColor(colors[0]);
						g2d.drawString(feature.getName(), geneStartX + totalWidthOfDiseaseColumns + 5, start);

						int circlestart = geneStartX;

						for (int disease = 0; disease < lociOnChrPerDisease.size(); disease++) {
							HashSet<Feature> f = lociOnChrPerDisease.get(disease);
							int circlex = circlestart + (disease * circlemargin) + (disease * fontsize);

							int circlewidth = fontheight;
							int circleheight = fontheight;
							int circley = start - circleheight;

							if (f.contains(feature)) {
								g2d.setColor(complement.getColor(disease));
								g2d.fillOval(circlex, circley, circlewidth, circleheight);
								g2d.setColor(complement.getDarkerColor(complement.getColor(disease), 0.2));
								g2d.drawOval(circlex, circley, circlewidth, circleheight);
							} else {
								g2d.setColor(complement.getDarkerColor(complement.getColor(disease), 0.2));
								g2d.drawOval(circlex, circley, circlewidth, circleheight);
							}
							g2d.setColor(colors[0]);
						}

//						} else if (dataset == 1) {
//						int geneStartX = startx - chrWidth - 10;
//						//System.out.println(f2.getStart() + "\t" + start);
//						int geneStrwidth = metrics.stringWidth(f2.getName());
//						g2d.drawLine(startx - 5, f2.getStart() - 5, startx - 30, start - 5);
//
//						g2d.drawString(f2.getName(), geneStartX - geneStrwidth, start);
//						}
					}


				}

			}
		}

		close();
	}

	public void heatmap(String[] gffLocusFiles, boolean onlySuggestedLoci) throws IOException {

		// plot the chromosomes
		int maxChrHeight = 1000;
		int chrWidth = 30;
		int nrRows = 1;

		Chromosome[] allChr = Chromosome.values();

		int ctr = 0;
		int maxSize = 0;
		int nrChr = 0;

		int betweenChrMarginX = 450;
		int betweenChrMarginY = 100;

		int fontsize = 15;

		for (Chromosome chr : allChr) {
			if (chr.getLengthB37() > 1) {
				nrChr++;
			}
			if (chr.getLengthB37() > maxSize) {
				maxSize = chr.getLengthB37();

			}
		}


		ArrayList<ArrayList<Feature>> associatedLociPerDisease = new ArrayList<ArrayList<Feature>>();
		HashSet<Feature> allLociHash = new HashSet<Feature>();

		GFFFile gff = new GFFFile();
		for (int f=0;f<gffLocusFiles.length;f++) {
			String gfffile = gffLocusFiles[f];
			ArrayList<Feature> lociForDisease = gff.readGFF(gfffile, onlySuggestedLoci);
			associatedLociPerDisease.add(lociForDisease);
			allLociHash.addAll(lociForDisease);
		}

		// sort
		ArrayList<Feature> allLoci = new ArrayList<Feature>();
		allLoci.addAll(allLociHash);
		Collections.sort(allLoci, new FeatureComparator(false));

		int circlemargin = 10;
		int circlesize = fontsize;

		int nrPerRow = nrChr / nrRows;
		DefaultTheme t = new DefaultTheme();
		Color[] colors = new Color[]{t.getColor(0), t.getColor(1)};


		g2d.setFont(new Font("Helvetica", Font.PLAIN, fontsize));


		double pixelsPerBp = (double) maxChrHeight / maxSize;
		System.out.println(pixelsPerBp);
		g2d.setFont(new Font("Helvetica", Font.BOLD, fontsize));
		FontMetrics metrics = g2d.getFontMetrics();


		for (int disease = 0; disease < associatedLociPerDisease.size(); disease++) {
			int starty = margin + (disease * circlemargin) + (disease * circlesize);
			int startx = margin;

			HashSet<Feature> diseaseLoci = new HashSet<Feature>();
			ArrayList<Feature> diseaseLociArr = associatedLociPerDisease.get(disease);
			diseaseLoci.addAll(diseaseLociArr);

			Chromosome prevChr = null;

			int q = 0;
			g2d.setColor(colors[0]);
//			for (int f = 0; f < allLoci.size(); f++) {
//				int circleX = startx + (f * circlesize) + (f * circlemargin);
//
//				Feature locus = allLoci.get(f);
//
//				if (prevChr == null) {
//					prevChr = locus.getChromosome();
//				} else {
//					if (prevChr.equals(locus.getChromosome())) {
//						// don't change color
//					} else {
//						q++;
//						int remainder = q % 2;
//						g2d.setColor(colors[remainder]);
//						prevChr = locus.getChromosome();
//					}
//				}
//
//				if (diseaseLoci.contains(locus)) {
//					g2d.fillOval(circleX, starty, circlesize, circlesize);
//				} else {
//					g2d.drawOval(circleX, starty, circlesize, circlesize);
//				}
//
//			}


		}

		int nrDiseases = gffLocusFiles.length;
		int totalWidthOfCircles = (nrDiseases * circlemargin) + (nrDiseases + circlesize);


		for (Chromosome chr : allChr) {
			// get cytobands
			if (chr.getLengthB37() > 1) {
				ArrayList<Feature> lociOnChr = new ArrayList<Feature>();
				for (int loc = 0; loc < allLoci.size(); loc++) {
					Feature f = allLoci.get(loc);
					if (f.getChromosome().equals(chr)) {
						lociOnChr.add(f);
					}
				}

				int chrNr = chr.getNumber();
				int startX = margin + (chrNr * betweenChrMarginX) + (chrNr * totalWidthOfCircles);

				int startY = margin;

				for (int disease = 0; disease < nrDiseases; disease++) {
					HashSet<Feature> diseaseLoci = new HashSet<Feature>();
					ArrayList<Feature> diseaseLociArr = associatedLociPerDisease.get(disease);
					diseaseLoci.addAll(diseaseLociArr);

					for (int locus = 0; locus < lociOnChr.size(); locus++) {


						int locusXpos = startX + (disease * circlemargin) + (disease * circlesize);
						int locusYPos = startY + (locus * circlemargin) + (locus * circlesize) - circlesize;
						Feature loc = lociOnChr.get(locus);
						if (diseaseLoci.contains(loc)) {
							g2d.fillOval(locusXpos, locusYPos, circlesize, circlesize);
						} else {
							g2d.drawOval(locusXpos, locusYPos, circlesize, circlesize);
						}
					}

				}

				// draw the locus names
				for (int locus = 0; locus < lociOnChr.size(); locus++) {
					Feature loc = lociOnChr.get(locus);
					String locusName = loc.getName();
					int locusXpos = startX - fontsize - metrics.stringWidth(locusName);
					int locusYPos = startY + (locus * circlemargin) + (locus * circlesize);

					g2d.drawString(locusName, locusXpos, locusYPos);


				}


			}
		}

		close();
	}

}
