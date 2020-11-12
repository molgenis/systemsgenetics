/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package nl.systemsgenetics.downstreamer.io;

import org.apache.poi.ss.usermodel.CellStyle;
import org.apache.poi.ss.usermodel.DataFormat;
import org.apache.poi.ss.usermodel.Font;
import org.apache.poi.ss.usermodel.HorizontalAlignment;
import org.apache.poi.ss.usermodel.Workbook;

/**
 *
 * @author patri
 */
public class ExcelStyles {

	private final CellStyle zscoreStyle;
	private final CellStyle largePvalueStyle;
	private final CellStyle smallPvalueStyle;
	private final CellStyle hlinkStyle;
	private final CellStyle boldStyle;
	private final CellStyle genomicPositionStyle;
	private final CellStyle boldGenomicPositionStyle;
	private final CellStyle rightAlignedText;

	public ExcelStyles(Workbook wb) {

		DataFormat format = wb.createDataFormat();

		//Also used for OR and AUC
		zscoreStyle = wb.createCellStyle();
		zscoreStyle.setDataFormat(format.getFormat("0.00"));

		largePvalueStyle = wb.createCellStyle();
		largePvalueStyle.setDataFormat(format.getFormat("0.0000"));

		smallPvalueStyle = wb.createCellStyle();
		smallPvalueStyle.setDataFormat(format.getFormat("0.00E+0"));

		hlinkStyle = wb.createCellStyle();
		Font hlinkFont = wb.createFont();
		hlinkFont.setUnderline(Font.U_SINGLE);
		hlinkStyle.setFont(hlinkFont);

		boldStyle = wb.createCellStyle();
		Font fontBold = wb.createFont();
		fontBold.setBold(true);
		boldStyle.setFont(fontBold);

		genomicPositionStyle = wb.createCellStyle();
		genomicPositionStyle.setDataFormat(format.getFormat("###,###,##0"));

		boldGenomicPositionStyle = wb.createCellStyle();
		fontBold.setFontHeightInPoints((short) 10);
		boldGenomicPositionStyle.setFont(fontBold);
		boldGenomicPositionStyle.setDataFormat(format.getFormat("###,###,##0"));

		rightAlignedText = wb.createCellStyle();
		rightAlignedText.setAlignment(HorizontalAlignment.RIGHT);

	}

	/**
	 * Also used for OR and AUC
	 * 
	 * @return 
	 */
	public CellStyle getZscoreStyle() {
		return zscoreStyle;
	}

	public CellStyle getLargePvalueStyle() {
		return largePvalueStyle;
	}

	public CellStyle getSmallPvalueStyle() {
		return smallPvalueStyle;
	}

	public CellStyle getHlinkStyle() {
		return hlinkStyle;
	}

	public CellStyle getBoldStyle() {
		return boldStyle;
	}

	public CellStyle getGenomicPositionStyle() {
		return genomicPositionStyle;
	}

	public CellStyle getBoldGenomicPositionStyle() {
		return boldGenomicPositionStyle;
	}

	public CellStyle getRightAlignedText() {
		return rightAlignedText;
	}
	
	

}
