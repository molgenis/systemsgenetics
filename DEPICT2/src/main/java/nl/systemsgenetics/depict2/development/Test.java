/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package nl.systemsgenetics.depict2.development;

import java.io.FileNotFoundException;
import java.io.IOException;
import org.apache.commons.math3.stat.regression.SimpleRegression;

/**
 *
 * @author patri
 */
public class Test {

	/**
	 * @param args the command line arguments
	 */
	public static void main(String[] args) throws IOException, FileNotFoundException {

		//System.out.println(DurationFormatUtils.formatDuration(1123456789, "H:mm:ss.S"));
//		ChiSquaredDistribution x = new ChiSquaredDistribution(1);
//		double y = x.inverseCumulativeProbability(0.5);
//		System.out.println(y);
		//System.out.println(new ChiSquaredDistribution(1).inverseCumulativeProbability(0.5));
//		Workbook enrichmentWorkbook = new XSSFWorkbook();
//
//		XSSFSheet databaseSheet = (XSSFSheet) enrichmentWorkbook.createSheet("PathwayDatabase");
//
//		XSSFTable table = databaseSheet.createTable(new AreaReference(new CellReference(0, 0), new CellReference(3, 2), SpreadsheetVersion.EXCEL2007));
//		table.setName("database1");
//		table.setDisplayName("database1");
//
//		table.setStyleName("TableStyleLight2");
//		table.getCTTable().getTableStyleInfo().setShowRowStripes(true);
//
//		table.getCTTable().addNewAutoFilter();
//
//		databaseSheet.createFreezePane(0, 1);
//		table_style.getStyle()
//		table_style.setShowColumnStripes(false); //showColumnStripes=0
//		table_style.setShowRowStripes(true); //showRowStripes=1   
//		XSSFRow headerRow = databaseSheet.createRow(0);
//		headerRow.createCell(0, CellType.STRING).setCellValue("LongHeader1");
//		headerRow.createCell(1, CellType.STRING).setCellValue("Header2");
//		headerRow.createCell(2, CellType.STRING).setCellValue("Header3");
//
//		for (int r = 1; r < 4; ++r) {
//			XSSFRow row = databaseSheet.createRow(r);
//			for (int c = 0; c < 3; ++c) {
//				row.createCell(c, CellType.NUMERIC).setCellValue(r * c);
//			}
//		}
//
//		for (int c = 0; c < 3; ++c) {
//			databaseSheet.autoSizeColumn(c);
//			System.out.println(databaseSheet.getColumnWidth(c));
//			databaseSheet.setColumnWidth(c, databaseSheet.getColumnWidth(c) + 500);//compensate for with auto filter
//		}
//
//		enrichmentWorkbook.write(new FileOutputStream(new File("C:\\UMCG\\Genetica\\Projects\\Depict2Pgs\\test.xlsx")));
		SimpleRegression regression = new SimpleRegression();

		regression.addData(10, 3);
		regression.addData(11, 4);
		regression.addData(12, 5);
		regression.addData(11, 4);
		regression.addData(9, 3);
		regression.addData(9, 2);

		System.out.println("Intercept: " + regression.getIntercept());
		System.out.println("Slope: " + regression.getSlope());
		
		double b0 = regression.getIntercept();
		double b1 = regression.getSlope();
		
		System.out.println(regression.predict(10));
		System.out.println(b0 + (b1 * 10));
		System.out.println(3 - regression.predict(10));
		

		System.out.println(regression.predict(11));
		System.out.println(4 - regression.predict(11));

	}

}
