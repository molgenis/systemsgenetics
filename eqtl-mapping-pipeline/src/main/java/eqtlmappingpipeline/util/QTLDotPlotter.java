/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package eqtlmappingpipeline.util;

import com.lowagie.text.DocumentException;
import eqtlmappingpipeline.metaqtl3.graphics.EQTLDotPlot;
import java.io.IOException;
import java.util.logging.Level;
import java.util.logging.Logger;

/**
 *
 * @author harmjan
 */
public class QTLDotPlotter {

    public void plot(String file) {
        System.out.println("Reading from: " + file);
        System.out.println("Writing to: " + file + "EQTLDotPlot.pdf");
        try {
            EQTLDotPlot edp = new EQTLDotPlot();
            edp.draw(file, file + "EQTLDotPlot.pdf", EQTLDotPlot.Output.PDF);
        } catch (IOException ex) {
            Logger.getLogger(QTLDotPlotter.class.getName()).log(Level.SEVERE, null, ex);
        } catch (DocumentException ex) {
            Logger.getLogger(QTLDotPlotter.class.getName()).log(Level.SEVERE, null, ex);
        }


    }
}
