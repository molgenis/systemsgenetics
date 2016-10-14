/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package eqtlmappingpipeline.chromosomeyexpressionplotter;

import eqtlmappingpipeline.metaqtl3.MetaQTL3;
import eqtlmappingpipeline.metaqtl3.containers.Settings;
import java.io.IOException;
import umcg.genetica.io.Gpio;

/**
 *
 * @author harmjan
 */
public class ChromosomeYExpressionPlotter extends MetaQTL3 {

    public void run(Settings settings) {
        try {
            m_settings = settings;
            String origOutputDir = m_settings.outputReportsDir;
            Gpio.createDir(m_settings.outputReportsDir);
            initialize(null, null, null, null, null, null, null, null, null, null, null, true, true, 0, true, false, null, null, null, null, null, true, true, null, null, null);
            ChromosomeYExpressionPlot cp = new ChromosomeYExpressionPlot();
            cp.draw(m_gg, origOutputDir);

        } catch (IOException e) {
            e.printStackTrace();
        } catch (Exception e) {
            e.printStackTrace();
        }
    }
}
