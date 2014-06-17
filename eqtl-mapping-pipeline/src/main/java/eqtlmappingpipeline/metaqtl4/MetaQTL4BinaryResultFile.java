/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package eqtlmappingpipeline.metaqtl4;

import java.io.IOException;
import umcg.genetica.io.bin.BinaryFile;

/**
 *
 * @author Harm-Jan
 */
public class MetaQTL4BinaryResultFile extends BinaryFile {

    public MetaQTL4BinaryResultFile(String loc, boolean mode) throws IOException {
        super(loc, mode);
    }

    public void readPermutationFile(double[] buffer) throws IOException {
        int bufferLen = buffer.length;
        for (int i = 0; i < bufferLen; i++) {
            buffer[i] = this.readDouble();
        }
    }

    public void writePermutationFile(double[] buffer) throws IOException{
        int bufferLen = buffer.length;
        for(int i=0; i<bufferLen; i++){
            this.writeDouble(buffer[i]);
        }
    }
    
    
}
