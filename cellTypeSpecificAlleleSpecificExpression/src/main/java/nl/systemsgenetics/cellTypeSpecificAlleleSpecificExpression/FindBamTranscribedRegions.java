/*
 * Copyright (C) 2016 adriaan
 *
 * This program is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program.  If not, see <http://www.gnu.org/licenses/>.
 */
package nl.systemsgenetics.cellTypeSpecificAlleleSpecificExpression;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.HashMap;
import org.jdom.IllegalDataException;

/**
 *
 * @author adriaan
 * 
 * This method will read in a .bai index (as an array of bytes)
 * and it will return a hashmap containing a boolean per binlocation.
 * 
 */
public class FindBamTranscribedRegions {
    
    
    public int num_reference_sequences;
    
    public HashMap<String, Boolean> bin_filled = new HashMap<String, Boolean>();
    public ArrayList<boolean[]> referenceList = new ArrayList<boolean[]>();
    
    
    
    public FindBamTranscribedRegions(byte[] tabix_bytes) throws IllegalDataException{
        
        int iArray = 0;
        String magic;
        
        if( (tabix_bytes[0] != 'B') || 
            (tabix_bytes[1] != 'A') || 
            (tabix_bytes[2] != 'I') || 
            (tabix_bytes[3] != '\1')){
            throw new IllegalDataException("The tabix magic string did not start with the expected BAI\1 string");
        }
        iArray +=4;
        
        num_reference_sequences = decode32IntFromArray(tabix_bytes, iArray);
        iArray +=4;
        
        //Do for every chromosome (reference)
        for(int refs=0; refs < num_reference_sequences; refs++){
            
            int empty=0;
            int full=0;
            
            
            int nBinsThisRef = decode32IntFromArray(tabix_bytes, iArray);
            iArray+=4;
            
            boolean[] tempList = new boolean[37451];
            for(int n=0; n < 37451; n++){
                tempList[n] = false;
            }
            
            //do for every bin in the chromosome (reference)
            for(int refBins=0; refBins < nBinsThisRef; refBins++ ){
                //no need to do anything for the bin 37450, as it will be processed correctly.
               
                int bin = (int) decode32uIntFromArray(tabix_bytes, iArray);
                iArray +=4;
                
                int n_chunks = decode32IntFromArray(tabix_bytes, iArray);
                iArray +=4;
                
                byte[] startValue = Arrays.copyOfRange(tabix_bytes, iArray, iArray+8);
                byte[] endValue   = Arrays.copyOfRange(tabix_bytes, iArray + 16*(n_chunks-1), iArray + 16*(n_chunks-1)+ 8);
                iArray += 16*n_chunks;
                
                //if the bin is empty, show it.
                if(Arrays.equals(startValue,endValue)){
                    tempList[bin] = true;
                    empty++;
                }else{
                    tempList[bin] = false;
                    full++;
                }
                
            }
            
            int num_intervals = decode32IntFromArray(tabix_bytes, iArray);
            iArray+=4;
            iArray+= num_intervals*8;
            
            referenceList.add(tempList);
            
            //System.out.printf("If I'm correct. Sequence %d contains %d empty and %d full bins.\n", refs, empty, full );
        }
        
        
       
        
    }
            
    
    //adapted from SO: stackoverflow.com/questions/23417810/encoding-int32-t-to-a-byte-array
    public static int decode32IntFromArray(byte[] array, int pos)
    {
        return (int)
                (
                        ((array[pos+3] & 0xFF) << 24) 
                    +   ((array[pos+2] & 0xFF) << 16)                 
                    +   ((array[pos+1] & 0xFF) << 8) 
                    +   ((array[pos+0] & 0xFF))
                );        
    }
    
    public static long decode32uIntFromArray(byte[] array, int pos)
    {
        return (long)
                (
                        ((array[pos+3] & 0xFF) << 24) 
                    +   ((array[pos+2] & 0xFF) << 16)                 
                    +   ((array[pos+1] & 0xFF) << 8) 
                    +   ((array[pos+0] & 0xFF))
                );        
    }
    

    //taken from the sam specification sheet.
    private static int reg2bin(int begin, int end)
    {
        --end;
        if (begin>>14 == end>>14) return ((1<<15)-1)/7 + (begin>>14);
        if (begin>>17 == end>>17) return ((1<<12)-1)/7 + (begin>>17);
        if (begin>>20 == end>>20) return ((1<<9)-1)/7 + (begin>>20);
        if (begin>>23 == end>>23) return ((1<<6)-1)/7 + (begin>>23);
        if (begin>>26 == end>>26) return ((1<<3)-1)/7 + (begin>>26);
        return 0;
    }
    
    public boolean bamHasOverlap(int reference, int bpPosition){
        return !(referenceList.get(reference)[reg2bin(bpPosition-1, bpPosition+1)]);
    }
    
}
