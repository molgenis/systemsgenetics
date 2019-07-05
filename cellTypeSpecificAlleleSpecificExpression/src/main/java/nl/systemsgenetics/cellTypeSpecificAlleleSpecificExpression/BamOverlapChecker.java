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

import htsjdk.samtools.SAMFileHeader;
import htsjdk.samtools.SAMRecordIterator;
import htsjdk.samtools.SAMSequenceDictionary;
import htsjdk.samtools.SAMSequenceRecord;
import htsjdk.samtools.SamReader;
import java.util.HashMap;
import java.util.List;

/**
 *
 * @author adriaan
 */
class BamOverlapChecker {
    
    
    
    int stepSize = 100000;
    HashMap<String, boolean[]> booleanMap;
    
    public BamOverlapChecker(SamReader bam_file){
        
        SAMFileHeader  header = bam_file.getFileHeader();
        SAMSequenceDictionary dict = header.getSequenceDictionary();
        List<SAMSequenceRecord> sequences = dict.getSequences();
       
        
        booleanMap = new HashMap<String, boolean[]>();
        
        for(SAMSequenceRecord sequence : sequences){
            int sequenceEnd = sequence.getSequenceLength();
            int arrayLength = (int) Math.ceil( (float) sequenceEnd / (float)stepSize );
            boolean[] tempArray;
            tempArray = new boolean[arrayLength];
            
            for(int i=0;i<arrayLength;i++){
                SAMRecordIterator bamQuery = bam_file.queryOverlapping(sequence.getSequenceName(), i*stepSize, (i+1)*stepSize);
                if(bamQuery.hasNext()){
                    tempArray[i] = true;
                }else{
                    tempArray[i] = false;
                }
                bamQuery.close();
            }
            booleanMap.put(sequence.getSequenceName(),tempArray );
            //System.out.println("Finished checking the bam for chromosome " + sequence.getSequenceName());
        }
        
    }
    
    public boolean bamHasOverlap(String sequenceName, int position){
        int indice = (int) Math.floor((double) position / (double) stepSize);
        return booleanMap.get(sequenceName)[indice];
    }
    
}
