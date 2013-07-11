/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package umcg.genetica.io.trityper.bin;

import java.io.DataOutputStream;
import java.io.File;
import java.io.FileOutputStream;
import java.io.IOException;
import java.io.RandomAccessFile;
import java.nio.ByteBuffer;
import java.util.zip.DataFormatException;
import java.util.zip.Deflater;
import java.util.zip.Inflater;

/**
 *
 * @author harm-jan
 */
public class BinaryGZipFloatMatrix {

    private DataOutputStream out;

    public static boolean W = true;
    public static boolean R = false;

    private final Deflater d = new Deflater(9);
    private long indexposition = 0;
    private long iterations = 0;
    private long ratio = 0;
    private long uncompressedSize = 0;

    private final Inflater inflater = new Inflater();
    private RandomAccessFile ra;

    private byte[] writebuffer = null;
    private byte[] tmpbytebuffer = null;
    private int nrOfBytesInWriteBuffer = 0;
    private int writebuffersize = 134217728;
    private int tmpbuffersize = 4096;

    public BinaryGZipFloatMatrix(String filename, boolean W) throws IOException {
        if(W){
            out = new DataOutputStream(new FileOutputStream(new File(filename+".ZScoreMatrix.dat")));
            writebuffer = new byte[writebuffersize];
            tmpbytebuffer = new byte[tmpbuffersize];
        } else {
            ra = new RandomAccessFile(filename,"r");

        }
    }

    public long write(Float[] zscorelist) throws IOException {
        long positionToReturn = indexposition;
        int nrBytesRequired  = (zscorelist.length*4);
        ByteBuffer buff       = ByteBuffer.allocate( nrBytesRequired );

        for(Float f: zscorelist){
            if(f == null){
                System.out.println("ERROR!");
            }

            buff.putFloat(f);
        }

        byte[] input = buff.array();
        d.setInput(input);
        d.finish();

        int compressedDataLength = tmpbuffersize;
        long compressedsize = 0;
        while(compressedDataLength == tmpbuffersize){
            compressedDataLength = d.deflate(tmpbytebuffer);
            // out.write(bytebuffer, 0, compressedDataLength);
            write(tmpbytebuffer, 0, compressedDataLength);
            compressedsize      += compressedDataLength;
        }

        d.reset();

        uncompressedSize+= nrBytesRequired;
        indexposition   += compressedsize;
        ratio           += (long) Math.floor( (double) input.length / compressedsize) ;

        iterations++;
        return positionToReturn;
    }

    public long write(byte[] input) throws IOException {
        long positionToReturn = indexposition;
//        int nrBytesRequired  = (zscorelist.length*4);
//        ByteBuffer buff       = ByteBuffer.allocate( nrBytesRequired );
//
//        for(Float f: zscorelist){
//            if(f == null){
//                System.out.println("ERROR!");
//            }
//            buff.putFloat(f);
//        }
//
//        byte[] input = buff.array();
        d.setInput(input);
        d.finish();

        int compressedDataLength = tmpbuffersize;
        long compressedsize = 0;
        while(compressedDataLength == tmpbuffersize){
            compressedDataLength = d.deflate(tmpbytebuffer);
            // out.write(bytebuffer, 0, compressedDataLength);
            write(tmpbytebuffer, 0, compressedDataLength);
            compressedsize      += compressedDataLength;
        }

        d.reset();

        uncompressedSize+= input.length;
        indexposition   += compressedsize;
        ratio           += (long) Math.floor( (double) input.length / compressedsize) ;

        iterations++;
        return positionToReturn;
    }

    public long writeDeflated(byte[] input) throws IOException {
        long positionToReturn = indexposition;
//        int nrBytesRequired  = (zscorelist.length*4);
//        ByteBuffer buff       = ByteBuffer.allocate( nrBytesRequired );
//
//        for(Float f: zscorelist){
//            if(f == null){
//                System.out.println("ERROR!");
//            }
//            buff.putFloat(f);
//        }
//
//        byte[] input = buff.array();

	write(input, 0, input.length);
//	d.setInput(input);
//        d.finish();
//
//        int compressedDataLength = tmpbuffersize;
//        long compressedsize = 0;
//        while(compressedDataLength == tmpbuffersize){
//            compressedDataLength = d.deflate(tmpbytebuffer);
//            // out.write(bytebuffer, 0, compressedDataLength);
//            write(tmpbytebuffer, 0, compressedDataLength);
//            compressedsize      += compressedDataLength;
//        }
//
//        d.reset();
//
//        uncompressedSize+= input.length;
        indexposition   += input.length;
//        ratio           += (long) Math.floor( (double) input.length / compressedsize) ;

        iterations++;
//	System.out.println(indexposition+"\tbytes written");
        return positionToReturn;
    }

    private void write(byte[] bytebuffer, int start, int end) throws IOException {

        int bytesToAdd = end - start;

        if(nrOfBytesInWriteBuffer + bytesToAdd >= writebuffersize){
            int nrOfBytesThatWillNotFitInWriteBuffer  = (nrOfBytesInWriteBuffer+bytesToAdd) - writebuffersize;
            int nrOfBytesThatWillFitInWriteBuffer     = bytesToAdd - nrOfBytesThatWillNotFitInWriteBuffer;

            // arraycopy(Object src, int srcPos, Object dest, int destPos, int length)
//            System.out.println(bytesToAdd+"\t0\t"+writebufferpointer+"\t"+toBuffer);

            System.arraycopy(bytebuffer, 0, writebuffer, nrOfBytesInWriteBuffer, nrOfBytesThatWillFitInWriteBuffer);
            nrOfBytesInWriteBuffer = writebuffersize;
            flush();
            System.arraycopy(bytebuffer, nrOfBytesThatWillFitInWriteBuffer, writebuffer, 0, nrOfBytesThatWillNotFitInWriteBuffer);
            nrOfBytesInWriteBuffer = nrOfBytesThatWillNotFitInWriteBuffer;
        } else {
            System.arraycopy(bytebuffer, 0, writebuffer, nrOfBytesInWriteBuffer, bytesToAdd);
            nrOfBytesInWriteBuffer += bytesToAdd;
        }


    }

    public void flush() throws IOException {
//        System.out.println("flushing: "+writebufferpointer);
        out.write(writebuffer, 0, nrOfBytesInWriteBuffer);
        nrOfBytesInWriteBuffer = 0;
    }

    public Float[] read(long index, long nextIndex, int numElems) throws IOException, DataFormatException {

        if(index > ra.length()){
            throw new IOException("IO Error: SNP ID out of scope!");
        }

        long seekloc = index;

        int compressedsize = 0;

        // if this is the last element, read to the end
        if(nextIndex == -1){
            long curpos  = seekloc;
            long nextpos = ra.length();
            compressedsize = (int) (nextpos - curpos); // TODO: watch out for buffer overruns here --> size might be greater than int.maxvalue
        } else {
            long curpos = seekloc;
            long nextpos = nextIndex;
            compressedsize = (int) (nextpos - curpos); // TODO: watch out for buffer overruns here --> size might be greater than int.maxvalue
        }
        if(seekloc + compressedsize > ra.length()){
            throw new IOException("IO Error: buffer to be loaded is out of scope");
        }

        byte[] buffer = new byte[compressedsize];
        ra.seek(seekloc);
        ra.read(buffer);
        inflater.setInput(buffer);
        inflater.finished();
        byte[] decompressed = new byte[numElems*4];
        inflater.inflate(decompressed);

        long actuallydecompressed = inflater.getBytesWritten();
        if(actuallydecompressed != numElems*4){
            throw new IOException("IO Error: uncompressed data does not correspond to the size requested\t"+actuallydecompressed+"\t"+numElems*4+"\t Index: "+index);
        }

        inflater.reset();



        ByteBuffer bytebuffer = ByteBuffer.wrap(decompressed);

        Float[] output = new Float[numElems];
        for(int i=0; i<numElems; i++){
            output[i] = bytebuffer.getFloat();
        }

        return output;
    }

    public synchronized byte[] readDeflated(long index, long nextIndex, int numElems) throws IOException {

        if(index > ra.length()){
            throw new IOException("IO Error: SNP ID out of scope!");
        }

        long seekloc = index;

        int compressedsize = 0;

        // if this is the last element, read to the end
        if(nextIndex == -1){
            long curpos  = seekloc;
            long nextpos = ra.length();
            compressedsize = (int) (nextpos - curpos); // TODO: watch out for buffer overruns here --> size might be greater than int.maxvalue
        } else {
            long curpos = seekloc;
            long nextpos = nextIndex;
            compressedsize = (int) (nextpos - curpos); // TODO: watch out for buffer overruns here --> size might be greater than int.maxvalue
        }
        if(seekloc + compressedsize > ra.length()){
            throw new IOException("IO Error: buffer to be loaded is out of scope");
        }

        byte[] buffer = new byte[compressedsize];
        ra.seek(seekloc);
        ra.read(buffer);


//	inflater.setInput(buffer);
//        inflater.finished();
//        byte[] decompressed = new byte[numElems*4];
//        inflater.inflate(decompressed);
//
//        long actuallydecompressed = inflater.getBytesWritten();
//        if(actuallydecompressed != numElems*4){
//            throw new IOException("IO Error: uncompressed data does not correspond to the size requested\t"+actuallydecompressed+"\t"+numElems*4+"\t Index: "+index);
//        }
//
//        inflater.reset();
//
//        ByteBuffer bytebuffer = ByteBuffer.wrap(decompressed);
//        Float[] output = new Float[numElems];
//        for(int i=0; i<numElems; i++){
//            output[i] = bytebuffer.getFloat();
//        }

        return buffer;
    }

    public void close() throws IOException {
        if(out!=null){
            if(nrOfBytesInWriteBuffer > 0){
                flush();
            }
            out.close();
            System.out.println("Average compression ratio: "+((double) uncompressedSize/indexposition));
        }
        if(ra!=null){
            ra.close();
        }
    }
}
