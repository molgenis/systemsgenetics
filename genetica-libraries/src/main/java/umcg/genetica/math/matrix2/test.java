package umcg.genetica.math.matrix2;

import org.apache.commons.compress.compressors.lz4.FramedLZ4CompressorInputStream;
import org.apache.commons.compress.compressors.lz4.FramedLZ4CompressorOutputStream;
import org.apache.commons.io.input.RandomAccessFileInputStream;
import org.apache.commons.io.output.CountingOutputStream;

import java.io.*;
import java.security.DigestOutputStream;
import java.security.MessageDigest;
import java.security.NoSuchAlgorithmException;

public class test {

    public static void main(String[] args) throws IOException, NoSuchAlgorithmException {

        final MessageDigest matrixFileMd5Digest = MessageDigest.getInstance("MD5");;

        CountingOutputStream matrixFileWriter = new CountingOutputStream(new DigestOutputStream(new BufferedOutputStream(new FileOutputStream(new File("D:\\test3.datg")), 262144), matrixFileMd5Digest));

        FramedLZ4CompressorOutputStream compressionFrame = new FramedLZ4CompressorOutputStream(matrixFileWriter);

        compressionFrame.write(1);
        compressionFrame.write(2);
        compressionFrame.write(3);
        compressionFrame.write(4);

        compressionFrame.finish();

        System.out.printf("Start of block 2: %d%n", matrixFileWriter.getByteCount());

        compressionFrame = new FramedLZ4CompressorOutputStream(matrixFileWriter);

        double[] rowData = new double[]{5,Double.NaN, Double.NEGATIVE_INFINITY, Double.POSITIVE_INFINITY, 0.00000001};

        byte[] rowBuffer = new byte[8 * rowData.length];

        int b = 0;
        for (int c = 0; c < rowData.length; ++c) {
            System.out.println(rowData[c]);
            long v = Double.doubleToLongBits(rowData[c]);
            rowBuffer[b++] = (byte) (v >>> 56);
            rowBuffer[b++] = (byte) (v >>> 48);
            rowBuffer[b++] = (byte) (v >>> 40);
            rowBuffer[b++] = (byte) (v >>> 32);
            rowBuffer[b++] = (byte) (v >>> 24);
            rowBuffer[b++] = (byte) (v >>> 16);
            rowBuffer[b++] = (byte) (v >>> 8);
            rowBuffer[b++] = (byte) (v);

        }
        compressionFrame.write(rowBuffer, 0, rowBuffer.length);
        compressionFrame.finish();

        System.out.printf("End of block 2: %d%n", matrixFileWriter.getByteCount());

        matrixFileWriter.getByteCount();

        matrixFileWriter.write(42);

        matrixFileWriter.close();

        System.out.println("reading");

        RandomAccessFile matrixFileReader = new RandomAccessFile(new File("D:\\test3.datg"), "r");
        FramedLZ4CompressorInputStream reader = new FramedLZ4CompressorInputStream(new RandomAccessFileInputStream(matrixFileReader, false), true);

        byte[] blockData = new byte[4];
        int readCount = reader.read(blockData, 0, blockData.length);
        System.out.println("bytes read: " + readCount);



        System.out.println(blockData[0]);
        System.out.println(blockData[1]);
        System.out.println(blockData[2]);
        System.out.println(blockData[3]);




        matrixFileReader.seek(23);
        BufferedInputStream reader2 = new BufferedInputStream(new FramedLZ4CompressorInputStream(new RandomAccessFileInputStream(matrixFileReader, false)));
        blockData = new byte[40];
        readCount = reader2.readNBytes(blockData, 0, blockData.length);

        System.out.println("bytes read: " + readCount);

    b = 0;
        double value = Double.longBitsToDouble(((long) blockData[b++] << 56)
                + ((long) (blockData[b++] & 255) << 48)
                + ((long) (blockData[b++] & 255) << 40)
                + ((long) (blockData[b++] & 255) << 32)
                + ((long) (blockData[b++] & 255) << 24)
                + ((blockData[b++] & 255) << 16)
                + ((blockData[b++] & 255) << 8)
                + ((blockData[b++] & 255)));

        System.out.println(value);

         value = Double.longBitsToDouble(((long) blockData[b++] << 56)
                + ((long) (blockData[b++] & 255) << 48)
                + ((long) (blockData[b++] & 255) << 40)
                + ((long) (blockData[b++] & 255) << 32)
                + ((long) (blockData[b++] & 255) << 24)
                + ((blockData[b++] & 255) << 16)
                + ((blockData[b++] & 255) << 8)
                + ((blockData[b++] & 255)));

        System.out.println(value);
        value = Double.longBitsToDouble(((long) blockData[b++] << 56)
                + ((long) (blockData[b++] & 255) << 48)
                + ((long) (blockData[b++] & 255) << 40)
                + ((long) (blockData[b++] & 255) << 32)
                + ((long) (blockData[b++] & 255) << 24)
                + ((blockData[b++] & 255) << 16)
                + ((blockData[b++] & 255) << 8)
                + ((blockData[b++] & 255)));

        System.out.println(value);
        value = Double.longBitsToDouble(((long) blockData[b++] << 56)
                + ((long) (blockData[b++] & 255) << 48)
                + ((long) (blockData[b++] & 255) << 40)
                + ((long) (blockData[b++] & 255) << 32)
                + ((long) (blockData[b++] & 255) << 24)
                + ((blockData[b++] & 255) << 16)
                + ((blockData[b++] & 255) << 8)
                + ((blockData[b++] & 255)));

        System.out.println(value);
        value = Double.longBitsToDouble(((long) blockData[b++] << 56)
                + ((long) (blockData[b++] & 255) << 48)
                + ((long) (blockData[b++] & 255) << 40)
                + ((long) (blockData[b++] & 255) << 32)
                + ((long) (blockData[b++] & 255) << 24)
                + ((blockData[b++] & 255) << 16)
                + ((blockData[b++] & 255) << 8)
                + ((blockData[b++] & 255)));



        System.out.println(value);

//        byte[] bytes = reader.readNBytes(4);
//        for(byte b2 : bytes){
//            System.out.println(b2);
//        }


    }
}
