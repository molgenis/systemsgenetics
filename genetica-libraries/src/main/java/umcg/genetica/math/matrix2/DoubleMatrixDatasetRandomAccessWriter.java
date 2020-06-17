package umcg.genetica.math.matrix2;

import org.apache.commons.io.input.CountingInputStream;
import org.apache.commons.io.output.CountingOutputStream;
import umcg.genetica.console.ProgressBar;
import umcg.genetica.io.text.TextFile;

import java.io.*;
import java.nio.ByteBuffer;
import java.nio.channels.Channels;
import java.nio.channels.FileChannel;
import java.util.ArrayList;
import java.util.LinkedHashMap;

public class DoubleMatrixDatasetRandomAccessWriter {

    FileChannel channel;
    private int nrCols;
    private int nrRows;
    private int buffersize = 32 * 1024;
    private CountingOutputStream counter;
    private DataOutputStream os;
    private int bytesPerRow;
    private long headerLen;
    private long currentPos;
    private ByteBuffer bytebuffer;
    private LinkedHashMap<String, Integer> hashCols;
    private LinkedHashMap<String, Integer> hashRows;

    public void initializeFullMatrix(ArrayList<String> rows, ArrayList<String> cols, String out) throws IOException {
        initialize(rows, cols, out);
        initializeWithNaNs();
    }

    public void open(String loc) throws IOException {
        channel = new RandomAccessFile(loc + ".dat", "rw").getChannel();
        CountingInputStream countertmp = new CountingInputStream(new BufferedInputStream(Channels.newInputStream(channel), buffersize));
        DataInputStream is = new DataInputStream(countertmp);
        nrRows = is.readInt();
        nrCols = is.readInt();
        hashRows = DoubleMatrixDataset.loadIdentifiers(loc + ".rows.txt");
        hashCols = DoubleMatrixDataset.loadIdentifiers(loc + ".cols.txt");

        headerLen = 8;
        currentPos = headerLen;
        bytesPerRow = 8 * nrCols;
        buffersize = bytesPerRow * 10;
        channel.position(currentPos);
        counter = new CountingOutputStream(new BufferedOutputStream(Channels.newOutputStream(channel), buffersize));
        os = new DataOutputStream(counter);

        System.out.println("Read header. current pos: " + channel.position());
        System.out.println("Header: " + headerLen);
    }

    public void initialize(ArrayList<String> rows, ArrayList<String> cols, String out) throws IOException {
        System.out.println("Initializing new output matrix: " + out);

        channel = new RandomAccessFile(out + ".dat", "rw").getChannel();
        counter = new CountingOutputStream(new BufferedOutputStream(Channels.newOutputStream(channel), buffersize));
        os = new DataOutputStream(counter);
        os.writeInt(rows.size());
        os.writeInt(cols.size());

        TextFile tf = new TextFile(out + ".rows.txt", TextFile.W);
        for (String s : rows) {
            tf.writeln(s);
        }
        tf.close();
        tf = new TextFile(out + ".cols.txt", TextFile.W);
        for (String s : cols) {
            tf.writeln(s);
        }
        tf.close();


        headerLen = 8;
        currentPos = 8;
        this.nrCols = cols.size();
        this.nrRows = rows.size();
        System.out.println("Matrix has size: " + this.nrRows + " rows x " + this.nrCols + " cols");
        bytesPerRow = cols.size() * 8;
        buffersize = bytesPerRow * 10;
    }


    private void initializeWithNaNs() throws IOException {

        ProgressBar pb = new ProgressBar(nrRows, "Initializing with NaNs: " + nrRows + " rows x " + nrCols + " cols");
        for (int row = 0; row < nrRows; row++) {
            for (int col = 0; col < nrCols; col++) {
                os.writeDouble(Double.NaN);
            }
            pb.iterate();
        }
        pb.close();
    }


    public void writeRow(double[] cols) throws IOException {
        if (cols.length != nrCols) {
            throw new IllegalArgumentException("Length of cols not equal: " + cols.length + " found, " + nrCols + " expected.");
        }
        for (int c = 0; c < cols.length; c++) {
            os.writeDouble(cols[c]);
        }
        currentPos += bytesPerRow;
    }

    public void write(double d) throws IOException {
        os.writeDouble(d);
        currentPos += 8;
    }

    public void writeRow(int row, double[] cols) throws IOException {
        long seekLoc = ((long) row * bytesPerRow) + headerLen;

        if (seekLoc > channel.size()) {
            throw new IllegalArgumentException("Seek location for row: " + row + ", " + seekLoc + " is outside file size: " + channel.size());
        }

        // if row is the next row, just write.
        if (seekLoc - currentPos == 0) {

            writeRow(cols);
        } else {
            // else, random access to new location

            channel.position(seekLoc);
            if (bytebuffer == null) {
                bytebuffer = ByteBuffer.wrap(new byte[bytesPerRow]);

            }
            channel.write(bytebuffer);

            // this is probably extremely slow?
            counter = new CountingOutputStream(new BufferedOutputStream(Channels.newOutputStream(channel), buffersize));
            os = new DataOutputStream(counter);
            currentPos = seekLoc + bytesPerRow;


        }
    }

    ByteBuffer singledouble;

    public void write(int row, int col, double val) throws IOException {
        long seekLoc = ((long) row * bytesPerRow) + headerLen + (col * 8);
//        System.out.println(row + "\t" + col + "\t" + seekLoc + "\t" + currentPos + "\t" + val);
        if (seekLoc - currentPos == 0) {

            os.writeDouble(val);
            os.flush();
            currentPos = seekLoc + 8;

        } else {
            if (seekLoc > channel.size()) {
                throw new IllegalArgumentException("Seek location for row: " + row + ", " + seekLoc + " is outside file size: " + channel.size());
            }

            if (singledouble == null) {
                singledouble = ByteBuffer.allocate(8);
            }

            singledouble.putDouble(val);
            singledouble.flip();

//            System.out.println("Seeking: " + seekLoc);
            channel.position(seekLoc);

            channel.write(singledouble);
            currentPos = seekLoc + 8;
            singledouble.compact();


            // this is probably extremely slow?
            counter = new CountingOutputStream(new BufferedOutputStream(Channels.newOutputStream(channel), buffersize));
            os = new DataOutputStream(counter);
        }

    }

    ByteBuffer blockbuffer;

    public void writeBlock(int startRow, int startCol, double[] vals) throws IOException {
        long seekLoc = ((long) startRow * bytesPerRow) + headerLen + (startCol * 8);
//		System.out.println(startRow + "\t" + startCol + "\t" + seekLoc + "\t" + vals[0]);
        channel.position(seekLoc);

        if (blockbuffer == null || blockbuffer.limit() != vals.length * 8) {
            blockbuffer = ByteBuffer.allocate(vals.length * 8);
        }

        for (int b = 0; b < vals.length; b++) {
            blockbuffer.putDouble(vals[b]);
        }
        blockbuffer.flip();

        channel.write(blockbuffer);
        blockbuffer.compact();
    }

    public void close() throws IOException {
        os.close();
        counter.close();
        channel.close();
    }

}
