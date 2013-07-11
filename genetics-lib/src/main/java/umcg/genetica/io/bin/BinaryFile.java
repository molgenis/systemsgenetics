/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package umcg.genetica.io.bin;

import java.io.*;

/**
 *
 * @author harm-jan
 */
public class BinaryFile {

    public static final boolean W = true;
    public static final boolean R = false;
    protected final DataOutputStream os;
    protected final DataInputStream is;
    protected final String loc;
    protected final boolean writeable;

    public BinaryFile(String loc, boolean mode) throws IOException {
        if (loc.trim().length() == 0) {
            throw new IOException("Could not find file: no file specified");
        }
        this.writeable = mode;
        this.loc = loc;

        if (writeable) {
            is = null;
            os = new DataOutputStream(new BufferedOutputStream(new FileOutputStream(loc)));
        } else {
            is = new DataInputStream(new BufferedInputStream(new FileInputStream(loc)));
            os = null;
        }
    }
    
    public BinaryFile(String loc, boolean mode, int buffersize) throws IOException {
        if (loc.trim().length() == 0) {
            throw new IOException("Could not find file: no file specified");
        }
        this.writeable = mode;
        this.loc = loc;

        if (writeable) {
            is = null;
            os = new DataOutputStream(new BufferedOutputStream(new FileOutputStream(loc), buffersize));
        } else {
            is = new DataInputStream(new BufferedInputStream(new FileInputStream(loc), buffersize));
            os = null;
        }
    }

    public void writeBytes(byte[] v) throws IOException {
        if (writeable) {
            os.write(v);
        } else {
            throw new IOException("File is read only.");
        }
    }

    public void writeInt(int v) throws IOException {
        if (writeable) {
            os.writeInt(v);
        } else {
            throw new IOException("File is read only.");
        }
    }

    public void writeString(String s) throws IOException {
        if (writeable) {
            os.writeChars(s);
        } else {
            throw new IOException("File is read only.");
        }
    }

    public void writeBool(boolean b) throws IOException {
        if (writeable) {
            os.writeBoolean(b);
        } else {
            throw new IOException("File is read only.");
        }
    }

    public void writeFloat(float f) throws IOException {
        if (writeable) {
            os.writeFloat(f);
        } else {
            throw new IOException("File is read only.");
        }
    }

    public void writeDouble(double d) throws IOException {
        if (writeable) {
            os.writeDouble(d);
        } else {
            throw new IOException("File is read only.");
        }
    }

    public void writeLong(long l) throws IOException {
        if (writeable) {
            os.writeLong(l);
        } else {
            throw new IOException("File is read only.");
        }
    }

    // read functions
    public int readInt() throws IOException, EOFException {
        if (writeable) {
            throw new IOException("File is write only.");
        } else {
            return is.readInt();
        }
    }

    public boolean readBool() throws IOException, EOFException {
        if (writeable) {
            throw new IOException("File is write only.");
        } else {
            return is.readBoolean();
        }
    }

    public String readString() throws IOException, EOFException {
        if (writeable) {
            throw new IOException("File is write only.");
        } else {
            return is.readUTF();
        }
    }

    public float readFloat() throws IOException, EOFException {
        if (writeable) {
            throw new IOException("File is write only.");
        } else {
            return is.readFloat();

        }
    }

    public double readDouble() throws IOException, EOFException {
        if (writeable) {
            throw new IOException("File is write only.");
        } else {
            return is.readDouble();
        }
    }

    public long readLong() throws IOException, EOFException {
        if (writeable) {
            throw new IOException("File is write only.");
        } else {
            return is.readLong();
        }
    }

    public void close() throws IOException {
        if (writeable) {
            os.close();
        } else {
            is.close();
        }
    }

    public void writeByte(byte b) throws IOException {
        if (writeable) {
            os.writeByte(b);
        } else {
            throw new IOException("File is read only.");
        }
    }

    public int read() throws IOException {
        return is.read();
    }

    public void write(int b) throws IOException {
        os.write(b);
    }
}
